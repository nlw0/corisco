#!/usr/bin/python
#coding:utf-8

# Copyright 2012, 2013 Nicolau Leal Werneck, Anna Helena Reali Costa and
# Universidade de São Paulo
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import code
import cStringIO

import numpy as np
from numpy import array, zeros
import scipy.ndimage

import corisco_aux

from quaternion import Quat
from filters import DERIVATIVE_FILTERS

import Image

## Colors used for plotting different labels (directions)
dir_colors=['#ea6949', '#51c373', '#a370ff', '#444444']

def aligned_quaternion(v):
  ## The largest component
  ll = np.argmax(np.abs(v))
  ee = np.zeros(3)
  ee[ll] = np.sign(v[ll])
  q = Quat(np.cross(v, ee))
  R = q.sqrt()

  if np.sign(v[ll]) > 0:
    if ll == 1:
      R = R*Quat(1,-1,-1,-1).normalize()
    if ll == 2:
      R = R*Quat(1,1,1,1).normalize()

  if np.sign(v[ll]) < 0:
    if ll == 0:
      R = R*Quat(0,0,1,0).normalize()
    if ll == 1:
      R = R*Quat(1,1,-1,1).normalize()
    if ll == 2:
      R = R*Quat(1,-1,-1,1).normalize()
  return R.inverse()

class Picture(object):

  def __init__(self, ipconf, image_file_str):

    self.read_image(image_file_str)

    ## Image dimensions
    self.Iheight, self.Iwidth = self.frame.shape[:2]

    self.edgels = []
    self.new_labels = []
    self.rect_edgels = []

    self.read_intrinsic_parameters(ipconf)

    ## Focal distance (default)
    self.fd = self.i_param[1]

    ## Principal point
    self.middlex = self.i_param[2]
    self.middley = self.i_param[3]

  def read_image(self, image_file_str):
    '''Read image, as a string, and convert to RGB floating point, no alpha.'''
    image_file = cStringIO.StringIO(image_file_str)
    self.frame = array(Image.open(image_file).convert('RGB'), dtype=float)[:,:,:3]

  def smooth(self, smooth_factor):
    '''Apply Gaussian smoothing to image.'''
    if smooth_factor > 0.0:
      for c in range(3):
        self.frame[:,:,c] = scipy.ndimage.gaussian_filter(self.frame[:,:,c], smooth_factor)

  def read_intrinsic_parameters(self, ipconf):
    if not 'projection_model' in ipconf:
      raise Exception("Camera model not specified.")
    else:
      model = ipconf['projection_model']

    if not (model in ['pinhole', 'harris', 'polar_equidistant', 'cylindrical_equidistant']):
      raise NotImplementedError

    elif model == 'pinhole':
      focal_distance = ipconf['focal_distance']
      pp = ipconf['principal_point']
      self.i_param = array([0.0, focal_distance, pp[0], pp[1]])
    elif model == 'harris':
      focal_distance = ipconf['focal_distance']
      pp = ipconf['principal_point']
      distortion = ipconf['distortion_coefficient']
      self.i_param = array([2.0, focal_distance, pp[0], pp[1], distortion])
    elif model == 'polar_equidistant':
      focal_distance = ipconf['focal_distance']
      pp = ipconf['principal_point']
      self.i_param = array([3.0, focal_distance, pp[0], pp[1]])
    elif model == 'cylindrical_equidistant':
      focal_distance = ipconf['focal_distance']
      pp = ipconf['principal_point']
      self.i_param = array([4.0, focal_distance, pp[0], pp[1]])

  def extract_edgels(self, gstep, glim, method):
    '''
    Extracts edgels from the image.

    Inputs:
      gstep - The grid spacing, or distance between liens and columns form the grid.
      glim - The minimum threshold value for the edge detector.
      method - selects between different gradient filters and edge detection methods:
        0. Sobel filter
        1. Scharr filter
        2. Shigeru filter
        3. Zernike moments, 5x5 window
        4. Zernike moments, 7x7 window
        5. Zernike moments, 9x9 window

    Outcomes:
      The extracted edgels are stored in self.edgels. self.labels is also initialized.

    '''

    try:
      self.derivative_filter, self.laplacean_filter = DERIVATIVE_FILTERS[method]
    except IndexError:
      raise Exception('Invalid edgel extraction method.')
      
    ## Calculate gradients
    gradx = zeros(self.frame.shape, dtype=np.float32)
    grady = zeros(self.frame.shape, dtype=np.float32)
    for c in range(3):
      scipy.ndimage.convolve(self.frame[:,:,c], self.derivative_filter,  gradx[:,:,c] )
      scipy.ndimage.convolve(self.frame[:,:,c], self.derivative_filter.T, grady[:,:,c] )

    if method >= 3 and method <= 5:
      ## Calculate laplacean
      A20 = zeros(self.frame.shape, dtype=np.float32)
      for c in range(3):
        scipy.ndimage.convolve(self.frame[:,:,c], self.laplacean_filter, A20[:,:,c] )

    ## Run edgel extraction procedure using the calculated gradients
    if method[:7] == 'zernike':
      self.edgels = corisco_aux.edgel_extractor_zernike(gstep, glim, gradx, grady, A20)
    else:
      self.edgels = corisco_aux.edgel_extractor(gstep, glim, gradx, grady)

    ## Array that contains the label of each edgel
    self.labels = np.ascontiguousarray( zeros(len(self.edgels), dtype=np.int32) )

  ##################################################################
  ## Calculate initial estimate by picking random edgels and
  ## calculating an orientation from them, in a direct,
  ## non-iterative way. This is a kind of random search in a space
  ## that depends on the observations, so it should be better than
  ## sweeping the whole space in a ergular way, like a complete
  ## idiot. This "adaptive" sampling of the parameter space is
  ## pretty much the same that is used in RANSAC or un
  ## J-linkage. But it's the same final target function that we are
  ## calculating at each point, not something else.
  def random_search(self, initial_trials, error_function_parameters):
    last_edgel = len(self.edgels) - 1

    bestv = np.Inf ## Smallest value found
    for k in xrange(initial_trials):
      ## Pick indices of the reference normals. Re-sample until we
      ## get a list of three different values.
      pk_a = np.random.random_integers(0, last_edgel)
      pk_b = np.random.random_integers(0, last_edgel)
      while pk_b == pk_a:
          pk_b = np.random.random_integers(0, last_edgel)
      pk_c = np.random.random_integers(0, last_edgel)
      while pk_c == pk_a or pk_c == pk_b:
          pk_c = np.random.random_integers(0, last_edgel)

      ## Get the normals with the first two chosen indices, and
      ## calculate a rotation matrix that has the x axis aligned to
      ## them.
      n_a = self.normals[pk_a]
      n_b = self.normals[pk_b]
      vp1 = np.cross(n_a, n_b)
      vp1 = vp1 * (vp1**2).sum()**-0.5
      q1 = aligned_quaternion(vp1)

      ## Pick a new random third norm, and find the rotation to align
      ## the y direction to this edgel.
      n_c = self.normals[pk_c]
      vaux = np.dot(n_c, q1.rot())
      ang = np.arctan2(vaux[1], -vaux[2])
      q2 = Quat(np.sin(ang/2),0,0) * q1 ## The resulting orientation

      ## Find the value of the target function for this sampled
      ## orientation.
      # newv = camori_angle_error(
      #   q2.q, self.edgels, self.i_param, error_function_parameters)
      newv = corisco_aux.angle_error_with_jacobians(
       q2.q, self.edgels, self.jacobians, error_function_parameters)

      ## If the value is the best yet, store solution.
      if newv <= bestv:
        bestv = newv
        bpk_a = pk_a
        bpk_b = pk_b
        bpk_c = pk_c
        qopt = q2
    return qopt

  def calculate_edgel_normals(self):
    self.normals = corisco_aux.calculate_normals(self.edgels, self.i_param)

  def calculate_edgel_jacobians(self):
    self.jacobians = corisco_aux.calculate_all_jacobians(self.edgels, self.i_param)

  def plot_edgels(self, ax, scale=1):
    ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
            (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
            '-',lw=1.0,alpha=1.0,color='#ff0000')

  def plot_edgels_lab(self, ax, scale=1):
    for lab in range(4):
      sel = np.nonzero(self.new_labels==lab)[0]
      ed = self.edgels[sel]
      ax.plot((ed[:,[0,0]] - scale*np.c_[-ed[:,3], ed[:,3]]).T,
              (ed[:,[1,1]] + scale*np.c_[-ed[:,2], ed[:,2]]).T, dir_colors[lab])

  def plot_vdirs(self, ax, spacing, Lorientation, labrange=[0,1,2]):
    '''
    Plot the vanishing point directions at various pixels. ax is a
    matplotlib axes, taken with "gca()". Spacing the separation
    between the points, and myR the rotation matrix.

    '''

    orientation = Quat(array(Lorientation))

    qq = spacing * 0.45 * array([-1, +1])
    bx = 0.+(self.Iwidth/2)%spacing
    by = 0.+(self.Iheight/2)%spacing
    qL = np.mgrid[bx:self.Iwidth:spacing,by:self.Iheight:spacing].T.reshape((-1,2))
    Nq = qL.shape[0]
    # vL = equidistant_vdirs(self.middlex, self.middley, self.fd,
    #                             orientation.q, array(qL, dtype=np.float32))
    vL = corisco_aux.calculate_vdirs(orientation.q, array(qL, dtype=np.float32), self.i_param)
    LL = zeros((3,Nq,4))
    for lab in labrange:
      for num in range(Nq):
        vx,vy = vL[lab,num]
        k,j = qL[num]
        LL[lab,num,:] = np.r_[k+vx*qq, j+vy*qq]
    for lab in range(3):
      ax.plot( LL[lab,:,:2].T, LL[lab,:,2:].T, dir_colors[lab], lw=3)

## Local variables:
## python-indent: 2
## end:
