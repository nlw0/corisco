#!/usr/bin/python
#coding:utf-8

# Copyright 2012 Nicolau Leal Werneck, Anna Helena Reali Costa and
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

import sys
import numpy as np
import scipy
import scipy.ndimage

import corisco_aux

from quaternion import Quat

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

class Picture:
  '''For perspective transform, "pinhole" cameras.'''

  def __init__(self, frame, i_param):

    try:
      # self.frame = scipy.array( pylab.flipud(pylab.imread(filename)), dtype=float)
      self.frame = frame
    except IOError:
      print "SOCORRO!, File not found"
      print filename
      self=None
      raise Exception('File not found')

    self.grid_inc = 1 ## The increment when using a grid...

    ## Image dimensions
    self.Iheight, self.Iwidth = self.frame.shape[:2]

    self.edgels = []
    self.new_labels = []
    self.rect_edgels = []

    self.i_param = i_param
    
    ## Focal distance (default)
    self.fd = i_param[1]

    ## Principal point
    self.middlex = i_param[2]
    self.middley = i_param[3]



  #############################################################################
  ## Call the edgel extraction procedure.
  def extract_edgels(self, gstep, glim, method=0):
    if method == 0:
      self.derivative_filter = self.sobel_filter
    elif method == 1:
      self.derivative_filter = self.scharr_filter
    elif method == 2:
      self.derivative_filter = self.shigeru_filter
    elif method == 3:
      self.derivative_filter = self.zernike_V11_5x5
      self.laplacean_filter = self.zernike_V20_5x5
    elif method == 4:
      self.derivative_filter = self.zernike_V11_7x7
      self.laplacean_filter = self.zernike_V20_7x7
    elif method == 5:
      self.derivative_filter = self.zernike_V11_9x9
      self.laplacean_filter = self.zernike_V20_9x9
    else:
      raise 'Invalid edgel extraction method.'
      
      
    ## Calculate gradients
    self.gradx = scipy.zeros(self.frame.shape, dtype=np.float32)
    self.grady = scipy.zeros(self.frame.shape, dtype=np.float32)
    for c in range(3):
      scipy.ndimage.convolve(self.frame[:,:,c], self.derivative_filter,  self.gradx[:,:,c] )
      scipy.ndimage.convolve(self.frame[:,:,c], self.derivative_filter.T, self.grady[:,:,c] )

    if method >= 3 and method <= 5:
      ## Calculate laplacean
      self.A20 = scipy.zeros(self.frame.shape, dtype=np.float32)
      for c in range(3):
        scipy.ndimage.convolve(self.frame[:,:,c], self.laplacean_filter, self.A20[:,:,c] )

    ## Run edgel extraction procedure using the calculated gradients
    if method >= 0 and method < 3:
      self.edgels = corisco_aux.edgel_extractor(gstep, glim, self.gradx, self.grady)
    elif method >= 3 and method <= 5:
      self.edgels = corisco_aux.edgel_extractor_zernike(gstep, glim, self.gradx,
                                                   self.grady, self.A20)
    else:
      raise 'Invalid edgel extraction method.'

    ## Store number of edgels
    self.Ned = self.edgels.shape[0]

    ## Array that contains the label of each edgel
    self.labels = np.ascontiguousarray( np.zeros(len(self.edgels), dtype=np.int32) )
  ##
  #############################################################################   

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
  def random_search(self, initial_trials, fp):
    ## Default M-estimator uses Tukey with 0.15 error variance.

    bestv = np.Inf ## Smallest value found
    for k in range(initial_trials):
      ## Pick indices of the reference normals. Re-sample until we
      ## get a list of three different values.
      pk_a = np.random.random_integers(0,self.Ned-1)
      pk_b = np.random.random_integers(0,self.Ned-1)
      while pk_b == pk_a:
          pk_b = np.random.random_integers(0,self.Ned-1)
      pk_c = np.random.random_integers(0,self.Ned-1)
      while pk_c == pk_a or pk_c == pk_b:
          pk_c = np.random.random_integers(0,self.Ned-1)

      ## Get the normals with the first two chosen indices, and
      ## calculate a rotation matrix that has the x axis aligned to
      ## them.
      n_a = self.normals[pk_a]
      n_b = self.normals[pk_b]
      vp1 = np.cross(n_a, n_b)
      vp1 = vp1 * (vp1**2).sum()**-0.5
      q1 = aligned_quaternion(vp1)
      # q1 = aligned_quaternion_from_dirs(n_a, n_b)

      ## Pick a new random third norm, and find the rotation to align
      ## the y direction to this edgel.
      n_c = self.normals[pk_c]
      vaux = np.dot(n_c, q1.rot())
      ang = np.arctan2(vaux[1], -vaux[2])
      q2 = Quat(np.sin(ang/2),0,0) * q1 ## The resulting orientation

      ## Find the value of the target function for this sampled
      ## orientation.
      # newv = camori_angle_error(
      #   q2.q, self.edgels, self.i_param, fp)
      newv = corisco_aux.angle_error_with_jacobians(
       q2.q, self.edgels, self.jacobians, fp)

      ## If the value is the best yet, store solution.
      if newv <= bestv :
        bestv = newv
        bpk_a = pk_a
        bpk_b = pk_b
        bpk_c = pk_c
        qopt = q2
    return qopt, bpk_a, bpk_b, bpk_c
  ##
  #############################################################################

  #############################################################################
  ## Plotting methods
  ##
  def plot_edgels(self, ax, scale=1):
    # ax.plot(self.edgels[:,0], self.edgels[:,1], 'r+', mew=1.0)

    # ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
    #         (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
    #         'r-')

    # ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
    #         (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
    #         '-',lw=3.0,alpha=0.25,color='#ff0000')
    ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
            (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
            '-',lw=1.0,alpha=1.0,color='#ff0000')
    
    # ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
    #         (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
    #         '-',lw=3.0,alpha=0.4,color='#aaccff')
    #         # '-',lw=3.0,alpha=0.7,color='#000044')
    # ax.plot((self.edgels[:,[0,0]] - scale*np.c_[-self.edgels[:,3], self.edgels[:,3]]).T,
    #         (self.edgels[:,[1,1]] + scale*np.c_[-self.edgels[:,2], self.edgels[:,2]]).T,
    #         '-',lw=1.0,alpha=1.0,color='#ff8844')
    #         # '-',lw=1.0,alpha=1.0,color='#ff8844')

  def plot_edgels_lab(self, ax, scale=1):
    for lab in range(4):
      sel = np.nonzero(self.new_labels==lab)[0]
      ed = self.edgels[sel]
      ax.plot((ed[:,[0,0]] - scale*np.c_[-ed[:,3], ed[:,3]]).T,
              (ed[:,[1,1]] + scale*np.c_[-ed[:,2], ed[:,2]]).T, dir_colors[lab])

  def plot_vdirs(self, ax, spacing, orientation, labrange=[0,1,2]):
    '''
    Plot the vanishing point directions at various pixels. ax is a
    matplotlib axes, taken with "gca()". Spacing the separation
    between the points, and myR the rotation matrix.

    '''

    qq = spacing*0.45*np.array([-1,+1])
    bx = 0.+(self.Iwidth/2)%spacing
    by = 0.+(self.Iheight/2)%spacing
    qL = np.mgrid[bx:self.Iwidth:spacing,by:self.Iheight:spacing].T.reshape((-1,2))
    Nq = qL.shape[0]
    # vL = equidistant_vdirs(self.middlex, self.middley, self.fd,
    #                             orientation.q, np.array(qL, dtype=np.float32))
    vL = corisco_aux.calculate_vdirs(orientation.q, np.array(qL, dtype=np.float32), self.i_param)
    LL = np.zeros((3,Nq,4))
    for lab in labrange:
      for num in range(Nq):
        vx,vy = vL[lab,num]
        k,j = qL[num]
        LL[lab,num,:] = np.r_[k+vx*qq, j+vy*qq]
    for lab in range(3):
      ax.plot( LL[lab,:,:2].T, LL[lab,:,2:].T, dir_colors[lab], lw=3)
  ##
  #############################################################################

  def calculate_edgel_normals(self):
    self.normals = corisco_aux.calculate_normals(self.edgels, self.i_param)

  def calculate_edgel_jacobians(self):
    self.jacobians = corisco_aux.calculate_all_jacobians(self.edgels, self.i_param)

  ## Directional derivative filters
  sobel_filter=scipy.array([
      [-1, 0, 1],
      [-2, 0, 2],
      [-1, 0, 1] ])/8.0

  scharr_filter=scipy.array([
      [ -3, 0, 3],
      [-10, 0,10],
      [ -3, 0, 3] ])/32.0

  ## Directional derivative filter
  shigeru_filter=-scipy.array([
      [ -0.003776, -0.010199, 0., 0.010199, 0.003776 ],
      [ -0.026786, -0.070844, 0., 0.070844, 0.026786 ],
      [ -0.046548, -0.122572, 0., 0.122572, 0.046548 ],
      [ -0.026786, -0.070844, 0., 0.070844, 0.026786 ],
      [ -0.003776, -0.010199, 0., 0.010199, 0.003776 ]
      ])

  zernike_V11_5x5 = np.array([
      [ -146.67,  -468.68,     0.  ,   468.68,   146.67],
      [ -933.33,  -640.  ,     0.  ,   640.  ,   933.33],
      [-1253.33,  -640.  ,     0.  ,   640.  ,  1253.33],
      [ -933.33,  -640.  ,     0.  ,   640.  ,   933.33],
      [ -146.67,  -468.68,     0.  ,   468.68,   146.67]])/19368.05


  zernike_V11_7x7 = np.array([
      [   0.  , -150.46, -190.2 ,    0.  ,  190.2 ,  150.46,    0.  ],
      [-224.25, -465.74, -233.24,    0.  ,  233.24,  465.74,  224.25],
      [-573.37, -466.47, -233.24,    0.  ,  233.24,  466.47,  573.37],
      [-689.99, -466.47, -233.24,    0.  ,  233.24,  466.47,  689.99],
      [-573.37, -466.47, -233.24,    0.  ,  233.24,  466.47,  573.37],
      [-224.25, -465.74, -233.24,    0.  ,  233.24,  465.74,  224.25],
      [   0.  , -150.46, -190.2 ,    0.  ,  190.2 ,  150.46,    0.  ]])/27314.0


  zernike_V11_9x9 = np.array([[ 0., -11.73, -109.17, -94.2, 0., 94.2, 109.17,
                             11.73, 0. ],
                           [ -16.09, -253.68, -219.48, -109.74, 0., 109.74, 219.48,
                              253.68, 16.09],
                           [-214.91, -329.22, -219.48, -109.74, 0., 109.74, 219.48,
                             329.22, 214.91],
                           [-379.52, -329.22, -219.48, -109.74, 0., 109.74, 219.48,
                             329.22, 379.52],
                           [-434.39, -329.22, -219.48, -109.74, 0., 109.74, 219.48,
                             329.22, 434.39],
                           [-379.52, -329.22, -219.48, -109.74, 0., 109.74, 219.48,
                             329.22, 379.52],
                           [-214.91, -329.22, -219.48, -109.74, 0., 109.74, 219.48,
                             329.22, 214.91],
                           [ -16.09, -253.68, -219.48, -109.74, 0., 109.74, 219.48,
                              253.68, 16.09],
                           [ 0., -11.73, -109.17, -94.2, 0., 94.2, 109.17,
                             11.73, 0. ]]) / 35236.70736

  zernike_V20_5x5 = 2.5 * np.array([
      [  176.  ,   595.07,   505.86,   595.07,   176.  ],
      [  595.07,  -490.67, -1002.67,  -490.67,   595.07],
      [  505.86, -1002.67, -1514.67, -1002.67,   505.86],
      [  595.07,  -490.67, -1002.67,  -490.67,   595.07],
      [  176.  ,   595.07,   505.86,   595.07,   176.  ]])/19368.05


  zernike_V20_7x7 = 3.5 * np.array([
      [   0.  ,  224.66,  393.73,  395.52,  393.73,  224.66,    0.  ],
      [ 224.66,  271.06, -127.72, -261.  , -127.72,  271.06,  224.66],
      [ 393.73, -127.72, -527.56, -660.84, -527.56, -127.72,  393.73],
      [ 395.52, -261.  , -660.84, -794.11, -660.84, -261.  ,  395.52],
      [ 393.73, -127.72, -527.56, -660.84, -527.56, -127.72,  393.73],
      [ 224.66,  271.06, -127.72, -261.  , -127.72,  271.06,  224.66],
      [   0.  ,  224.66,  393.73,  395.52,  393.73,  224.66,    0.  ]])/27314.0

  zernike_V20_9x9 = 4.5 * np.array([[ 0. , 19.03, 200.93, 278.78, 290.06, 278.78, 200.93,
                             19.03, 0. ],
                           [ 19.03, 274.72, 148.35, 2.03, -46.74, 2.03, 148.35,
                             274.72, 19.03],
                           [ 200.93, 148.35, -95.51, -241.83, -290.61, -241.83, -95.51,
                             148.35, 200.93],
                           [ 278.78, 2.03, -241.83, -388.15, -436.93, -388.15, -241.83,
                             2.03, 278.78],
                           [ 290.06, -46.74, -290.61, -436.93, -485.7 , -436.93, -290.61,
                             -46.74, 290.06],
                           [ 278.78, 2.03, -241.83, -388.15, -436.93, -388.15, -241.83,
                             2.03, 278.78],
                           [ 200.93, 148.35, -95.51, -241.83, -290.61, -241.83, -95.51,
                             148.35, 200.93],
                           [ 19.03, 274.72, 148.35, 2.03, -46.74, 2.03, 148.35,
                             274.72, 19.03],
                           [ 0. , 19.03, 200.93, 278.78, 290.06, 278.78, 200.93,
                             19.03, 0. ]])/ 35236.70736

## Local variables:
## python-indent: 2
## end:
