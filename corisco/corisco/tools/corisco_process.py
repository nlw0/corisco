#!/usr/bin/python
#coding:utf-8

# Copyright 2011 Nicolau Leal Werneck, Anna Helena Reali Costa and
# Universidade de São Paulo
#
# Licensed under theglim Apache License, Version 2.0 (the "License");
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

###############################################################################
## Open a (pinhole camera model) picture, find its orientation and
## extract the rectified.
##
## Changes that must be performed as soon as possible: read the
## intrinsic parameters from somewhere, and make it easy to switch to
## the equirectangular model.

import sys
import os
import cProfile
import code

# if __name__ == '__main__':
#     if sys.argv[0][-7:] == '-nop.py':
#         PlotStuff=False
#     else:
#         PlotStuff=True
#         import matplotlib
#         from pylab import *
#         from plot_aux import *
#         if sys.argv[0][-7:] == '-nox.py':
#             matplotlib.use('Agg') 



import scipy.io

from numpy import *
import numpy as np

from corisco import Picture
from corisco import corisco_aux
from corisco.quaternion import Quat, random_quaternion

import simplejson

import scipy.optimize

import Image

set_printoptions(precision=7)

import filtersqp

######################################################################
## Functions used in the optimization

## Conic constraint
def val_c(x):
    return np.linalg.norm(x)
def grad_c(x):
    return x/np.linalg.norm(x)
def hess_c(x):
    nx = np.linalg.norm(x)
    return (nx**2 * np.identity(x.shape[0]) - np.outer(x,x)) / nx**3

## Target function
def val_f(x,*fargs):
    return corisco_aux.angle_error(x, fargs[0],fargs[1],fargs[2])
def grad_f(x,*fargs):
    return corisco_aux.angle_error_gradient(x, fargs[0],fargs[1],fargs[2])
def hess_f(x,*fargs):
    return corisco_aux.angle_error_hessian(x, fargs[0],fargs[1],fargs[2])
##
######################################################################


def utime():
    '''Pick user time from os to measure the code performance.'''
    return os.times()[0]

def main():

    PlotStuff=True
    import matplotlib
    from pylab import *
    from corisco.plot_aux import *
    if sys.argv[0][-7:] == '-nox.py':
        matplotlib.use('Agg') 

    tt_total_init = utime()
    if PlotStuff:
        rc('text',usetex=False)

    ## Avoid zero divide warnins...
    np.seterr(divide='ignore')

    if PlotStuff:
        ## Plot stuff immediately
        ion()

    #################################################################
    ## Load image and initialize pic object
    tt_init = utime()

    ## Sets filename from input argument
    if len(sys.argv) < 3:
        print sys.argv[0], '<job_file.json> <frame_number>'
        raise Exception('Insufficient number of parameters')

    finput = open(sys.argv[1])
    job_params = simplejson.load(finput)
    finput.close()

    fileroot = job_params['root_directory']

    framenum = int(sys.argv[2])
    filename = fileroot+'/frames/'+job_params['filename_format']%framenum

    im = Image.open(filename)
    frame = array(im.convert('RGB'), dtype=float)
    imr = array(im.convert('RGB'), dtype=float)
    imr = imr[:,:,:3] #remove alpha channel

    # Smooth out
    if ("gaussian_smoothing_factor" in job_params.keys() and 
        job_params["gaussian_smoothing_factor"] > 0):
        for c in range(3):
            imr[:,:,c] = scipy.ndimage.gaussian_filter( imr[:,:,c], \
                          double(job_params["gaussian_smoothing_factor"]))

    ## Creates picture object
    if not 'projection_model' in job_params:
        raise "Missing camra model in job file"
    else:
        model = job_params['projection_model']

    if not (model in ['pinhole', 'harris', 'polar_equidistant', 'cylindrical_equidistant']):
        raise NotImplementedError

    elif model == 'pinhole':
        ## Intrinsic parameters
        focal_distance = job_params['focal_distance']
        p_point = array(job_params['principal_point'])
        i_param = array([0.0, focal_distance, p_point[0], p_point[1]])

    elif model == 'harris':
        ## Intrinsic parameters
        focal_distance = job_params['focal_distance']
        p_point = array(job_params['principal_point'])
        distortion = job_params['distortion_coefficient']
        i_param = array([2.0, focal_distance, p_point[0], p_point[1], distortion])
        
    elif model == 'polar_equidistant':
        ## Intrinsic parameters
        focal_distance = job_params['focal_distance']
        p_point = array(job_params['principal_point'])
        i_param = array([3.0, focal_distance, p_point[0], p_point[1]])

    elif model == 'cylindrical_equidistant':
        ## Intrinsic parameters
        focal_distance = job_params['focal_distance']
        p_point = array(job_params['principal_point'])
        i_param = array([4.0, focal_distance, p_point[0], p_point[1]])

    sqp_funcs = (val_c, grad_c, hess_c, val_f, grad_f, hess_f)
    pic = Picture(imr, i_param)


    tt_image_initialization = utime()-tt_init
    ##
    ##################################################################

    ## Edgel extractor parameters
    gmethod = job_params['edge_detection_method']
    gspc = job_params['grid_spacing']
    glim = job_params['gradient_threshold']

    do_robust = job_params['do_robust']
    initial_trials = job_params['initial_trials']

    do_optimization = job_params['do_optimization']
    optimization_tolerance = job_params['optimization_tolerance']

    if 'do_decimation' in job_params.keys():
        do_decimation = job_params['do_decimation']
        dec_t = job_params['decimator_threshold']
        dec_d = job_params['decimator_distance']
        dec_m = job_params['decimator_method']
        dec_l = job_params['decimator_lower']
    else:
        do_decimation = 0

    if 'do_multiple_initializations' in job_params.keys():
        do_multiple_initializations = job_params['do_multiple_initializations']
    else:
        do_multiple_initializations = 0

    if 'do_second_optimization' in job_params.keys():
        do_second_optimization = job_params['do_second_optimization']
    else:
        do_second_optimization = 0


    ##################################################################
    ## Extract the edgels from the image using the grid mask
    tt_init = utime()
    pic.extract_edgels(gspc, glim, method=gmethod)
    tt_edgel_extraction = utime()-tt_init
    ##
    ##################################################################

    ##################################################################
    ## Caculate the edgel (interpretation plane) normals.
    pic.calculate_edgel_normals()
    ##
    ##################################################################

    ##################################################################
    ## Filter the edgels to remove the lonely ones. Impose a kind of
    ## curvature restriction in the edgels, but still not a full
    ## countour extraction...  Uses a kd-tree, but still slow.
    tt_init = utime()
    if do_decimation and dec_t > 0:
        pic.decimate_edgels(dec_t, dec_d, dec_l, dec_m)
    tt_directional_filter = utime() - tt_init
    ##
    ##################################################################

    ##################################################################
    ## Calculate Jacobians at each edgel
    tt_init = utime()
    if do_robust == 1:
        pic.calculate_edgel_jacobians()
    tt_jacobians = utime() - tt_init
    ##
    ##################################################################

    ##################################################################
    ## Calculate initial estimate
    tt_init = utime()
    if do_robust == 1:
        fp = array(job_params["fp_robust"])
        qini,bpk_a,bpk_b,bpk_c = pic.random_search(initial_trials, fp)
        # aa = cProfile.runctx('qini,bpk_a,bpk_b,bpk_c = pic.random_search(initial_trials, fp)', globals(), locals(), "prof.txt")
        qini = qini.canonical()
    else:
        initial_trials = 0
    tt_initial_estimate = utime() - tt_init
    ##
    ##################################################################

    #################################################################
    ## Perform optimization procedure to find orientation even more
    ## accurately.
    tt_init = utime()

    ## No optimization, the estimate is just whatever came from the
    ## random search.
    if do_optimization == 0:
        qopt = Quat(qini.q)
        sqp_iters = 0
    ## Each case below runs non-linear optimization somehow.

    ## This is just a final optimization after the random search.
    elif do_robust == 1 and do_multiple_initializations == 0:
        fp = array(job_params["fp_optimization"])
        # sqp_funcs = (val_c, grad_c, hess_c, val_f, grad_f, hess_f)
        # args_f = int_parms + (pic.edgels,) + fp
        args_f = (pic.edgels, i_param, fp)

        filterSQPout = filtersqp.filterSQP(qini.q,.0,1e-3,sqp_funcs,args_f,
                                           delta_tol=optimization_tolerance)
        # filterSQPout = filtersqp.filterSQP(qini.q,0.0,0.05,sqp_funcs,args_f,
        #                                    delta_tol=optimization_tolerance)

        xo, err, sqp_iters,Llam,Lrho = filterSQPout

        qopt = Quat(xo)

    ## If there was no initial random search, the only current option
    ## is to run the optimization on multiple initial guesses, then
    ## pick up the best solution found.
    elif do_robust == 0 and do_multiple_initializations == 1:
        xini = array(job_params['initial_guesses'])   ## Read initial points form job file
        errs=[]
        Niter=0

        ## Run optimization from multiple initial points, using a less
        ## selective error function.
        fp = tuple(job_params["fp_optimization"])
        # sqp_funcs = (val_c, grad_c, hess_c, val_f, grad_f, hess_f)
        args_f = int_parms + (pic.edgels,) + fp

        for nk,xi in enumerate(xini):
            xo, err, sqp_iters = filtersqp.filterSQP(Quat(xi).q,.0,.2,sqp_funcs,args_f,
                                                     delta_tol=optimization_tolerance)
            xo = Quat(xo).canonical()
            errs.append((xo,err))
            if PlotStuff:
                print xi, xo, err

        ## Get best estimation from the multiple optimizations
        qini = (errs[argmin([ee[1] for ee in errs])][0]).canonical()

        ## This is to run a second optmization with different loss
        ## function parameters. The idea is that this second
        ## optmization would be able to use a more restrictive
        ## function, while in the first optimization we would be more
        ## concerned with global convergence.
        if do_second_optimization == 1:
            ## Optimize again using that best estimate as initial point
            fp = tuple(job_params["fp_second_optimization"])
            args_f2 = int_parms + (pic.edgels,) + fp
            xo, err, sqp_iters = filtersqp.filterSQP(qini.q,.0,.2,sqp_funcs,args_f2,
                                          delta_tol=optimization_tolerance)
            qopt = Quat(xo)
        else:
            qopt = qini

    else:
        raise BaseException('Invalid Parameters')

    tt_filtersqp = utime() - tt_init
    ##
    ##################################################################

    tt_total = utime() - tt_total_init

    if PlotStuff:
        figure(1,figsize=(6,4.5))
        title('Extracted edgels')
        imshow(.25+pic.frame*.75/255.)
        pic.plot_edgels(gca(),gspc*0.35)
        axis([0,pic.Iwidth,pic.Iheight,0])
        # savefig('edg.pdf',dpi=300)

        figure(11,figsize=(6,4.5))
        title(u'Estimated orientation')
        imshow(.2+pic.frame*.8/255.)
        #pic.plot_vdirs(gca(), 63, qopt)
        pic.plot_vdirs(gca(), 47, qopt)
        axis([0,pic.Iwidth, pic.Iheight, 0])
        # savefig('ori.pdf',dpi=300)

    qini = qini.canonical()
    qopt = qopt.canonical()

    print framenum,
    print 4*' % 10.8f'%(qini.q[0], qini.q[1], qini.q[2], qini.q[3]),
    print 4*' % 10.8f'%(qopt.q[0], qopt.q[1], qopt.q[2], qopt.q[3]),
    print pic.edgels.shape[0], initial_trials, sqp_iters,
    print 5*' %.2f'%(tt_edgel_extraction, tt_directional_filter,
                     tt_initial_estimate, tt_filtersqp, tt_total)


    code.interact()
