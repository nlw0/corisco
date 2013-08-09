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

import argparse

import sys
import os
import cProfile
import code

import json

from numpy import array
import numpy as np

from corisco import Picture
from corisco import corisco_aux
from corisco.quaternion import Quat, random_quaternion

import simplejson

import scipy.optimize

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
    ## Avoid zero divide warnings
    np.seterr(divide='ignore')

    ## Command-line argument parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', type=str)
    parser.add_argument('frame_number', type=int)
    parser.add_argument('--do_plot', default=False, action='store_true')

    parser.add_argument('--edge_detection_method', type=str, default='shigeru')
    parser.add_argument('--grid_spacing', type=int, default=4)
    parser.add_argument('--gradient_threshold', type=float, default=20.0)
    parser.add_argument('--ransac_iterations', type=int, default=1000)
    parser.add_argument('--optimization_tolerance', type=float, default=1e-12)
    parser.add_argument('--smooth', type=float)
    parser.add_argument('--error_function_parameters', type=float, nargs=3, default=[3.0, 1.0, 0.15])

    parser.add_argument('--directions_plot_spacing', type=int, default=47)
    parser.add_argument('--save_graphics', default=False, action='store_true')


    args = parser.parse_args()

    tt_total = utime()

    ## Edgel extraction parameters
    gmethod = args.edge_detection_method
    gspc = args.grid_spacing
    glim = args.gradient_threshold
    ## Optimization parameters
    ransac_iterations = args.ransac_iterations
    optimization_tolerance = args.optimization_tolerance
    ## Parameters from the error function. Default is Tukey bisquare with s=0.15
    error_function_parameters = array(args.error_function_parameters)

    #################################################################
    ## Load image and initialize pic object
    tt_init = utime()

    with open(args.input_file) as finput:
        job_params = simplejson.load(finput)

    ## Creates picture object
    pic = Picture(job_params, args.frame_number)

    if args.smooth is not None:
        pic.smooth(args.smooth)

    ## Extract the edgels from the image using the grid mask
    tt_edgel_extraction = utime()
    pic.extract_edgels(gspc, glim, method=gmethod)
    tt_edgel_extraction = utime() - tt_edgel_extraction

    ## Caculate the edgel normals (interpretation plane), and Jacobians.
    pic.calculate_edgel_normals()
    pic.calculate_edgel_jacobians()

    ## Calculate initial estimate
    tt_initial_estimate = utime()
    qini = pic.random_search(ransac_iterations, error_function_parameters)
    qini = qini.canonical()
    tt_initial_estimate = utime() - tt_initial_estimate

    ## Perform second-stage continuous optimization procedure, based on FilterSQP.
    tt_filtersqp = utime()

    sqp_funcs = (val_c, grad_c, hess_c, val_f, grad_f, hess_f)
    args_f = (pic.edgels, pic.i_param, error_function_parameters)
    try:
        filterSQPout = filtersqp.filterSQP(
            qini.q, .0, 1e-3, sqp_funcs, args_f, delta_tol=optimization_tolerance
            )
    except:
        print '*** Numerical error for input:'
        print {'set': args.input_file,
               'frame': args.frame_number,
               'grid': {'size': gspc,
                        'lim': glim,
                        'method': gmethod
                        },
               'ransac_itr': ransac_iterations,
               'smooth': args.smooth,
               'optol': optimization_tolerance
               }
        raise SystemExit

    xo, err, sqp_iters,Llam,Lrho = filterSQPout
    qopt = Quat(xo)

    tt_filtersqp = utime() - tt_filtersqp

    tt_total = utime() - tt_total

    first_orientation_estimate = qini.canonical().q.tolist()
    final_orientation_estimate = qopt.canonical().q.tolist()

    output_data = {
        'input': {'set': args.input_file,
                  'frame': args.frame_number,
                  'grid': {'size': gspc,
                           'lim': glim,
                           'method': gmethod
                           },
                  'ransac_itr': ransac_iterations,
                  'smooth': args.smooth,
                  'optol': optimization_tolerance
                  },
        'time': {
            'total': tt_total,
            'edgel': tt_edgel_extraction,
            'ransac': tt_initial_estimate,
            'sqp': tt_filtersqp
            },
        'Nedgels': pic.edgels.shape[0],
        'sqp_itr': sqp_iters,
        'ransac_ori_est': first_orientation_estimate,
        'ori_est': final_orientation_estimate,
        }

    print json.dumps(output_data)

    if args.do_plot or args.save_graphics:
        import matplotlib            
        if args.save_graphics:
            matplotlib.use('Agg') 

        import pylab as plt

        if not args.save_graphics:
            plt.ion()
        plt.figure(1,figsize=(6,4.5))
        plt.title('Extracted edgels')
        plt.imshow(.25+pic.frame*.75/255.)
        pic.plot_edgels(plt.gca(), gspc * 0.35)
        plt.axis([0, pic.Iwidth, pic.Iheight, 0])
        if args.save_graphics:
            plt.savefig('edg.pdf',dpi=300)

        plt.figure(11,figsize=(6,4.5))
        plt.title(u'Estimated orientation')
        plt.imshow(.2 + pic.frame * .8 / 255.)
        pic.plot_vdirs(plt.gca(), args.directions_plot_spacing, qopt)
        plt.axis([0, pic.Iwidth, pic.Iheight, 0])
        
        if args.save_graphics:
            plt.savefig('ori.pdf',dpi=300)

        if not args.save_graphics:
            code.interact()
