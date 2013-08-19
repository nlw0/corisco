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

from corisco import Picture
from corisco import corisco_aux
from corisco.quaternion import Quat, random_quaternion
import filtersqp
from hashlib import md5
import numpy as np
from numpy import array, seterr
from os import times

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
    return times()[0]

def estimate_orientation(process_args, image_file_str):
    seterr(divide='ignore')  ## Avoid zero divide warnings

    tt_total = utime()

    ## Edgel extraction parameters
    gmethod = process_args['grid']['method']
    gspc = process_args['grid']['size']
    glim = process_args['grid']['lim']
    ## Optimization parameters
    ransac_iterations = process_args['ransac_itr']
    optimization_tolerance = process_args['op_tol']
    ## Parameters from the error function. Default is Tukey bisquare with s=0.15
    error_function_parameters = array(process_args['err_func'])
    intrinsic_parameters = process_args['int_parm']

    #################################################################
    ## Load image and initialize pic object
    tt_init = utime()

    ## Creates picture object
    pic = Picture(intrinsic_parameters, image_file_str)

    if process_args['smooth'] is not None:
        pic.smooth(process_args['smooth'])

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


    filterSQPout = filtersqp.filterSQP(
        qini.q, .0, 1e-3, sqp_funcs, args_f, delta_tol=optimization_tolerance
        )

    # try:
    #     filterSQPout = filtersqp.filterSQP(
    #         qini.q, .0, 1e-3, sqp_funcs, args_f, delta_tol=optimization_tolerance
    #         )
    # except:
    #     print '*** Numerical error for input:'
    #     print {'img_md5': get_hash(image_file_str),
    #            'proc_args': process_args}
    #     raise SystemExit

    xo, err, sqp_iters,Llam,Lrho = filterSQPout
    qopt = Quat(xo)

    tt_filtersqp = utime() - tt_filtersqp

    tt_total = utime() - tt_total

    first_orientation_estimate = qini.canonical().q.tolist()
    final_orientation_estimate = qopt.canonical().q.tolist()

    output_data = {
        'input': {'img_md5': get_hash(image_file_str),
                  'proc_args': process_args},
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

    return output_data, pic

def get_hash(image_file_str):
    return md5(image_file_str).hexdigest()
