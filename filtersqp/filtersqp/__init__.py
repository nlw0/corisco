#!/usr/bin/python
#coding: utf-8

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

## This module contains the procedure filterSQP, that implements the
## nonlinear programming technique of the same name, developed by
## Fletcher and Leyffer
## (http://www.springerlink.com/content/qqj37x00y79ygdl8/)
##
## The current implementation only supports a single constraint.


import numpy as np
from numpy import array, copy, diag, dot, nonzero, r_, c_, zeros
from numpy.linalg import norm

from filtersqp.trust_rootfind import find_step_size, calculate_distance
from fletcherfilter import FletcherFilter

def svd_eig(M):
    su, ss, sv = np.linalg.svd(dot(M.T, M))
    lam = diag(dot(dot(sv.T, M), sv.T))
    return lam, sv

def trust_region_step(G, g, rho, rho_tol):
    Nv = g.shape #number of dimensions

    ## Eigen decomposition of the matrix.
    lam, v = np.linalg.eig(G)
    ## lam, v = svd_eig(G) ## This is WRONG! (?)
    # print '==lam',lam
    # print '==v',v
    ln = lam.min()
    lni = nonzero(lam==ln)[0]
    lok = nonzero(1-(lam==ln))[0]

    ## Gradient projectin on the eigenvectors.
    alpha = dot(v.T, g)
    # print 'New QP problem:'
    # print 'lam:', lam
    # print 'rho:', rho
    # print lni, lok
    # print 'alpha:', alpha

    ## A plane. Just use the gradient...
    if (np.abs(lam)<1e-10).all():
        # print 'case 0'
        if norm(g) > 0:
            return -rho * g / norm(g)
        else:
            #print "joga pro santo"
            return rho * v[0]

    ## Test for case i -> Newton step (nu=0)
    if ln > 0:
        # print 'case 1'
        ## Matrix is positive definite.
        ## Test norm for the nu=0 case.
        beta = -alpha / lam
        nd = norm(beta)
        if nd <= rho + rho_tol:
        ## If norm is smaller than radius, this is "case i". The
        ## Newton step is feasible, so we just return it.
            delta = dot(v,beta)
            return delta

    if ln <= 0:
        ## This may be case iii or ii. Test alpha.
        if (np.abs(alpha[lni]) < 1e-10).all():
            # print "case 3"
            ## This is case iii. First we test if the saddle point is
            ## outside the trust region, in which case we perofrm the
            ## normal calculation in the end. If not, we calculate the
            ## saddle point position, and return a vector stemming
            ## from it up to the trust region boundary.
            beta = zeros(Nv)
            if lok.shape[0]>0:
                beta[lok] = -alpha[lok] / (lam[lok] - ln)
            nd = norm(beta) ## this is h_bar from Fletcher's book, pg. 104
            if nd <= rho + rho_tol:
                # print "saddle in"
                ## If norm is smaller than radius, this is "case iii"
                ## with hk>=h_bar. We return the newton step plus a
                ## vector in the direction of the smallest eigenvalue.
                db = dot(v,beta) ## delta_bar from Fletcher's book
                mv = v[:,lni[0]]
                # print db, mv
                ## Solve to find the sum that touches the region border.
                # print '----', rho, nd
                delta = db + mv * np.sqrt((rho + rho_tol) ** 2 - nd ** 2)
                # print "Return some crazy vector", delta
                return delta
        #     print "saddle out", nd, rho
        # print 'case 2'

    ##################################################################
    ## Perform an iterative root-finding procedure (secant method) to
    ## find nu.
    nu_opt = find_step_size(alpha, lam, rho, rho_tol)

    beta = -alpha / (lam + nu_opt)
    delta = dot(v,beta)
    return delta

def SQP_step(x, lam, rho, filt,
             cons_val, cons_grad, cons_hess,
             func_val, func_grad, func_hess, func_args):
    '''One step from the trust-region LM SQP algorithm.'''

    ## Calculate the local gradient and Hessian matrix
    g = func_grad(x,*func_args)
    G = func_hess(x,*func_args)
    vcon = cons_val(x)
    gcon = cons_grad(x)
    Gcon = cons_hess(x)

    ## The coefficients from the constraint function linearized around
    ## x. Our constraint is the unit circle, and the simplest way to
    ## implement the constraint function is to use a paraboloid. Thus
    ## the linearization is simply x itself, and the null space is the
    ## hyperplane orthogonal to the current solution estimate.
    #A = r_[[ 2 * x ]]
    A = np.r_[[ gcon ]]

    ## Find Y amd Z, the pseudo-inverse and nullspace of the linearized
    ## constraint coefficients, A.
    Au,As,Avt = np.linalg.svd(A)
    Y = Au * As**-1 * Avt[0]
    Z = Avt[1:]

    ## The value of the constraint function at the point.
    # c_x = dot(x,x) - 1 #quadratic constraint
    c_x = vcon - 1 #conic constraint
    b = -c_x

    ## Find out the Hessian of the Lagrangian
    W = G - lam * Gcon # The constraint hessian (con from "cone")
    ## Calculate the reduced Hessian
    M = dot(dot(Z,W),Z.T)
    ## Calculate the reduced gradient
    m = dot(Z,g + dot(dot(G,Y.T)[:,0],b))
    
    ## Find the restricted step using a trust-region solver.
    y = trust_region_step(M, m, rho, rho * 1e-1)
    if y.shape[0]==1:
        y = np.r_[[y]]
        delta_null = dot(y.T, Z)[0]
    else:
        delta_null = dot(y.T, Z)
    ## # delta = (dot(Y.T,b) + dot(Z.T,y))[:,0]
    # print '-----', y, delta_null
    nd = norm(delta_null)

    ##################################################################
    #### Perform second-order correction (SOC).
    ## Calculate new estimate for c, and consequently b
    c_hat = c_x + dot(gcon,delta_null)\
            + .5 * dot(dot(delta_null, Gcon), delta_null)
    b = -c_hat
    ## Repeat previous steps with new b, and consequently m, to find
    ## the displacement in the cosntraint null space.
    m = dot(Z,g + dot(dot(G, Y.T)[:,0], b))
    y = trust_region_step(M, m, rho, rho * 1e-1)
    if y.shape[0]==1:
        y = r_[[y]]
        delta_null = dot(y.T, Z)[0]
    else:
        delta_null = dot(y.T, Z)
    nd = norm(delta_null)
    delta_q = .5 * dot(dot(y,M),y) + dot(m, y)
    ####
    ##################################################################

    ## Now we calculate the delta in the direction of the "range" of
    ## the local approcximation to the constraint space, and then the
    ## total displacement.
    delta_range = dot(b, Y[0])
    delta = delta_null + delta_range
    ## Calculate the new value of lambda from the result.
    newlam = dot(Y,g)
    newx = x + delta

    filt_pass = True

    ## Use filter to test point, and update filter
    pt = array([func_val(newx,*func_args), np.abs(cons_val(newx)-1)])
    # print pt
    # print filt.values[filt.valid]
    # print filt.dominated(pt, delta_q)
    # print 70*'-'
    if not filt.dominated(pt):
        filt_pass = True
        # rho = min(0.8, rho * 2.0)
        # rho = min(0.8, nd * 2.0)
        filt.add(pt, delta_q, newlam)
    else:
        filt_pass = False
        rho = nd *.25
        newlam = lam
        delta = 0

    return delta, newlam, rho, filt_pass, pt

def filterSQP(x0, lam0, rho0, funcs, args_f, delta_tol=1e-15):
    x = copy(x0)
    delta = zeros(2)
    rho = rho0
    lam = lam0
    val_c, grad_c, hess_c, val_f, grad_f, hess_f = funcs
    Niter = 100
    filt = FletcherFilter()

    Llam = [lam0]
    Lrho = [rho0]
    for it in range(0,Niter):
        delta, lam, rho, filt_pass, pt = SQP_step(x, lam, rho, filt,
                                                 val_c, grad_c, hess_c,
                                                 val_f, grad_f, hess_f, args_f)
        Llam.append(lam)
        Lrho.append(rho)
        if not filt_pass:
            continue
        x = x + delta        

        if np.abs(val_c(x)-1) < 1e-15 and np.abs(delta).sum() < delta_tol:
            break
    return x, val_f(x, *args_f), it, Llam, Lrho
