#!/usr/bin/python2.7
#coding: utf-8

## Old silly test script, should turn into something more interesting
## some day.

import argparse
import code
import pylab as plt

from filtersqp.fletcherfilter import FletcherFilter

def main():
    plt.ion()

    fil = FletcherFilter()
    Niter = 12
    logp = plt.zeros((Niter,2))
    for k in range(Niter):
        while True:
            #print k
            p = plt.rand(2)
            if not fil.dominated(p):
                break
        logp[k] = p
        fil.add(p, 0.0, 0.0)
        ff = fil.values[fil.valid]
        ff = plt.r_[[[1e-6,1]], ff[plt.argsort(ff[:,0])], [[1,1e-6]]]
        ww = plt.zeros((ff.shape[0] * 2 - 1, 2))
        ww[::2] = ff
        ww[1::2,0] = ff[1:,0]
        ww[1::2,1] = ff[:-1,1]
        plt.loglog(ww[:,0], ww[:,1], '-')
    plt.loglog(logp[:,0], logp[:,1], 'ys-', lw=2)
    plt.axis([0,1,0,1])
    plt.axis('equal')
    plt.grid()
        
    code.interact()
