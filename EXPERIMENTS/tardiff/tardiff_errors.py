#coding: utf-8

import argparse

import json

## Analysis from the Apa dataset --- two sets of 24 images from the same
## buildings, taken with two different cameras and with 90 degree
## rotations at each of the 12 locations.

from pylab import *
import sys

from corisco.quaternion import *
from corisco.quaternion import Quat, matrix_to_quaternion as mtq

def fixit(M):
    M = dot(diag([1,1,1]), M)
    if dot(cross(M[0], M[1]) , M[2]) < 0:
        return dot(M, diag([1,1,-1]))
    else:
        return M

def fix_quaternion(q, r):
    return (q * r.inverse()).canonical() * r

def fix_reference_quaternions_apa(Lr):
    ref = Quat(cos(pi/8),0,sin(pi/8),0)
    for n in range(0,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q
    ref = Quat(cos(pi/8),0,sin(pi/8),0) * Quat(cos(pi/4),0,0,-sin(pi/4))
    for n in range(1,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--transformed_reference', type=str)
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--image_set', type=str)

    args = parser.parse_args()

    image_set = args.image_set

    Lref = loadtxt(args.transformed_reference)

    Nc = len(Lref)

    ## Read the estimated orientations
    data_gen = [json.loads(y) for y in open(args.input_file).readlines()]

    qq = [
        (a,
         mtq(fixit(array(b['matrix'])).T),
         b['time']['edge'] + b['time']['est'])
        for a, b in enumerate(data_gen)
        ]

    print qq

    ## Get orientations ordered by frame number
    Lest = array([Quat(x).q for _, x, _ in sorted(qq)])
    times = array([t for _, _, t in sorted(qq)])

    assert len(Lref) == len(Lest)

    # ## Fix and normalize reference quaternions.
    if image_set == 'apast':
        fix_reference_quaternions_apa(Lest)

    err = sort([(Quat(a) / Quat(b)).canonical().angle() for a,b in zip(Lref, Lest)])

    print (6*'{:>8s}').format('time', 'µ', 'std', 'Q1', 'Q2', 'Q3')
    print (6*'{:>8.3f}').format(mean(times), mean(err), std(err), err[Nc/4], err[Nc/2], err[Nc*3/4])

    ion()
    plot(sort(err), mgrid[0.0:len(err)]/len(err), '-+')

    title('Errors from set {}'.format(args.image_set))
    xlabel('Error [degrees]')
    ylabel('Cumulative frequency')
    grid()

    figure()
    subplot(2,2,1);
    plot(Lref[:,3], 'b-+'); plot(Lest[:,3], 'r-+');
    subplot(2,2,2); plot(Lref[:,1], 'b-+'); plot(Lest[:,1], 'r-+');
    subplot(2,2,3);plot(Lref[:,2], 'b-+');plot(Lest[:,2], 'r-+');
    subplot(2,2,4);plot(Lref[:,3], 'b-+');plot(Lest[:,3], 'r-+');
