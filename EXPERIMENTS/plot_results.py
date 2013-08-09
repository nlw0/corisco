#coding: utf-8



import json

## Analysis from the Apa dataset --- two sets of 24 images from the same
## buildings, taken with two different cameras and with 90 degree
## rotations at each of the 12 locations.

from pylab import *
import sys

from corisco.quaternion import *
from corisco.quaternion import Quat, matrix_to_quaternion as mtq


def error_from_lists(q1, q2):    
    return (Quat(*q1) / Quat(*q2)).canonical().angle()

def get_errors(exps, refs, imgset, it, gs, smooth=None):

    rr = {x['frame']: x['ori_ref']
          for x in refs
          if x['set'] == imgset}
    
    # return [ (x['ori_est'], rr[x['input']['frame']])
    result =  [ error_from_lists(x['ori_est'], rr[x['input']['frame']])
                for x in exps
                if x['input']['set'] == imgset
                and x['input']['ransac_itr'] == it
                and x['input']['grid']['size'] == gs
                and (smooth is None or x['input']['smooth'] == smooth)
                ]
    if not result:
        raise Exception
    else:
        return result


def get_errorbar(qq):
    N = len(qq)
    sqq = sorted(qq)
    return sqq[0], sqq[N/4], sqq[N/2], sqq[N*3/4], sqq[-3]


if __name__ == '__main__':

    with open('bundle.out') as fp:
        aa = fp.readlines()

    Nc = int(aa[1].split(' ')[0])

    Lori = zeros((Nc, 3,3))
    for n in range(Nc):
        for j in range(3):
            Lori[n,j] = array([float(x) for x in aa[3+n*5+j].split(' ')])
    
    Lq = zeros((Nc,4))
    for n in range(Nc):
        ori = Quat(mtq(Lori[n]))
        Lq[n] = ori.normalize().q * array([1,-1,1,1])

    ## Read the Corisco output
    it = 200

    with open('transformed_reference_apa.dat') as fp:
        refs = json.load(fp)

    with open('solutions_apa.dat') as fp:
        exps = [json.loads(x) for x in fp.readlines()]

    qq = []
    for gs in 2** mgrid[0:8]:
        all_errors = (get_errors(exps, refs, 'set-apa09.json', it, gs, smooth=1.0) +
                      get_errors(exps, refs, 'set-apa08.json', it, gs))

        print all_errors
        qq.append(get_errorbar(all_errors))

        

    qq = array(qq)



    ion()
    plot(qq, 'r-o')
    grid()

