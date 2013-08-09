#coding: utf-8

from pylab import *
import sys

def get_int_parm(Lpin, Nc):
    '''Get the internal parameters from the bundler output. The first
    half images in one camera (Sony), and the other half is the other
    camera (N8).'''

    f1 = mean(Lpin[:Nc/2,0])
    f2 = mean(Lpin[Nc/2:,0])
    fs1= std(Lpin[:Nc/2,0])
    fs2= std(Lpin[Nc/2:,0])
    k1 = 1e9 * mean(Lpin[:Nc/2,1])/f1**2
    k2 = 1e9 * mean(Lpin[Nc/2:,1])/f2**2
    ks1 = 1e9 * std(Lpin[:Nc/2,1])/f1**2
    ks2 = 1e9 * std(Lpin[Nc/2:,1])/f2**2

    ## Mean and standard deviation from focal length and quadratic distortion coefficient.
    print f1, fs1, k1, ks1
    print f2, fs2, k2, ks2



if __name__ == '__main__':

    with open('bundle.out') as fp:
        aa = fp.readlines()

    Nc = int(aa[1].split(' ')[0])

    Lpin = zeros((Nc, 3))
    for n in range(Nc):
        Lpin[n] = array([float(x) for x in aa[2+n*5].split(' ')])

    ## Get internal parameters from the Bundler output
    get_int_parm(Lpin, Nc)

    ion()
    figure()
    plot(sort(Lpin[:,0]), 'b+')
    grid()
    figure()
    plot(sort(Lpin[:,1]), 'r+')
    grid()
