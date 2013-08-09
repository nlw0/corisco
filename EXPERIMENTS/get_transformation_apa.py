#coding: utf-8



import json

## Analysis from the Apa dataset --- two sets of 24 images from the same
## buildings, taken with two different cameras and with 90 degree
## rotations at each of the 12 locations.

from pylab import *
import sys

from corisco.quaternion import *
from corisco.quaternion import Quat, matrix_to_quaternion as mtq

def fix_quaternion(q,r):
    return (q * r.inverse()).canonical() * r

def solve_transformation(Lq, Lr):
    Nc = Lq.shape[0]
    M = zeros((4*Nc,4))
    b = zeros(4*Nc)

    for nn in range(Nc):
        n = nn
        M[nn*4+0] = array([ Lq[n,0],-Lq[n,1],-Lq[n,2],-Lq[n,3]])
        M[nn*4+1] = array([ Lq[n,1], Lq[n,0],-Lq[n,3], Lq[n,2]])
        M[nn*4+2] = array([ Lq[n,2], Lq[n,3], Lq[n,0],-Lq[n,1]])
        M[nn*4+3] = array([ Lq[n,3],-Lq[n,2], Lq[n,1], Lq[n,0]])
        b[nn*4:nn*4+4] = Lr[n]

    lres = lstsq(M,b)
    return Quat(lres[0]).normalize()


def load_json_dump(filename, imgset, it, gs):
    with open(filename) as fp:
        data = [json.loads(x) for x in fp.readlines()]
    
    return [[x['input']['frame']] + x['ori_est']
            for x in data
            if x['input']['set'] == imgset
            and x['input']['ransac_itr'] == it
            and x['input']['grid']['size'] == gs]

def fix_reference_quaternions(Lr):
    ref = Quat(cos(pi/8),0,sin(pi/8),0)
    for n in range(0,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q
    ref = Quat(cos(pi/8),0,sin(pi/8),0) * Quat(cos(pi/4),0,0,-sin(pi/4))
    for n in range(1,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q


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
    it = 10000
    gs = 1

    # b1 = array(load_json_dump('solutions_apa.dat', 'set-apa09.json', it, gs))
    # b2 = array(load_json_dump('solutions_apa.dat', 'set-apa08.json', it, gs))
    b1 = array(load_json_dump('solutions_apa.dat', 'set-apa09.json', it, gs))
    b2 = array(load_json_dump('solutions_apa.dat', 'set-apa08.json', it, gs))

    ## Get orientations ordered by frame number
    Lr = r_[b1[argsort(b1[:,0]), 1:], b2[argsort(b2[:,0]), 1:]]

    ## fix and normalize reference quaternions.
    fix_reference_quaternions(Lr)

    Le = zeros((Nc,4))
    La = zeros(Nc)

    ## Find transformation to take alternative orientations to Corisco
    ## reference frame.
    tt = solve_transformation(Lq, Lr)
    Lr_hat = zeros((Nc,4))
    for n in range(Nc):
        Lr_hat[n] = (tt * Quat(Lq[n])).q

    dd = [{"set": s, "frame": fr}
          for s in ["set-apa09.json", "set-apa08.json"]
          for fr in range(24)]

    #print zip(dd, Lr_hat)
    tref = [dict(meta.items() + [("ori_ref", ori.tolist())])
           for meta, ori in zip(dd, Lr_hat)]

    with open('transformed_reference_apa.dat', 'w') as fp:
        json.dump(tref, fp, indent=4)

    ## Calculate errors between each estimated orientation and
    ## transformed reference.
    for n in range(Nc):
        Le[n] = ( Quat(Lr_hat[n]) / Quat(Lr[n]) ).canonical().q
        La[n] = Quat(Le[n]).angle()

    ion()

    figure(4,figsize=(8*.8,6*.8))
    title(u'Distribuição acumulada de erros')
    xlabel('Erro [graus]')
    ylabel(u'Frequência acumulada')
    plot(sort(La), mgrid[0.0:Nc]/Nc, '-+')
    # semilogx(sort(La), mgrid[0.0:Nc]/Nc, '-+')
    grid()
    # savefig('apa_cdf.pdf')


    ## Plot reference orientations
    figure(5, figsize=(10,6))
    suptitle(u'Referência transformada e orientações estimadas')

    subplot(1,2,1)
    title(u'Conjunto paisagem')
    plot(Lr_hat[::2], 'b-+')
    plot(Lr[::2], 'r-+')
    ylim(-1,1)
    xlabel(u'Índice da imagem')
    ylabel(u'Componentes do quaterno')
    grid()

    subplot(1,2,2)
    title(u'Conjunto retrato')
    plot(Lr_hat[1::2],'b-+')
    plot(Lr[1::2], 'r-+')
    ylim(-1,1)
    xlabel(u'Índice da imagem')
    grid()
