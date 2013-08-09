#coding: utf-8


## Analysis from the Apa dataset --- two sets of 24 images from the same
## buildings, taken with two different cameras and with 90 degree
## rotations at each of the 12 locations.

from pylab import *
import sys

from quaternion import *

from quaternion import matrix_to_quaternion as mtq




def fix_quaternion(q,r):
    return (q * r.inverse()).canonical() * r

def solve_transformation(Lq, Lr):
    Nc = Lq.shape[0]
    M = zeros((4*Nc,4))
    b = zeros(4*Nc)

    for nn in range(Nc):
        n = nn
        # M[nn*4+0] = array([ Lq[n,0],-Lq[n,1],-Lq[n,2],-Lq[n,3]])
        # M[nn*4+1] = array([ Lq[n,1], Lq[n,0], Lq[n,3],-Lq[n,2]])
        # M[nn*4+2] = array([ Lq[n,2],-Lq[n,3], Lq[n,0], Lq[n,1]])
        # M[nn*4+3] = array([ Lq[n,3], Lq[n,2],-Lq[n,1], Lq[n,0]])
        M[nn*4+0] = array([ Lq[n,0],-Lq[n,1],-Lq[n,2],-Lq[n,3]])
        M[nn*4+1] = array([ Lq[n,1], Lq[n,0],-Lq[n,3], Lq[n,2]])
        M[nn*4+2] = array([ Lq[n,2], Lq[n,3], Lq[n,0],-Lq[n,1]])
        M[nn*4+3] = array([ Lq[n,3],-Lq[n,2], Lq[n,1], Lq[n,0]])
        b[nn*4:nn*4+4] = Lr[n]

    #print 'Parâmetros para o mapeamento:'
    lres = lstsq(M,b)
    #print Quat(lres[0]).normalize()
    tt = Quat(lres[0]).normalize()
    return tt


def get_transformed_reference(C, D):
    N = dot(dot(inv(dot(C.T,C)), C.T), Lr)
    return dot(C,N), N
    




def get_int_parm(Lpin, Nc):
    ion()
    plot(sort(Lpin[:,0]), '+')
    grid()

    f1 = mean(Lpin[:Nc/2,0])
    f2 = mean(Lpin[Nc/2:,0])
    fs1= std(Lpin[:Nc/2,0])
    fs2= std(Lpin[Nc/2:,0])
    k1 = 1e9 * mean(Lpin[:Nc/2,1])/f1**2
    k2 = 1e9 * mean(Lpin[Nc/2:,1])/f2**2
    ks1 = 1e9 * std(Lpin[:Nc/2,1])/f1**2
    ks2 = 1e9 * std(Lpin[Nc/2:,1])/f2**2

    print f1,fs1, k1,ks1
    print f2,fs2, k2, ks2



if __name__ == '__main__':

    aa = open(sys.argv[1]).readlines()


    it = int(sys.argv[2])
    gs = int(sys.argv[3])


    
    Nc = int(aa[1].split(' ')[0])

    Lpin = zeros((Nc, 3))
    for n in range(Nc):
        Lpin[n] = array([float(x) for x in aa[2+n*5].split(' ')])

    Lori = zeros((Nc, 3,3))
    for n in range(Nc):
        for j in range(3):
            Lori[n,j] = array([float(x) for x in aa[3+n*5+j].split(' ')])

    ## Get internal parameters from the Bundler output
    get_int_parm(Lpin, Nc)
    # myxtplk()


    #trans = Quat(1.0,0,0,0) # Quat( 0.9377, -0.0749,  0.3391, -0.0091)
    #trans = Quat( 0.9389, -0.0302,  0.3427,  0.0070)
    trans = Quat( 1.0,0,0,0)
    # trans = Quat(cos(pi/8),0,sin(pi/8),0)
    #trans = Quat( 0.9389, -0.0345,  0.3425,  0.0056)  # 10000 001
    #trans = Quat( 0.9414, -0.0524,  0.3314, -0.0347)
    #trans = Quat( 0.9345, -0.0153,  0.3552, -0.0179)
    
    Lq = zeros((Nc,4))
    for n in range(Nc):
        ori = trans * Quat(mtq(Lori[n]))
        Lq[n] = ori.normalize().q * array([1,-1,1,1])





    ## Read the Corisco output
    b1 = loadtxt('exps-apa09/it%d_gspc%03d.json.out'%(it,gs))
    b2 = loadtxt('exps-apa08/it%d_gspc%03d.json.out'%(it,gs))
    b1 = b1[argsort(b1[:,0])]
    b2 = b2[argsort(b2[:,0])]


    Lr = zeros((Nc,4))
    Lr[:Nc/2] = b1[:,5:9]
    Lr[Nc/2:] = b2[:,5:9]


    ref = Quat(cos(pi/8),0,sin(pi/8),0)
    for n in range(0,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q

    ref = Quat(cos(pi/8),0,sin(pi/8),0) * Quat(cos(pi/4),0,0,-sin(pi/4))
    for n in range(1,Nc,2):
        ori = fix_quaternion(Quat(Lr[n]), ref)
        Lr[n] = ori.normalize().q






    Le = zeros((Nc,4))
    La = zeros(Nc)

    ## Find transformation to take reference to Corisco reference frame.
    # Lr_hat, N = get_transformed_reference(Lq, Lr)
    # for n in range(Nc):
    #     Lr_hat[n] = Quat(Lr_hat[n]).normalize().q

    tt = solve_transformation(Lq, Lr)
    Lr_hat = zeros((Nc,4))
    for n in range(Nc):
        Lr_hat[n] = (tt * Quat(Lq[n])).q


    ## Calculate errors between each estimated orientation and
    ## transformed reference.
    for n in range(Nc):
        Le[n] = ( Quat(Lr_hat[n]) / Quat(Lr[n]) ).canonical().q
        # Le[n] = ( Quat(Lr[n]) / Quat(Lr_hat[n]) ).canonical().q
        #Le[n] = (Quat(Lq[n]) * Quat(Lr[n]).inverse()).canonical().q
        La[n] = Quat(Le[n]).angle()


    ion()


    ## Plot reference orientations
    figure(1, figsize=(10,6))
    suptitle(u'Referência transformada e orientações estimadas')

    subplot(2,2,1)
    title(u'Conjunto paisagem (Bundler)')
    plot(Lr_hat[::2], '-+')
    # plot(Lq[::2], '-+')
    # plot(Lq[:Nc/2:2], '-+')
    # plot(Lq[Nc/2::2], '-x')
    grid()
    ylim(-1,1)
    # xlabel(u'Índice da imagem')
    ylabel(u'Componentes do quaterno')

    legend(['a','b','c','d'], ncol=2, loc='lower center')


    subplot(2,2,2)
    title(u'Conjunto retrato (Bundler)')
    plot(Lr_hat[1::2],'-+')
    # plot(Lq[1::2],'-+')
    # plot(Lq[1:Nc/2:2], '-+')
    # plot(Lq[1+Nc/2::2], '-x')
    grid()
    ylim(-1,1)
    # xlabel(u'Índice da imagem')

    # savefig('apa_ref.pdf')


    ## Plot estimated orientations
    # figure(2, figsize=(10,3))
    # suptitle(u'Estimativas (Corisco)')
    subplot(2,2,3)
    title(u'Conjunto paisagem (Corisco)')
    plot(Lr[::2], '-+')
    # plot(Lr[:Nc/2:2], '-+')
    # plot(Lr[Nc/2::2], '-x')
    grid()
    ylim(-1,1)
    xlabel(u'Índice da imagem')
    ylabel(u'Componentes do quaterno')

    subplot(2,2,4)
    title(u'Conjunto retrato (Corisco)')
    plot(Lr[1::2], '-+')
    # plot(Lr[1:Nc/2:2], '-+')
    # plot(Lr[1+Nc/2::2], '-x')
    grid()
    ylim(-1,1)
    xlabel(u'Índice da imagem')

    # savefig('apa_est.pdf')
    savefig('apa_params.pdf')





    # myzzzplkz()
    ## Plot residual orientations
    figure(3, figsize=(10,6))
    subplot(1,2,1)
    plot(Le[:Nc/2:2], '-+')
    plot(Le[Nc/2::2], '-x')
    grid()
    ylim(-1,1)

    subplot(1,2,2)
    plot(Le[1:Nc/2:2], '-+')
    plot(Le[1+Nc/2::2], '-x')
    grid()
    ylim(-1,1)


    figure(4,figsize=(8*.8,6*.8))
    title(u'Distribuição acumulada de erros')
    xlabel('Erro [graus]')
    ylabel(u'Frequência acumulada')
    #plot(sort(La), mgrid[0.0:Nc]/Nc, '-+')
    semilogx(sort(La), mgrid[0.0:Nc]/Nc, '-+')
    grid()
    savefig('apa_cdf.pdf')





