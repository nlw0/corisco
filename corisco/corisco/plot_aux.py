from pylab import *

def pplot(x,*a,**kw):
    if len(x.shape)>1:
        plot(x[:,0], x[:,1], *a, **kw)
    else:
        plot(x[0], x[1], *a, **kw)

## Colors used for plotting different kels (directions)
dir_colors=['#ea6949', '#51c373', '#a370ff', '#444444']

def plot_solution(Lines,edgels,normals,Q,labels=[]):
    figure(1)
    mul = 0.02
    for k in [3,0,1,2]:
        for p in Lines[k]:
            px,py,ux,uy = p[:,:4].T
            plot(px,py, '-+', color=dir_colors[k])
            plot(px+mul*ux,py+mul*uy, '+', color=dir_colors[k])

    axis('equal')
    xlim(-1.1, 1.1)
    ylim( 1.1,-1.1)
    grid()

    t = mgrid[-pi:pi+0.01:0.01]

    figure(2)
    mul = 0.02
    for px,py,ux,uy in edgels:
            plot([px-mul*uy,px+mul*uy],[py+mul*ux,py-mul*ux], 'r-' )
    axis('equal')

    xlim(-1.1, 1.1)
    ylim( 1.1,-1.1)
    grid()


    r_est = Q.rot().T

    figure(3)
    plot_spheres(r_est, normals,labels)

def plot_spheres(q, normals, labels=[]):
    r = q.rot().T
    t = mgrid[-pi:pi+0.01:0.01]
    circ = c_[sin(t), cos(t)]

    cx = dot(c_[zeros(t.shape[0]), sin(t), cos(t)], r)
    cy = dot(c_[sin(t), zeros(t.shape[0]), cos(t)], r)
    cz = dot(c_[sin(t), cos(t), zeros(t.shape[0])], r)

    suptitle('Observed normals and fitted vanishing points on the Gaussian sphere')

    subplot(2,2,1)
    if labels==[]:
        pplot(normals[:,[0,1]], 'x', alpha=0.5)
    else:
        for lab in [3,0,1,2]:
            pplot(normals[labels==lab][:,[0,1]], 'x', alpha=0.1, color=dir_colors[lab])
    pplot(circ, 'k-');    pplot(cx[:,[0,1]], 'k--',color=dir_colors[0]);
    pplot(cy[:,[0,1]], 'k--',color=dir_colors[1]);    pplot(cz[:,[0,1]], 'k--',color=dir_colors[2])
    plot([0,0],[-1,1],'k-');    plot([-1,1],[0,0],'k-')
    for k in range(3):
        plot([0,r[k,0]],[0,r[k,1]], '-', color=dir_colors[k])
    axis('equal')
    axis([-1.2,1.2,1.2,-1.2])
    grid()

    subplot(2,2,2)
    if labels==[]:
        pplot(normals[:,[2,1]], 'x', alpha=0.5)
    else:
        for lab in [3,0,1,2]:
            pplot(normals[labels==lab][:,[2,1]], 'x', alpha=0.1, color=dir_colors[lab])
    pplot(circ, 'k-');    pplot(cx[:,[2,1]], 'k--',color=dir_colors[0]);
    pplot(cy[:,[2,1]], 'k--',color=dir_colors[1]);    pplot(cz[:,[2,1]], 'k--',color=dir_colors[2])
    plot([0,0],[-1,1],'k-');    plot([-1,1],[0,0],'k-')
    for k in range(3):
        plot([0,r[k,2]],[0,r[k,1]], '-', color=dir_colors[k])
    axis('equal')
    axis([1.2,-1.2,1.2,-1.2])
    grid()
    
    subplot(2,2,3)
    if labels==[]:
        pplot(normals[:,[0,2]], 'x', alpha=0.5)
    else:
        for lab in [3,0,1,2]:
            pplot(normals[labels==lab][:,[0,2]], 'x', alpha=0.1, color=dir_colors[lab])
    pplot(circ, 'k-');    pplot(cx[:,[0,2]], 'k--',color=dir_colors[0]);
    pplot(cy[:,[0,2]], 'k--',color=dir_colors[1]);    pplot(cz[:,[0,2]], 'k--',color=dir_colors[2])
    plot([0,0],[-1,1],'k-');    plot([-1,1],[0,0],'k-')
    for k in range(3):
        plot([0,r[k,0]],[0,r[k,2]], '-', color=dir_colors[k])
    axis('equal')
    axis([-1.2,1.2,-1.2,1.2])
    grid()



def plot_edgels(ax,pic):
    ax.imshow(pic.frame/260., interpolation='nearest')
    pic.plot_edgels(ax,2.5)
    ax.axis([0,pic.Iwidth,pic.Iheight,0])

