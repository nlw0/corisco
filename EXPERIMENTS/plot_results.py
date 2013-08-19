#coding: utf-8



import json

## Analysis from the Apa dataset --- two sets of 24 images from the same
## buildings, taken with two different cameras and with 90 degree
## rotations at each of the 12 locations.

from pylab import *
import sys

from corisco.quaternion import *
from corisco.quaternion import Quat, matrix_to_quaternion as mtq

def draw_errbox(ax, x, y_input, width=1.0, outliers=0, draw_mean=False, lw=3):
    y = sorted(y_input)

    w = width * 0.5
    wm = w*0.8

    Ny = len(y)
    vmean = mean(y)
    vmin = y[0] ## Minimum value
    vb = y[Ny / 4] ## First quartile
    vm = y[Ny / 2] ## Median
    vt = y[Ny * 3 / 4] ## Third quartile
    vmax = y[-1-outliers] ## Maximum, ignoring outliers


    LL=[]

    LL.append(Line2D([vb, vt, vt, vb, vb],
                     [x-w, x-w, x+w, x+w, x-w], lw=lw, color='#0000ff') )
    LL.append(Line2D([vmax, vmax],[x-wm, x+wm], lw=lw, color='#0000ff') )
    LL.append(Line2D([vmin, vmin],[x-wm, x+wm], lw=lw, color='#0000ff') )
    LL.append(Line2D([vt, vmax],[x, x], lw=lw, ls='--', color='#0000ff') )
    LL.append(Line2D([vb, vmin],[x, x], lw=lw, ls='--', color='#0000ff') )

    LL.append(Line2D([vm, vm],[x-w, x+w], lw=lw, color='#ff0000') )
    LL.append(Line2D([vm],[x],color='#ff0000', mew=2, ms=7, marker='+') )
    if draw_mean:
        LL.append(Line2D([vmean],[x],color='#ff0000', mew=2, ms=7, marker='x') )

    if outliers > 0:
        LL.append(Line2D(y[-outliers:], x * ones(outliers),
                         color='#000000', ls='', mew=1.5, ms=4, marker='x') )

    ## Actually plot all the desired lines.
    for ll in LL:
        ax.add_line(ll)

def get_errors(exps, rr, imgset, it, gs, smooth=None):
    rr = refs[imgset]

    ## Create dictionary with all estiamted references, selecting just
    ## the specified images.
    qq = {x['input']['frame']: Quat(*x['ori_est'])
          for x in exps
          if x['input']['set'] == imgset
          and x['input']['ransac_itr'] == it
          and x['input']['grid']['size'] == gs
          and (smooth is None or x['input']['smooth'] == smooth)}
    
    assert rr.keys() == qq.keys()

    return [(qq[f] / rr[f]).canonical().angle()
            for f in rr.keys()]

def get_times(exps, imgset, it, gs, smooth=None):
    ## Create dictionary with all estiamted references, selecting just
    ## the specified images.
    times = [
            (x['time']['total'], x['time']['ransac'])
            for x in exps
            if x['input']['set'] == imgset
            and x['input']['ransac_itr'] == it
            and x['input']['grid']['size'] == gs
            and (smooth is None or x['input']['smooth'] == smooth)
            ]
    return times

def load_reference_orientations(filename):
    with open(filename) as fp:
        ref_db = json.load(fp)

    image_sets = {x['set'] for x in ref_db}

    ## Create a dictionary with reference orientation for each frame.
    return {imgs: {x['frame']: Quat(*x['ori_ref'])
                   for x in ref_db
                   if x['set'] == imgs}
            for imgs in image_sets}

if __name__ == '__main__':

    fset = 'ApaSt'
    it = int(sys.argv[1])
            
    font = {'family' : 'bitstream vera sans', # 'droid sans',
            'weight' : 'normal',
            'size'   : 14}

    matplotlib.rc('font', **font)

    ref_db_filename = 'transformed_reference_apa.dat'
    refs = load_reference_orientations(ref_db_filename)

    with open('solutions_apa.dat') as fp:
        exps = [json.loads(x) for x in fp.readlines()]
    
    result_errors = {
        gs: get_errors(exps, refs, it, gs, smooth=1.0) +
             get_errors(exps, refs, 'set-apa08.json', it, gs))
        for gs in 2 ** mgrid[0:8]
    }
    
    result_times = {
        gs: (get_times(exps, 'set-apa09.json', it, gs, smooth=1.0) +
             get_times(exps, 'set-apa08.json', it, gs))
        for gs in 2 ** mgrid[0:8]
    }

    tt = array([(array(result_times[gs])[:,1].mean(),
                 array(result_times[gs])[:,0].mean(),
                 array(result_times[gs])[:,0].mean() - array(result_times[gs])[:,0].min(),
                 array(result_times[gs])[:,0].max() - array(result_times[gs])[:,0].mean())
                for gs in 2 ** mgrid[0:8]
                ])

    ion()
    figure(2, figsize=(10,5))

    suptitle(u'Corisco performance on {}, $C_r$={} RANSAC iterations'.format(fset, it))

    ## Angular errors graphic
    ax = axes([0.07,0.1,0.6,0.79])
    title(u'Error distribution')
    ylabel(u'Grid spacing $C_g$ [pixels]')
    xlabel(u'Error [degrees]')

    for gs in mgrid[0:8]:
        draw_errbox(ax, gs, result_errors[2**gs], width=0.4, outliers=4)

    yticks(mgrid[:8], 2**mgrid[:8])
    grid()

    xlim(0,20)

    ## Time graphic
    bx=axes([0.69,0.1,0.3,0.79], sharey=ax)
    bx.set_xscale('log', basex=10)
    title(u'Process duration')
    xlabel('Time [seconds]')

    l2=plot(tt[:,0], mgrid[:8], 'bd--')[0]
    l1=plot(tt[:,1], mgrid[:8], 'ro-')[0]
    errorbar(tt[:,1], mgrid[:8], xerr=tt[:,2:4].T, lw=2, color='r')

    legend([l1,l2], ['Total','RANSAC'], loc='lower right', numpoints=1)

    yticks(mgrid[0:8], [])
    xticks([0.1,1,10,100,1000],[0.1,1,10,100,1000])
    grid()

    if it == 10000:
        xlim(1.8, 500)
    elif it == 1000:
        xlim(0.18, 50.0)
    elif it == 200:
        xlim(0.18, 50.0)


    ylim(7.3, -0.3)

    savefig('yahyah.pdf', dpi=100)


    for gs in 2**array([0,2,5]):
        err = array(result_errors[gs])
        print gs
        print 'mu: {} std{}'.format(err.mean(), err.std())

