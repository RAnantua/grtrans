# Run Grtrans with rrjet model
# The rrjet model is defined in "fluid_model_rrjet.py"
# NOTE -- currently the power law emissivity is very slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python

import numpy as np
import grtrans_batch as gr
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filt

ang=20
name = 'rrjet'+str(ang)
mu = np.cos(ang*np.pi/180.)
size  = 200.
uout = 1./(5*size)
npix = 100
ngeo = 1000

cmperMpc = 3.086e24
MBH = 6.7e9
DTOBH = 16.7*cmperMpc
psize_rg = 2*size/npix
cmperrg = 147708.8 * MBH
psize_cm = psize_rg * cmperrg
psize_rad = psize_cm / DTOBH

RF = 43.e9
cfun = 'jet'
cfun2 =  'seismic'
RERUN = True

def main():
    # run grtrans
    x=gr.grtrans()

    x.write_grtrans_inputs(name+'.in', oname=name+'.out',
                           fname='RRJET',phi0=0.,
                           betaeconst=1.e-4, ximax=10., 
                           nfreq=1,fmin=RF,fmax=RF,
                           gmin=10., p2=2.1, p1=2.1,
                           ename='POLSYNCHPL',
                           nvals=4,
                           spin=0.,standard=1,
                           uout=uout,
                           mbh=MBH,
                           mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,
                           nmu=1,mumin=mu,mumax=mu,
                           gridvals=[-size,size,-size,size],
                           nn=[npix,npix,ngeo],
                           hindf=1,hnt=1,
                           muval=1.)

    if RERUN:
        x.run_grtrans()

    # load image
    x.read_grtrans_output()
    x.convert_to_Jy(DTOBH)

    display_grtrans_image(x)
    #x.disp_grtrans_image()
    #x.disp_pol_map()


def display_grtrans_image(x,nvec=20,veccut=0,blur_kernel=1.1):
    plt.close('all')

    I_im = x.ivals[:,0,0].reshape(npix,npix)
    Q_im = x.ivals[:,1,0].reshape(npix,npix)
    U_im = x.ivals[:,2,0].reshape(npix,npix)
    V_im = x.ivals[:,3,0].reshape(npix,npix)

    I_im = filt.gaussian_filter(I_im, (blur_kernel, blur_kernel))
    Q_im = filt.gaussian_filter(Q_im, (blur_kernel, blur_kernel))
    U_im = filt.gaussian_filter(U_im, (blur_kernel, blur_kernel))
    V_im = filt.gaussian_filter(V_im, (blur_kernel, blur_kernel))

    # convert to Tb
    factor = 3.254e13/(RF**2 * psize_rad**2)
    I_im *= factor
    Q_im *= factor
    U_im *= factor
    V_im *= factor
    
    # Polarization Vectors
    P_im = np.abs(Q_im + 1j*U_im)
    m_im = P_im/I_im


    thin = npix//nvec
    mask = I_im > veccut * np.max(I_im)
    mask2 = mask[::thin, ::thin]

    m = m_im[::thin, ::thin][mask2]
    x = (np.array([[i for i in range(npix)] for j in range(npix)])[::thin, ::thin])
    x = x[mask2]
    y = (np.array([[j for i in range(npix)] for j in range(npix)])[::thin, ::thin])
    y = y[mask2]
    a = (-np.sin(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    a = a[mask2]
    #a = m*a
    b = ( np.cos(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    b = b[mask2]
    #b = m*b

    P_im[np.logical_not(mask)]=0.
    m_im[np.logical_not(mask)]=0.

    # ticks
    xticks = ticks(npix, 2*size/npix)
    yticks = ticks(npix, 2*size/npix)

    # display Stokes I 
    plt.figure(0)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes Q 
    plt.figure(1)
    im = plt.imshow(Q_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes Q, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes U
    plt.figure(2)
    im = plt.imshow(U_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes U, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes V 
    plt.figure(3)
    im = plt.imshow(V_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes V, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display P
    plt.figure(4)
    im = plt.imshow(P_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("P, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display m
    plt.figure(5)
    im = plt.imshow(m_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("P/I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display I with pol ticks
    plt.figure(6)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.01*npix, units='x', pivot='mid', color='k', angles='uv', 
           scale=1.0/thin)
    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.005*npix, units='x', pivot='mid', color='w', angles='uv', 
           scale=1.1/thin)

    plt.show()

def ticks(axisdim, psize, nticks=8):
    """Return a list of ticklocs and ticklabels
       psize should be in desired units
    """

    axisdim = int(axisdim)
    nticks = int(nticks)
    if not axisdim % 2: axisdim += 1
    if nticks % 2: nticks -= 1
    tickspacing = float((axisdim-1))/nticks
    ticklocs = np.arange(0, axisdim+1, tickspacing) - 0.5
    ticklabels= np.around(psize * np.arange((axisdim-1)/2.0, -(axisdim)/2.0, -tickspacing), decimals=1)

    return (ticklocs, ticklabels)

if __name__=='__main__':
    main()


