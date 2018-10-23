import grtrans_batch as gr
import matplotlib.pyplot as plt
import numpy as np

x=gr.grtrans()
#x.compile_grtrans()

mu=np.cos(30*np.pi/180.)

############################################################################################
#Image

x=gr.grtrans()
x.write_grtrans_inputs('koralm87.in', oname='koralm87.out',
                       sigcut=25,
                       fname='KORAL',
                       nfreq=1,fmin=230.e9,fmax=230.e9,
                       ename='SYNCHTHAV',
                       nvals=1,
                       spin=0.9375,standard=1,
                       uout=0.01,
                       mbh=4.e6, 
                       mdotmin=1.57e15, mdotmax=1.57e15, nmdot=1,
                       nmu=1,mumin=mu,mumax=mu,
                       gridvals=[-30,30,-30,30],
                       nn=[256,256,750],
                       hhfile='koralm87_001',hdfile='koralm87_',
                       hindf=1,hnt=1,
                       muval=1.)       
x.run_grtrans()
x.read_grtrans_output()
x.convert_to_Jy(8*gr.cmperkpc)
#x.ivals /= np.sum(x.ivals)
#x.ivals = np.sqrt(x.ivals**2)**.33
x.disp_grtrans_image(0)
#plt.clim(vmin=0,vmax=.1)
#plt.clim(vmin=0,vmax=.0009)
plt.clim(vmin=0,vmax=0.03)
plt.set_cmap('hot')
plt.savefig('koralm87_2dim.png')

############################################################################################
#Spectrum
x=gr.grtrans()
x.write_grtrans_inputs('koralm87_spec.in', oname='koralm87_spec.out',
                       sigcut=25,
                       fname='KORAL',
                       nfreq=25,fmin=1.e8,fmax=1.e24,
                       ename='SYNCHTHBREMS',
                       nvals=1,
                       spin=0.9375,standard=1,
                       uout=0.01,
                       mbh=4.e6, 
                       mdotmin=1.57e15, mdotmax=1.57e15, nmdot=1,
                       nmu=1,mumin=mu,mumax=mu,
                       gridvals=[-30,30,-30,30],
                       nn=[50,50,500],
                       hhfile='koralm87_001',hdfile='koralm87_',
                       hindf=1,hnt=1,
                       fscalefac=1.,
                       muval=1.)       
x.run_grtrans()
x.read_grtrans_output()
x.convert_to_lum()
spec  = x.spec[0]*x.freqs
plt.plot(x.freqs,spec,'r-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\\nu$ (Hz)', size=20)
plt.ylabel('$\\nu L_{\\nu}$ (erg s$^{-1}$)', size=20)
plt.xlim([1.e8,1.e24])
plt.ylim([1.e28,1.e37])
plt.savefig('koralm87_2dspec.png')

