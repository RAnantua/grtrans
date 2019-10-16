# Run Grtrans with rrjet model
# The rrjet model is defined in "fluid_model_rrjet.py"
# NOTE -- currently the power law emissivity is very slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python

import numpy as np
import grtrans_batch as gr

ang=20
name = 'rrjet'+str(ang)
mu = np.cos(ang*np.pi/180.)
size  = 200.
uout = 1./(5*size)

x=gr.grtrans()

x.write_grtrans_inputs(name+'.in', oname=name+'.out',
                       fname='RRJET',phi0=0.,
                       betaeconst=1.e-4, ximax=10., 
                       nfreq=1,fmin=43.e9,fmax=43.e9,
                       gmin=10., p2=2.1, p1=2.1,
                       ename='POLSYNCHPL',
                       nvals=4,
                       spin=0.,standard=1,
                       uout=uout,
                       mbh=6.7e9, 
                       mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,
                       nmu=1,mumin=mu,mumax=mu,
                       gridvals=[-size,size,-size,size],
                       nn=[100,100,1000],
                       hindf=1,hnt=1,
                       muval=1.)

x.run_grtrans()
x.read_grtrans_output()
x.disp_grtrans_image()
#x.disp_pol_map()
