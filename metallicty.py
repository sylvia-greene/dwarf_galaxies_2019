import numpy as np
import astropy.io.fits as pyfits
import pylab
pylab.ion()
import pandas as pd
import os

new_dataframe = pd.DataFrame(
    { }
    )

pylab.ion()


infile = '/afs/hep.wisc.edu/home/kbechtol/public/des/external_data/s5/dr1/s5cat_1700d_dr1_4.fits'



reader = pyfits.open(infile)
data = reader[1].data
reader.close()

# Select good quality stellar measurements
# Based on Section 5.3 of S5 overview paper
cut = (data['good_star'] == 1) \
      & (data['sn_1700d'] > 5.) \
      & (data['feh_std'] < 0.5)

# Dereddened photometry following DES DR1 prescription
mag = {}
mag['g'] = data['decam_g'] - 3.186 * data['ebv']
mag['r'] = data['decam_r'] - 2.140 * data['ebv']
mag['i'] = data['decam_i'] - 1.569 * data['ebv']
mag['z'] = data['decam_z'] - 1.196 * data['ebv']

pylab.figure()
pylab.scatter(mag['g'][cut] - mag['r'][cut],
              mag['r'][cut] - mag['z'][cut],
              c=data['feh50'][cut],
              s=2, vmin=-3.5, vmax=0.5)
pylab.colorbar(label='[Fe/H]')
pylab.xlim(0.2, 1.1)
pylab.ylim(0.0, 0.6)
pylab.xlabel('g - r (mag)')
pylab.ylabel('r - z (mag)')

d= (1/data.parallax)
pylab.figure()
pylab.scatter(mag['g'][cut] - mag['r'][cut],
              mag['r'][cut] - mag['z'][cut], 
              c=d[cut], s=2,
              vmin=0.,vmax=15.)
pylab.colorbar(label='distance')
pylab.xlim(0.2, 1.1)
pylab.ylim(0.0, 0.6)
pylab.xlabel('g - r (mag)')
pylab.ylabel('r - z (mag)')



pylab.figure()
pylab.scatter(mag['g'][cut] - mag['r'][cut],
              mag['r'][cut] - mag['z'][cut], 
              c=data.pmra[cut], s=2,
              vmin=-13.,vmax=15.)
pylab.colorbar(label='pmra')
pylab.xlim(0.2, 1.1)
pylab.ylim(0.0, 0.6)
pylab.xlabel('g - r (mag)')
pylab.ylabel('r - z (mag)')



pylab.figure()
band=( mag['g'] - mag['r'])
pylab.scatter(data.pmra[cut],data.pmdec[cut], c=band[cut], s=2., vmin=0., vmax=.9)
pylab.colorbar(label='g-r(mag)')
pylab.xlabel('pmra')
pylab.ylabel('pmdec')
pylab.xlim(-6.,10.)
pylab.ylim(-8.,4.5)




pylab.figure()
pylab.scatter(mag['g'],data.pmra_error, s=2.)
pylab.xlabel('g(mag)')
pylab.ylabel('pmra error')
pylab.scatter(mag['g'][cut], data.pmra_error[cut], c='red', s=2.)



pylab.figure()
pylab.scatter(data.pmra,data.pmdec, c=data.pmra_error, vmin=0., vmax=1., s=2)
pylab.colorbar(label='pmra error')
pylab.xlabel('pmra')
pylab.ylabel('pmdec')
pylab.xlim(-7.5,10.5)
pylab.ylim(-8.5,4.)

pylab.figure()
pylab.scatter(data.pmra[cut],data.pmdec[cut], c=data.pmra_error[cut], vmin=0., vmax=1., s=2)
pylab.colorbar(label='pmra error')
pylab.xlabel('pmra')
pylab.ylabel('pmdec')
pylab.xlim(-5.5,15.)
pylab.ylim(-10.,4.)

pylab.figure()
pylab.scatter(data.pmra[cut],data.pmdec[cut], c=data.pmdec_error[cut], vmin=0., vmax=.8, s=2.)
pylab.colorbar(label='pmdec error')
pylab.xlabel('pmra')
pylab.ylabel('pmdec')
pylab.xlim(-6.,10.)
pylab.ylim(-8.5,4.5)

pylab.figure()
pylab.scatter(data.pmra,data.pmdec, c=data.feh50, vmin=-3.5, vmax=.5, s=2.)
pylab.colorbar(label='[Fe/H]')
pylab.xlabel('pmra')
pylab.ylabel('pmdec')
pylab.xlim(-7.5,10.)
pylab.ylim(-10.,4.5)

