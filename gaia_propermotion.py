import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits

infile= '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match.fits'

reader= pyfits.open(infile)
data= reader[1].data
d=1/(data.PARALLAX)
reader.close()
parallax_ratio = data.PARALLAX / data.PARALLAX_ERROR
cut=(parallax_ratio > 1.)
pylab.figure()
pylab.scatter(data.PMRA[cut],data.PMDEC[cut],
c=d[cut], vmin=0., vmax=5.)
pylab.colorbar(label='distance')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.xlim(-20.,40.)
pylab.ylim(-30.,10.)



ratio_ra_error= (data.PMRA/data.PMRA_ERROR)
pylab.figure()
pylab.scatter(d,ratio_ra_error)
pylab.xlabel('distance')
pylab.ylabel('ratio PMRA/PMRA error')
c='r'





pylab.figure()
pylab.scatter(data.PMRA,data.PMRA_ERROR)
pylab.xlabel('PMRA')
pylab.ylabel('PMRA Error')
pylab.xlim(-30.,80.)
pylab.ylim(0.,4.)

parallax_ratio= data.PARALLAX / data.PARALLAX_ERROR

cut= (parallax_ratio > 1.)
d= data.PARALLAX
pylab.figure()
'''
pylab.scatter(data.WAVG_MAG_PSF_G_DERED[~cut] - data.WAVG_MAG_PSF_R_DERED[~cut], data.WAVG_MAG_PSF_G_DERED[~cut] - data.WAVG_MAG_PSF_I_DERED[~cut], c='0.5', s=2.)

pylab.scatter(data.WAVG_MAG_PSF_R_DERED[cut] - data.WAVG_MAG_PSF_G_DERED[cut], data.WAVG_MAG_PSF_R_DERED[cut] - data.WAVG_MAG_PSF_Z_DERED[cut], c= d[cut], vmin=0., vmax=5., s=2.)
pylab.colorbar(label='distance')
pylab.xlabel('r-g (mag)')
pylab.ylabel('g-i (mag)')
pylab.xlim(-2.,0.)
pylab.ylim(0.,2.)
'''


pylab.figure()
pylab.scatter(data.WAVG_MAG_PSF_G_DERED[cut] - data.WAVG_MAG_PSF_R_DERED[cut], data.WAVG_MAG_PSF_R_DERED[cut] - data.WAVG_MAG_PSF_Z_DERED[cut], c=d[cut], vmin=0., vmax=5., s=2.)
pylab.xlabel('g-r (mag)')
pylab.ylabel('r-z (mag)')
pylab.colorbar(label='distance')
pylab.xlim(-.5,2.)
pylab.ylim(-.5,2.5)    

'''
pylab.figure()
pylab.scatter(data.WAVG_MAG_PSF_G_DERED[cut] - data.WAVG_MAG_PSF_R_DERED[cut], data.WAVG_MAG_PSF_R_DERED[cut] - data.WAVG_MAG_PSF_I_DERED[cut], c=d[cut], vmin=0., vmax=5., s=2.)
pylab.xlabel('g-r (mag)')
pylab.ylabel('r-i (mag)')
pylab.colorbar(label='distance')


pylab.figure()
bins= np.arange(-99999.,99999.,10000.)
pylab.hist(data.RADIAL_VELOCITY,bins=bins)
pylab.xlabel('Radial Velocity')

pylab.figure()
pylab.scatter(data.RA,data.DEC,c=data.WAVG_MAG_PSF_R_DERED - data.WAVG_MAG_PSF_G_DERED, vmin=-3.5, vmax=5.5, s=2.)
pylab.colorbar(label='g-r(mag)')
pylab.xlabel('RA')
pylab.ylabel('DEC')




pylab.figure()
bins= np.arange(10.,25.,.5)
pylab.hist(data.WAVG_MAG_PSF_G_DERED, bins=bins)
bins= np.arange(10.,25.,.5)
pylab.xlabel('G band magnitude')
pylab.title('Proper motion no cut')


pylab.figure()
PM_ratio = (data.PMRA / data.PMRA_ERROR)
cut= (PM_ratio > 1.)
bins=np.arange(10.,25.,.5)

pylab.hist(data.WAVG_MAG_PSF_G_DERED[cut], bins=bins)
bins=np.arange(10.,25.,.5)
pylab.xlabel('G band magnitude')
pylab.title('Proper motion cut')
pylab.ylim(0.,175.)




pylab.figure()
c=PM_ratio
cut_color=((data.WAVG_MAG_PSF_G_DERED - data.WAVG_MAG_PSF_R_DERED) < 1.)
cut_1= (c > 1.)
cut_2= (c > 3.)
bins= np.arange (14.,22.,0.5)
pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_color], bins=bins, histtype='step', lw=2, color='black', label='all')

pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_1 & cut_color],
           bins=bins,
           histtype='step', lw=2, color='blue', label='PM SNR  > 1')
pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_2 & cut_color],
           bins=bins,
           histtype='step', lw=2, color='red', label='PM  SNR > 3')
pylab.xlabel('Magnitude')
pylab.title('cut data PM SNR       > 1')
pylab.ylim(0.,200.)
pylab.legend(loc='upper left')



pylab.figure()
pylab.scatter(data.PMRA,data.PMDEC, c=data.PMRA_ERROR, vmin=0., vmax=5., s=2.)
pylab.colorbar(label='PMRA Error')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.xlim(-10.,40.)
pylab.ylim(-25.,10.)

pylab.figure()
pylab.scatter(data.PMRA, data.PMDEC, c=data.PMDEC_ERROR, vmin=0., vmax=5., s=2.)
pylab.colorbar(label='PMDEC ERROR')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.xlim(-20.,60.)
pylab.ylim(-40.,15.)


pylab.figure()
mag= data.WAVG_MAG_PSF_G_DERED - data.WAVG_MAG_PSF_R_DERED

pylab.scatter(data.PMRA,data.PMDEC, c=mag, vmin=0., vmax=1., s=2.)
pylab.colorbar(label='g -r (mag)')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.xlim(-5.,20.)
pylab.ylim(-15.,5.)
'''
