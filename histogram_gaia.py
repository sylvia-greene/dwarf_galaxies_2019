import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits

save = False

infile= '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match.fits'

reader= pyfits.open(infile)
data= reader[1].data

reader.close()

d= 1/(data.PARALLAX)

pylab.figure()
bins= np.arange(14.,22.,0.5)
pylab.hist(data['WAVG_MAG_PSF_I_DERED'],bins=bins)
pylab.xlabel('MAG_I')


pylab.figure()
bins= np.arange(14.,22.,0.5)
pylab.hist(data['WAVG_MAG_PSF_R_DERED'],bins=bins)
pylab.xlabel('WAVG_MAG_PSF_R_DERED')

pylab.figure()
bins=np.arange(14.,22.,0.5)
pylab.hist(data['WAVG_MAG_PSF_Y_DERED'],bins=bins)
pylab.xlabel('MAG_Y')           

pylab.figure()
bins=np.arange(14.,22.,0.5)
pylab.hist(data['WAVG_MAG_PSF_G_DERED'],bins=bins)
pylab.xlabel('MAG_G')

pylab.figure()
bins=np.arange(14.,22.,0.5)
pylab.hist(data['WAVG_MAG_PSF_Z_DERED'],bins=bins)
pylab.xlabel('MAG_Z')


pylab.figure()
# I think that I did this histogram wrong
pylab.hist(d,bins=bins)
bins= np.arange(0.,15.,.05)
pylab.xlabel('distance')
pylab.ylabel('Frequency')

pylab.figure()
pylab.scatter(d,data.PARALLAX)
pylab.xlabel('Distance')
pylab.ylabel('parallax')
pylab.xlim(0.,8.)
pylab.ylim(0.,10.)

pylab.figure()
pylab.scatter(data.PARALLAX,data.PARALLAX_ERROR)
pylab.xlabel('PARALLAX')
pylab.ylabel('PARALLAX_ERROR')
pylab.xlim(-2.,8.)
pylab.ylim(0.,2.)



pylab.figure()
pylab.scatter(data.WAVG_MAG_PSF_I_DERED,data.PARALLAX_ERROR)
pylab.xlabel('Magnitude')
pylab.ylabel('Parallax Error')
pylab.xlim(14.,21.)
pylab.ylim(0.,3.)

pylab.figure()
pylab.scatter(d,data.WAVG_MAG_PSF_I_DERED)
pylab.xlabel('Distance')
pylab.ylabel('Magnitude')
pylab.xlim(-30.,30.)
pylab.ylim(15.,21.)


r=data.PARALLAX/data.PARALLAX_ERROR
pylab.figure()
pylab.scatter(r,d)
pylab.xlabel('Ratio Parallax/parallax error')

pylab.ylabel('distance')
pylab.xlim(-1.,14.)
pylab.ylim(-10.,50.)

pylab.figure()
pylab.scatter(data.PMDEC,d)
pylab.xlabel('PMDEC')
pylab.ylabel('Distance')

pylab.figure()
pylab.scatter(data.PMRA,d)
pylab.xlabel('PMRA')
pylab.ylabel('distance')

pylab.figure()
bins=np.arange(-10.,10.,1.)
pylab.hist(data['PMRA'],bins=bins)
pylab.xlabel('Proper motion right ascension')


pylab.figure()
pylab.scatter(data.PARALLAX,data.PARALLAX_ERROR)
pylab.xlabel('PARALLAX')
pylab.ylabel('PARALLAX_ERROR')














