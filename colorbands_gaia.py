import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits

infile= '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match.fits'

reader= pyfits.open(infile)
data= reader[1].data

reader.close()

d= 1/data.PARALLAX
pylab.figure()
cut= (data.PARALLAX/data.PARALLAX_ERROR) > 1.
cut_1= d > 5.
pylab.scatter(data.WAVG_MAG_PSF_G_DERED[cut]           - data.WAVG_MAG_PSF_R_DERED[cut], data.WAVG_MAG_PSF_R_DERED[cut]          - data.WAVG_MAG_PSF_Z_DERED[cut], s=2.    )

pylab.xlabel('g-r (mag)')
pylab.ylabel('r-z (mag)')
pylab.xlim(0.2, 1.1)
pylab.ylim(0.0, 0.6)


pylab.figure()
mag= data.WAVG_MAG_PSF_G_DERED[cut] - data.WAVG_MAG_PSF_R_DERED[cut]
pylab.scatter(data.PMRA[cut],data.PMDEC[cut], c=mag, vmin=0., vmax=1.2, s=2.)
pylab.colorbar(label='g-r (mag')
pylab.xlabel("PMRA")
pylab.ylabel("PMDEC")


pylab.figure()
pylab.scatter(data.PMRA_ERROR, data.WAVG_MAG_PSF_G_DERED)
pylab.xlabel('PMRA Error')
pylab.ylabel('g (mag)')
pylab.xlim(0.,4.)
pylab.ylim(14.,23.)
