import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits
import pylab


infile= '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match.fits'

reader= pyfits.open(infile)
data= reader[1].data

reader.close()
"""
pylab.figure()
bins= np.arange(14.,22.,0.5)
pylab.hist(data.WAVG_MAG_PSF_R_DERED,bins=bins)
pylab.title('All data')
pylab.xlabel('Magnitude')
"""

pylab.figure()
parallax_ratio=(data.PARALLAX/data.PARALLAX_ERROR)
c=parallax_ratio
#cut_color = np.tile(True, len(data))
cut_color = ((data.WAVG_MAG_PSF_G_DERED - data.WAVG_MAG_PSF_R_DERED) < 1.)
cut_1 = (c > 1.)
cut_3 = (c > 3.)
bins=np.arange(14.,22.,0.5)
pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_color],
           bins=bins,
           histtype='step', lw=2, color='black', label='All')
pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_1 & cut_color],
           bins=bins,
           histtype='step', lw=2, color='blue', label='Parallax SNR > 1')
pylab.hist(data.WAVG_MAG_PSF_R_DERED[cut_3 & cut_color],
           bins=bins,
           histtype='step', lw=2, color='red', label='Parallax SNR > 3')
pylab.xlabel('Magnitude')
pylab.title('cut data parallax ratio > 1')
pylab.ylim(0.,200.)
pylab.legend(loc='upper left')

