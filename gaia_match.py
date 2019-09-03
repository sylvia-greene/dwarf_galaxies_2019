import numpy as np
import astropy.io.fits as pyfits
import pylab

pylab.ion()

infile = '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match.fits'

print('Hello!')
print(infile)

reader = pyfits.open(infile)
data = reader[1].data
reader.close()

# To get the list of column names
#data.names
# Number of rows
#len(data)
# Access a particular column
# data['GAIA_RA']

pylab.figure()
pylab.scatter(data['GAIA_RA'], data['GAIA_DEC'], 
              c=data['PARALLAX'],
              vmin=0.)
pylab.colorbar(label='PARALLAX')
pylab.xlabel('GAIA_RA')
pylab.ylabel('GAIA_DEC')
# Reverses the RA axis to follow astro convention
pylab.xlim(pylab.xlim()[::-1])

bins = np.arange(14., 23., 0.5)
pylab.figure()
pylab.hist(data['WAVG_MAG_PSF_R_DERED'], bins=bins)
pylab.xlabel('r (mag)')
pylab.ylabel('Counts')

