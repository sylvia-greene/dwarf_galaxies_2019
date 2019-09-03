import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits
import healpy 

nside = 16
ra= 53.92
dec= -54.05

pix = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)

datadir = '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match_bulk_nest16/'

infile = '{}/y3_gold_gaia_match_bulk_nest16_{:05d}.fits'.format(datadir, pix)

reader= pyfits.open(infile)
data = reader[1].data

d= 1/data.PARALLAX
cut= (d > 8.) 

pylab.figure()
pylab.scatter(data.PMRA[cut],data.PMDEC[cut], c= data.SOF_PSF_MAG_G_CORRECTED[cut]      - data.SOF_PSF_MAG_R_CORRECTED[cut]  , vmin=0., vmax=2., s=10, zorder=1)
pylab.colorbar(label='g-r(mag)')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra = 53.92 degrees dec=-54.05, location of known dwarf galaxy, with distance cut')

pylab.scatter(data.PMRA[~cut],data.PMDEC[~cut], c='0.75', s=1., zorder=0)
pylab.xlim(0.,10.)
pylab.ylim(-8.,6.)

pylab.figure()
pylab.scatter(data.PMRA,data.PMDEC, c= d     , vmin=0., vmax=15., s=2.)
pylab.colorbar(label='Distance')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra = 53.92 degrees dec=-54.05, location of known dwarf galaxy')

pylab.figure()
pylab.scatter(data.RA[cut],data.DEC[cut], c= data.PMRA[cut], vmin=0., vmax=8., s=2.)
pylab.colorbar(label="PMRA")
pylab.xlabel('ra')
pylab.ylabel('dec')
pylab.title('tile ra = 53.92 degrees dec= -54.05 degrees, known dwarf galaxy location')

pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED - data.SOF_PSF_MAG_R_CORRECTED, data.SOF_PSF_MAG_R_CORRECTED - data.SOF_PSF_MAG_Z_CORRECTED, c = data.PMRA, vmin=0., vmax=15., s=2.)
pylab.colorbar(label='PMRA')
pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z(mag)")
pylab.title("tile ra=53.92 dec=-54.05, known dwarf galxaxy")

pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED      - data.SOF_PSF_MAG_R_CORRECTED     , data.SOF_PSF_MAG_R_CORRECTED -data.SOF_PSF_MAG_Z_CORRECTED , c=d, vmin= 0, vmax= 15., s=2)
pylab.colorbar(label='distance')
pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z  ag)")
pylab.title("tile ra=53.92 dec=-54.05, known dwarf galxaxy")

nside = 16
ra= 20.
dec= -50.

pix = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)

datadir = '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match_bulk_nest16/'

infile = '{}/y3_gold_gaia_match_bulk_nest16_{:05d}.fits'.format(datadir, pix)

reader= pyfits.open(infile)
data = reader[1].data
d= 1/data.PARALLAX
cut= (d > 8.) 

pylab.figure()
pylab.scatter(data.RA[cut],data.DEC[cut], c= data.PMRA[cut], vmin=0., vmax=8., s=2.)
pylab.colorbar(label="PMRA")
pylab.xlabel('ra')
pylab.ylabel('dec')
pylab.title('Tile ra = 20 degrees dec = -50 degrees')
d= 1/data.PARALLAX
cut= (d > 8.)

pylab.figure()
pylab.scatter(data.PMRA[cut],data.PMDEC[cut], c= data.SOF_PSF_MAG_G_CORRECTED[cut] - data.SOF_PSF_MAG_R_CORRECTED[cut], vmin=0., vmax=2., s=10., zorder=1)
pylab.colorbar(label='g-r(mag)')
pylab.scatter(data.PMRA[~cut],data.PMDEC[~cut], c='0.75',s=1., zorder=0)
pylab.xlim(0.,10.)
pylab.ylim(-8.,6.)


pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra= 20 degrees dec=-50degrees with cut')


d= 1/data.PARALLAX

pylab.figure()
pylab.scatter(data.PMRA,data.PMDEC, c= d     , vmin=0., vmax=15., s=2.)
pylab.colorbar(label='distance')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra= 20 degrees dec=-50degrees')


pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED - data.SOF_PSF_MAG_R_CORRECTED, data.SOF_PSF_MAG_R_CORRECTED - data.SOF_PSF_MAG_Z_CORRECTED, c = data.PMRA, vmin=0., vmax=15., s=2.)
pylab.colorbar(label='PMRA')

pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z(mag)")
pylab.title('ra=20 degrees dec=-50 degrees')

pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED   - data.SOF_PSF_MAG_R_CORRECTED     , data.SOF_PSF_MAG_R_CORRECTED     -data.SOF_PSF_MAG_Z_CORRECTED     , c=d     , vmin=0., vmax= 15., s=2)
pylab.colorbar(label='distance')
pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z(mag)")
pylab.title("tile ra=20 dec= -50")




nside = 16
ra= 350.

dec= -50.

pix = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)

datadir = '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match_bulk_nest16/'

infile = '{}/y3_gold_gaia_match_bulk_nest16_{:05d}.fits'.format(datadir, pix)

reader= pyfits.open(infile)
data = reader[1].data

d= 1/data.PARALLAX
cut= (d > 8.) 

pylab.figure()
pylab.scatter(data.PMRA[cut],data.PMDEC[cut], c= data.SOF_PSF_MAG_G_CORRECTED[cut]      - data.SOF_PSF_MAG_R_CORRECTED[cut]  , vmin=0., vmax=2., s=10., zorder=1)
pylab.colorbar(label='g-r(mag)')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra = 350   degrees dec=-50   ,')
pylab.scatter(data.PMRA[~cut],data.PMDEC[~cut], c='0.75', s=1., zorder=0)
pylab.xlim(0.,10.)
pylab.ylim(-8.,6.)


pylab.figure()
pylab.scatter(data.PMRA,data.PMDEC, c= d     , vmin=0., vmax=15., s=2.)
pylab.colorbar(label='Distance')
pylab.xlabel('PMRA')
pylab.ylabel('PMDEC')
pylab.title('ra = 350   degrees dec=-50')
pylab.figure()
pylab.scatter(data.RA,data.DEC)
pylab.xlabel('ra')
pylab.ylabel('dec')
pylab.title('tile ra = 350   degrees dec= -50    n')

pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED - data.SOF_PSF_MAG_R_CORRECTED, data.SOF_PSF_MAG_R_CORRECTED - data.SOF_PSF_MAG_Z_CORRECTED, c = data.PMRA, vmin=0., vmax=15., s=2.)
pylab.colorbar(label='PMRA')
pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z(mag)")
pylab.title("tile ra=350   dec=-50   , ")


pylab.figure()
pylab.scatter(data.SOF_PSF_MAG_G_CORRECTED - data.SOF_PSF_MAG_R_CORRECTED     , data.SOF_PSF_MAG_R_CORRECTED     -data.SOF_PSF_MAG_Z_CORRECTED     , c=d     , vmin=0., vmax= 15., s=2)
pylab.colorbar(label='distance')
pylab.xlabel("g-r(mag)")
pylab.ylabel("r-z(mag)")
pylab.title("tile ra=30 dec= -50")

