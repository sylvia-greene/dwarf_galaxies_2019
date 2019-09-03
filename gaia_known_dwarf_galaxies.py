import numpy as np
import pylab
pylab.ion()
import astropy.io.fits as pyfits
import healpy 
import scipy
from scipy import stats 

import ugali.utils.projector

############################################################

def getData(ra, dec, nside=16):

    datadir = '/afs/hep.wisc.edu/home/kbechtol/public/des/projects/mw_substructure/y3_gold/data/y3_gold_gaia_match_bulk_nest16/'
    
    pix_center = healpy.ang2pix(nside, ra, dec, nest=True, lonlat=True)
    pix_neighbors = healpy.get_all_neighbours(nside, ra, dec, nest=True, lonlat=True)
    pix_all = np.concatenate([[pix_center], pix_neighbors])

    data_array = []
    for pix in pix_all:
        infile = '{}/y3_gold_gaia_match_bulk_nest16_{:05d}.fits'.format(datadir, pix)
        reader= pyfits.open(infile)
        data = reader[1].data
        reader.close()
        data_array.append(data)

    data = np.concatenate(data_array)
    return data
        
############################################################




def p_value(ra, dec, name=None):
    nside = 16

    data = getData(ra, dec)
    
   
    
    angsep = ugali.utils.projector.angsep(ra, dec, data['RA'], data['DEC']) # degrees
    cut_angsep = (angsep < 0.1)
    cut_PMRA = (data['PMRA'] > -50.)
    cut_PMDEC = (data['PMDEC'] > -50.)
 

    m1= data['PMRA'][cut_angsep & cut_PMRA & cut_PMDEC]
    m2= data['PMDEC'][cut_angsep & cut_PMRA & cut_PMDEC] 
    xmin = m1.min()
    xmax= m1.max()
    ymin = m2.min()
    ymax = m2.max()

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    z = np.reshape(kernel(positions).T, X.shape)

    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots()
    #ax.imshow(np.rot90(z), cmap=plt.cm.gist_earth_r,
          #extent=[xmin, xmax, ymin, ymax])
    #ax.plot(m1, m2, 'k.')
    #ax.scatter(X.flatten()[np.argmax(z)], Y.flatten()[np.argmax(z)], c='red', s=20, zorder=999)

    #ax.set_title("kernel density for %s" %name)
    #ax.set_xlabel("PMRA")
    #ax.set_ylabel("PMDEC")

    #plt.show()


              
    #x= np.linspace(-1., 3., num=50)
    #t = .8845*x**3 - 1.263*x**2 +.7366*x - 0.02765
    #pylab.figure()
    #pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep]- data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_I_CORRECTED'][~cut_angsep], c='.75', s= 1, label='All stars')

    #pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_I_CORRECTED'][cut_angsep], c = data['PMRA'][cut_angsep], s = 10, label = 'Dwarf candidates')
    #pylab.plot(x, t, c='red', zorder=999)
    
    
    #pylab.legend()
    #pylab.colorbar(label='PMRA')
    #pylab.xlabel('g-r (mag)')
    #pylab.ylabel('r -i(mag)')
    #pylab.title('color magnitude for %s' %name)

 

    #pylab.figure()
    #pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep], data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep], c='.75', s=1, label='All stars')
    #pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep], data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep], c= data['PMRA'][cut_angsep], s=10, zorder=999, label='Dwarf Candidates')
    #pylab.legend()
  
    #pylab.colorbar(label='PMRA')
    #pylab.xlabel('g-r (mag)')
    #pylab.ylabel('g')
    #pylab.title('%s' %name)
    
    

    highest_density_x = X.flatten()[np.argmax(z)]
    highest_density_y = Y.flatten()[np.argmax(z)]

    cond = np.where(
    np.sqrt((m1 - highest_density_x)**2
            +(m2 - highest_density_y)**2)
     <= 2.5)
    good_x = m1[cond]
    good_y = m2[cond]

    cond = np.where(
    np.sqrt((data['PMRA'] - highest_density_x)**2
           +(data['PMDEC'] - highest_density_y)**2)
    <= 2.5)
    new_x = data['PMRA'][cond]
    new_y = data['PMDEC'][cond]

    black_points_in_circle = len(good_x)
    total_black_points = len(m1)
    red_points_in_circle = len(new_x)
    total_points = len(data['PMRA'])

    n_observed = black_points_in_circle
    n_predicted = (total_black_points * red_points_in_circle)/total_points

    pvalue = scipy.stats.poisson.sf(mu = n_predicted, k= n_observed)

    #print(pvalue)
    #if name is not None:
        #print("pvalue for %s" % name)

    return pvalue




    

#random positons

n = 100
ra_array = np.random.uniform(0., 360., size=n)
dec_array = np.degrees(np.arcsin(np.random.uniform(-1., 1., size=n)))
p_array = np.empty(n)

for i in range(n):
     #random_ra = np.random.uniform(0,360)
     #random_dec = np.random.random()
     #name = ('ra %s, %s' %(ra_array, dec_array) )
     #try:
     #    p_value(random_ra, random_dec)
     try:
         p_array[i] = p_value(ra_array[i], dec_array[i])
     except FileNotFoundError:
         p_array[i] = -9999.
         continue

cut = (p_array > 0.)
     
pylab.figure()  
pylab.scatter(ra_array[cut], dec_array[cut], c=p_array[cut])
pylab.xlabel('RA')
pylab.ylabel('DEC')
pylab.colorbar(label='p-value')
    

                   

     
#ret II
#name= 'ret II'
#p_value(53.92, -54.05)


#ra = 53.92
#dec = -54.05
#data = getData(ra, dec)
#angsep = ugali.utils.projector.angsep(ra, dec, data['RA'], data['DEC']) # degrees

#cut_angsep = (angsep < 0.1)

#x = data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep]
#y=  data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_I_CORRECTED'][cut_angsep]

#z= np.polyfit(x, y, 3)
#f= np.poly1d(z)

#xmin = x.min()
#xmax = x.max()
#x_new = np.linspace([xmin], [xmax], 10.)
#y_new = f(x_new)

#print(f)

#tuc II
name = 'Tucana II'
p_value(343.06, -58.57)


#tuc III
name = 'Tucana III'
p_value(359.15, -59.6)

#hor I
name = 'Horologium I'
p_value(43.87, -54.11)
'''
#Grus I
name = 'Grus I'
p_value(344.1765, -50.1633)

#Grus II
name= 'Grus II'
p_value(331.02, -46.44)

#Pictor I
name= 'Pictor I'
p_value(70.9475, -50.283)

#Tucana V
name = 'Tucana V'
p_value(354.35, -63.27)

#Indus II
name = 'Indus II'
p_value(309.72, -46.16)

#Phoenix II
nside = 16
ra = 354.9975
dec = -54.4060
data = getData(ra, dec)
name = 'Phoenix II'
p_value(354.9975, -54.4060)

#Cetus II
name = 'Cetus II'
p_value(19.47, -17.42)

#Columba I
name = 'Columba I'
p_value(82.857, -28.0435)

#Eridanus II
name = 'Eridanus II'
p_value(56.0838, -43.5338)

#Eridanus III
name = 'Eridanus III'
p_value(35.6888, -52.2847)

#Horologium II
name = 'Horologium II'
p_value(49.1338, -50.0181)

#Tucana IV
name = 'Tucana IV'
p_value(.73, -60.850)

####################
#color magnitude for Ret II

pylab.figure()
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][~cut_angsep], c='.75', s=1., label="other stars")
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][cut_angsep], c=data['PMRA'][cut_angsep], s=30., label="dwarf candidates")
pylab.colorbar(label="PMRA")
pylab.legend()
pylab.title('Color magnitude for Ret II')
pylab.xlabel('g -r (mag)')
pylab.ylabel('r - z (mag)')

x = data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep]
y=  data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][cut_angsep]

z= np.polyfit(x, y, 5)
f= np.poly1d(z)

xmin = x.min()
xmax = x.max()
x_new = np.linspace([xmin], [xmax], 10.)
y_new = f(x_new)

print(f)

pylab.plot(x_new, y_new)
pylab.xlim(-1.,2.5)
pylab.ylim(-.5, 3)

x = np.linspace(-1., 3., num=30)
t = (-0.07518*x**5 + 0.8382*x**4 - 0.8168*x**3 - 0.1145*x**2 + 0.8078*x - 0.10030)

pylab.figure()
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][~cut_angsep], c='.75', s=1., label="other stars")
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][cut_angsep], c=data['PMRA'][cut_angsep], s=30., label="dwarf candidates")
pylab.colorbar(label="PMRA")
pylab.title('Color magnitude Tuc III')
pylab.plot(x, t, markersize=20, color='blue', zorder=999)
pylab.xlabel('g -r (mag)')
pylab.ylabel('r -z (mag)')
pylab.xlim(-.5,2.)
pylab.ylim(-1.,3)
pylab.legend()


#########################


nside= 16
ra = 354.65054
dec= -54.405307
data = getData(ra, dec)

angsep = ugali.utils.projector.angsep(ra, dec, data['RA'], data['DEC']) # degrees

cut_angsep = (angsep < 0.1)
cut_PMRA = (data['PMRA'] > -50.)
cut_PMDEC = (data['PMDEC'] > -50.)

x = np.linspace(-1., 3., num=30)
t = (-0.07518*x**5 + 0.8382*x**4 - 0.8168*x**3 - 0.1145*x**2 + 0.8078*x - 0.10030)

pylab.figure()
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][~cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][~cut_angsep], c='.75', s=1., label="other stars")
pylab.scatter(data['SOF_PSF_MAG_G_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep], data['SOF_PSF_MAG_R_CORRECTED'][cut_angsep] - data['SOF_PSF_MAG_Z_CORRECTED'][cut_angsep], c=data['PMRA'][cut_angsep], s=30., label="dwarf candidates")
#pylab.plot(x,t, c='black')
pylab.colorbar(label="PMRA")
pylab.title('Color magnitude Phe II')
pylab.plot(x, t, markersize=20, color='blue', zorder=999)
pylab.xlabel('g -r (mag)')
pylab.ylabel('r -z (mag)')
pylab.xlim(-.5,2.)
pylab.ylim(-1.,3)
pylab.legend()
'''










































