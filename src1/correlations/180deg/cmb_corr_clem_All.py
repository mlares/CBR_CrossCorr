import numpy as np
import matplotlib.pyplot as plt

import healpy as hp

import matplotlib as mp
mp.use('Agg')

import copy
import random

import time

from joblib import Parallel, delayed

import pickle

plt.close()

# LECTURA MAPAS TEMPERATURA CMB_____________________

filename='../../Data/cmb_maps/COM_CMB_IQU-smica_2048_R3.00_full.fits'
hp_data, hp_header = hp.read_map(filename, h=True,field=(0,3),nest=True)#,

# CONSTRUYO LAS MÁSCARAS ----------------------------------------------------------
hp_data_sel=copy.copy(hp_data[0])
hp_mask=copy.copy(hp_data[1])

#........................................................

# CALCULO LA CORRELACION USANDO PARES RANDOM:


#........................................................

#PARAMETROS:

nside=2048

first=0.03/57.2958  #<--- primer bin, de grados a radianes desde 100 segundos

last=180./57.2958   #<--- ultimo bin, de grados a radianes

nbins=3600

#npairs=3000000        #<--- numero de pares usados para estimar la correlacion en cada parche
npairs=500000

#nran=500
nran=10

listpix=np.array(list(range(hp.nside2npix(nside))))
#........................................................


def PCorrelation_highres(rseed,npairs,nside,hp_data_sel,hp_mask,listpix,first,last,nbins):
    
    random.seed((rseed+42)*3)

    lbin=(last-first)/float(nbins)

    theta_corr=np.arange(start=first,stop=last,step=lbin) 

    DTT=np.zeros(len(theta_corr))
    DnTT=np.zeros(len(theta_corr))

    for j in range(npairs):

        pix1=listpix[int(np.random.random()*(len(listpix)-1))]
        pix2=listpix[int(np.random.random()*(len(listpix)-1))]

        if((hp_mask[pix1]==1)and(hp_mask[pix2]==1)):

            v1=hp.pix2vec(nside,pix1,nest=True)
            v2=hp.pix2vec(nside,pix2,nest=True)

            ad1=hp.rotator.angdist(v1,v2)

            tt1=hp_data_sel[pix1]*hp_data_sel[pix2]*(10**(12))

            binang1=int(((ad1-first)/(last-first))*(nbins-1))

            if((binang1>-1) and (binang1<nbins)):

                DTT[binang1]=DTT[binang1]+tt1
                DnTT[binang1]=DnTT[binang1]+1
    
    ang_corr1=(DTT/DnTT)

    aa=np.array([theta_corr,ang_corr1,DTT,DnTT])

    return(aa)


#........................................................

Tstart=time.time()
aalist=Parallel(n_jobs=56,verbose=49, backend="threading")(delayed(PCorrelation_highres)
(rseed,npairs,nside,hp_data_sel,hp_mask,listpix,first,last,nbins)for rseed in range(nran))
Tend=time.time()

print(Tend-Tstart,'sec','pairs',npairs)

#........................................................

xx=aalist[0][0]
yyc=xx*0.
yyn=xx*0.

for i in range(nbins-1):
  for j in range(nran-1):

    yyc[i]=yyc[i]+aalist[j][2][i]
    yyn[i]=yyn[i]+aalist[j][3][i]

yy=yyc/yyn

#........................................................

jack=yy*0.
call=yy*0.
nall=yy*0.


for i in range(nbins-1):

    for j in range(nran-1):

      call[i]=aalist[j][2][i]+call[i]
      nall[i]=aalist[j][3][i]+nall[i]


for i in range(nbins-1):

  for j in range(nran-1):

    if(np.isnan(aalist[j][1][i])==False ):

      jack[i]=jack[i]+(((call[i]-aalist[j][2][i])/(nall[i]-aalist[j][3][i]))-yy[i])**2

jack=np.sqrt(jack*((nran-1)/nran))



#........................................................

plt.close()
for i in range(nran-1):

  plt.plot(aalist[i][0]*57,aalist[i][1],'--',color='black',alpha=0.2)
  

plt.axhline(y=0., color='b',alpha=0.1, linestyle='-')
plt.fill_between(xx*57, yy-jack, yy+jack,alpha=0.2,edgecolor='darksalmon', facecolor='lightsalmon')

plt.plot(xx*57,yy,color='darksalmon')

plt.savefig('/mnt/is0/heliana/corr_All.pdf')

#........................................................

fname = '/mnt/is0/heliana/corr_All.txt'
with open( fname, "wb" ) as f:
    pickle.dump(aalist, f)

#........................................................


