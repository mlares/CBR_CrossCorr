#----------------LECTURA DEL CATÁLOGO---------------------------- 

import numpy as np
import matplotlib.pyplot as plt
plt.close()

from io import StringIO
s = StringIO("1,1.3,abcde")
D = np.genfromtxt('../Data/2mass/2mrs_1175_done.dat',unpack=True,names=True,skip_header=9,delimiter='',
dtype=None)

D = D[D["ba"]<1.] #hay solo un valor choto y es 9.999


#SELECCION DE GALAXIAS

#splist=["1","2","3","4","5","6","7","8"] #<------ ESPIRALES opciones para el primer caracter de la columna
splist=["4","5","6"] #<------Sc opciones para el primer caracter de la columna
#eplist=["5","6","7"] #<------ ELIPTICAS opciones para el primer caracter de la columna

spty=np.zeros(len(D))

for i in range(len(D)):
#for i in range(len(D)-len(D)+5):
      
      sp=0

      #for j in range(len(eplist)):
      for j in range(len(splist)):


            #if((str(D['type'][i])[2]=='-')&(str(D['type'][i])[3]==eplist[j])):
            if((str(D['type'][i])[2]==splist[j])):

                  #print('yes!',str(D['type'][i])[2],i)
            
                  sp=1

      spty[i]=sp


Dsp=D[spty==1]

Dsp=Dsp[(Dsp['v']/300000.>0.001)&(Dsp['v']/300000.<0.01)]

#SELECCION DE GALAXIAS ESPIRALES FACE ON ¿?
DspF=Dsp[Dsp["ba"]>0.5]


#SELECCION DE GALAXIAS ESPIRALES EDGE ON ¿?
DspE=Dsp[Dsp["ba"]<0.5]

Dsel=DspE

#mapa de chai
#filename='/home/heliana/CMB_analysis/data/caimap/CMB_Lensing/lensmap512_10arcmin_y2.fits'  



#----------------CORRELACIONES CON CMB----------------------------

import healpy as hp


import cmb_functions as fnc

import imp
imp.reload(fnc)

import matplotlib as mp
mp.use('Agg')

# LECTURA MAPA_____________________

filename='../Data/cmb_maps/COM_CMB_IQU-smica_2048_R3.00_full.fits'
hp_data, hp_header = hp.read_map(filename, h=True,field=(0,3))#,

filename='../Data/cmb_maps/COM_CMB_IQU-nilc_2048_R3.00_full.fits'
hp_data1, hp_header1 = hp.read_map(filename, h=True,field=(0,3))#,

filename='../Data/cmb_maps/COM_CMB_IQU-sevem_2048_R3.01_full.fits'
hp_data2, hp_header2 = hp.read_map(filename, h=True,field=(0,3))#,

filename='../Data/cmb_maps/COM_CMB_IQU-commander_2048_R3.00_full.fits'
hp_data3, hp_header3 = hp.read_map(filename, h=True,field=(0,3))#,

#filename='../Data/cmb_maps/COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits'
#hp_data, hp_header = hp.read_map(filename)#,

#filename='../Data/cmb_maps/lensmap512_10arcmin_y2.fits'
#hp_data_sel, hp_header = hp.read_map(filename, h=True)#,field=(0,3))#, nest=None)


#filename='../Data/cmb_maps/lensmask512_10arcmin_y2.fits'
#hp_mask, hp_maskheader = hp.read_map(filename, h=True)#,field=(0,3))#, nest=None)


#0=smica
#1=nilc
#2=sevem
#3=commander


# PARAMETERS_______________________

nside=2048
#nside=512

fac=6

rprof=np.arange(start=0.2,stop=40.8,step=0.6) 

vec=hp.ang2vec(np.pi/2.-Dsel['b']*np.pi/180.,Dsel['l']*np.pi/180.,lonlat=False)
#nran=150

# PERFILES_________________________ 

hp_data_sel=hp_data[0]
hp_mask=hp_data[1]

Kr,serr=fnc.profdata(nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask)

hp_data_sel=hp_data1[0]
hp_mask=hp_data1[1]

Kr1,serr1=fnc.profdata(nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask)


hp_data_sel=hp_data2[0]
hp_mask=hp_data2[1]

Kr2,serr2=fnc.profdata(nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask)


hp_data_sel=hp_data3[0]
hp_mask=hp_data3[1]

Kr3,serr3=fnc.profdata(nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask)


# PERFILES RANDOM___________________ 

#Kr_ran=fnc.profran(nran,nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask)


# PLOT =) ____________________________

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

#fig = plt.figure(1)
#ax = fig.add_subplot(111)

#for i in range(nran):

#        plt.plot(rprof, Kr_ran[i,:], '--',color='black',alpha=0.2)

#plt.axhline(y=0., color='b',alpha=0.1, linestyle='-')
#plt.fill_between(rprof, Kr-serr, Kr+serr,alpha=0.2,edgecolor='darksalmon', facecolor='lightsalmon')
#plt.plot(rprof, Kr, '-',color='maroon',alpha=0.4)
#plt.plot(rprof, Kr, 'o',color='maroon')

#plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=13)
#plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=13)
#plt.setp(ax.tick_params(direction='in'))

#plt.xlabel('r[rad]',size=15)
#plt.ylabel('Mean~~$\Delta T~[10^6]$',size=15)
plt.axhline(y=0., color='b',alpha=0.1, linestyle='-')
plt.fill_between(rprof, Kr-serr, Kr+serr,alpha=0.2,edgecolor='lightgray', facecolor='lightgray')
plt.fill_between(rprof, Kr1-serr1, Kr1+serr1,alpha=0.2,edgecolor='lightgray', facecolor='lightgray')
plt.fill_between(rprof, Kr2-serr2, Kr2+serr2,alpha=0.2,edgecolor='lightgray', facecolor='lightgray')
plt.fill_between(rprof, Kr3-serr3, Kr3+serr3,alpha=0.2,edgecolor='lightgray', facecolor='lightgray')
plt.plot(rprof, Kr, '-',label='smica')
plt.plot(rprof, Kr1, '-',label='nilc')
plt.plot(rprof, Kr2, '-',label='sevem')
plt.plot(rprof, Kr3, '-',label='commander')
plt.legend()
#0=smica
#1=nilc
#2=sevem
#3=commander


#plt.savefig('Temp_profiles_2048_L05_maps_r3n.pdf')
plt.savefig('Temp_profiles_2048_L05_maps_r40n_Sc.pdf')

#___________________________________

