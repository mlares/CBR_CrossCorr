# --------------------- IMPORT MODULES ------------------------------------ : #

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import healpy as hp
import time

# -------------- USEFUL DEFINITIONS AND PARAMETERS ------------------------- : #

c_luz = 2.99792458e+5			# Velocidad de la luz en km x s^{-1}.

arcsec_to_rad = np.pi/( 180.0*3600.0 )	# Convertion factor between [arcseconds] and [radians].
arcmin_to_rad = np.pi/( 180.0*60.0 )	# Convertion factor between [arcminutes] and [radians].
deg_to_rad = np.pi/( 180.0 )		# Convertion factor between [degree] and [radians].
rad_to_arcsec = ( 180.0*3600.0 )/np.pi  # Convertion factor between [radians] and [arcseconds].
rad_to_arcmin = ( 180.0*60.0 )/np.pi    # Convertion factor between [radians] and [arcminutes].
rad_to_deg    = ( 180.0 )/np.pi         # Convertion factor between [radians] and [degree].
K_to_muK = 10**6 			# Convertion factor between [Kelvin] and [microKelvin].

Nside = 2048				# Nisde resolution of the CMB maps.
Samp_choice = 0.03 	# 0.02 , 0.04, etc.		# This parameter chooses the sample
Rext_cut = True			# If True it cut galaxies with 10**(r_ext) > 120.0 arcsec.
#colormap = 'viridis' 
#colormap = 'plasma' 
colormap = 'inferno' 
#colormap = 'magma'
#colormap = 'cividis'
x_Rext = 50.0
Tmax = 0.0005*K_to_muK
Tmin = -0.0005*K_to_muK

print('###############################################################')
print('# USEFUL DEFINITIONS AND PARAMETERS                           #')
print('###############################################################')
print('')
print('Speed of light in km x s^{-1} =', c_luz)
print('Nisde resolution of the CMB maps = ', Nside)
print('Color pallete CMB maps = ', colormap)
print('')
print('###############################################################')
print('')
print('')
print('###############################################################')
print('# GALAXY SAMPLES IN USE			                             #')
print('###############################################################')
print('')
print('Upper limit in z =', Samp_choice)
print('They have cut in galaxies larger that r_ext = 120[arsec]?', Rext_cut)
print('')
print('###############################################################')



# ---- CMB Maps data files ------------------------------: 

#fileCMB = '/home/ezequiel/Documentos/Planck+Data+Product+Release_2019/'
fileCMB = '/home/zboero/Documentos/Planck+Data+Product+Release_2019/'

SMICA=(fileCMB+'CMB_Maps/COM_CMB_IQU-smica_2048_R3.00_full.fits')
SEVEM=(fileCMB+'CMB_Maps/COM_CMB_IQU-sevem_2048_R3.01_full.fits')
NILC=(fileCMB+'CMB_Maps/COM_CMB_IQU-nilc_2048_R3.00_full.fits')
Commander=(fileCMB+'CMB_Maps/COM_CMB_IQU-commander_2048_R3.00_full.fits')

CMBmap =[SMICA, SEVEM, NILC, Commander]
idMap = ['SMICA', 'SEVEM', 'NILC', 'Commander']	# ID's for the maps

fileCMB_Maps = '/home/zboero/Projects/CMB/Maps/'

# ---- 2MASS galaxy files -------------------------------: 

#fileProject = '/home/ezequiel/Projects/CMB/Anisotropy_signal_around_ScGalaxies/'	#Local PC path
fileProject = '/home/zboero/Projects/CMB/Anisotropy_signal_around_ScGalaxies/'		#Path in cluster
file_TwoMASS = fileProject+'data/2MASS/2mrs_1175_done.dat'

Columns = ['l', 'b', 'r_iso', 'r_ext', 'b/a', 'type', 'v']
# We use Gal Coord (l,b) since it is appropriated at the moment to work with CMB maps.
hd= 9					# Commented lines of header in the 2MASS catalog file

id_Sample = ['Sa_sample', 'Sb_sample', 'Sc_sample']
# -------------------------- FUNCTIONS  -------------------------------------: #

def v2z(x):				# It makes the conversion Velocity in Km/s to redshift.
	return x/c_luz

def z2v(x):				# It makes the conversion redshift to Velocity in Km/s.
    return x*c_luz

def Rd_map(filename):			# It reads the CMB maps.
    CMB_T, header = hp.read_map( filename, nest=False, h=True, field=(0,3) )
    T_data = CMB_T[0]           	# Temperature (numpy) array 
    T_mask = CMB_T[1]               	# Mask (numpy) array with 0's and 1's.
    return( T_data, T_mask )

def Rd_cat(filename, hd, Col):		# It reads the catalog of galaxies (2MASS) as a DataFrame
    Gxs = pd.read_table(filename, sep="\s+", header=hd, usecols= Col)


def CMB_MwPlot(T_data, Tmin, Tmax, idMap, id_Sample):
	plt.close()
	plt.figure()
    
	hp.mollview( T_data*K_to_muK, coord='G',unit='Temperature [$\mu$K]', xsize=800,\
				title='CMB Temperature Map ('+str(idMap)+') + '+str(id_Sample)+' mask',\
				cbar=True, cmap=colormap, min=Tmin, max=Tmax, badcolor='black',bgcolor='white')	
	plt.show()
	plt.savefig(fileCMB_Maps+'graficos/'+str(idMap)+'_Tmask_'+str(id_Sample)+'.png')
	plt.close()

				
def CMB_Hist(T_data, T_mask, idMap):
	T_data = T_data[T_mask == 1]
	sigma = np.std(T_data)*K_to_muK	
	
	sigmashort=float('{0:.5f}'.format(sigma))
	plt.close()
	plt.figure()
	plt.hist( T_data*K_to_muK, range=(-6*sigma, 6*sigma), bins=5000, label='$\sigma=$'+str(sigmashort)+'$[\mu K]$' )
	plt.title('CMB Temperature distribution ('+str(idMap)+')', loc='center') 
	plt.xlabel('$Temperature [\mu K]$', fontsize='x-large')
	plt.ylabel('Number of pixels',fontsize='x-large')    
	plt.legend()
	plt.savefig(fileCMB_Maps+'graficos/Hist_'+str(idMap)+'.png')
	
	plt.close()	



# -------------------------- INPUT DATA ------------------------------------: #

# ---- We load the samples saved in pickle files --------: 
# ---- 2MASS galaxy files -------------------------------:
# ---- Sa Sample
path_out_Sa = fileProject+'data/Galaxy_Samples/Sa_sample.pickle'		# Path to the Filename
with open(path_out_Sa, 'rb') as f_Sa:
    Sa_sample = pickle.load(f_Sa)

# ---- Sb Sample
path_out_Sb = fileProject+'data/Galaxy_Samples/Sb_sample.pickle'		# Path to the Filename
with open(path_out_Sb, 'rb') as f_Sb:
    Sb_sample = pickle.load(f_Sb)

# ---- Sc Sample
path_out_Sc = fileProject+'data/Galaxy_Samples/Sc_sample.pickle'		# Path to the Filename
with open(path_out_Sc, 'rb') as f_Sc:
    Sc_sample = pickle.load(f_Sc)

# ---- We set the final sample --------------------------:

if Samp_choice == 0.04:		# --- 0.003 < z < 0.04
	Sa_0003_004 = Sa_sample[ ( Sa_sample['v'] < z2v(0.04) ) & ( Sa_sample['v'] > z2v(0.003) ) ]
	Sb_0003_004 = Sb_sample[ ( Sb_sample['v'] < z2v(0.04) ) & ( Sb_sample['v'] > z2v(0.003) ) ]
	Sc_0003_004 = Sc_sample[ ( Sc_sample['v'] < z2v(0.04) ) & ( Sc_sample['v'] > z2v(0.003) ) ]
	
	if Rext_cut == True:
		Sa_0003_004 = Sa_0003_004[ ( 10**(Sa_0003_004['r_ext']) < 120  ) ]
		Sb_0003_004 = Sb_0003_004[ ( 10**(Sb_0003_004['r_ext']) < 120  ) ]
		Sc_0003_004 = Sc_0003_004[ ( 10**(Sc_0003_004['r_ext']) < 120  ) ]
#	else:
#		continue
		
	Sa_sample = Sa_0003_004
	Sb_sample = Sb_0003_004
	Sc_sample = Sc_0003_004

elif Samp_choice == 0.03:	# --- 0.003 < z < 0.03
	Sa_0003_003 = Sa_sample[ ( Sa_sample['v'] < z2v(0.03) ) & ( Sa_sample['v'] > z2v(0.003) ) ]
	Sb_0003_003 = Sb_sample[ ( Sb_sample['v'] < z2v(0.03) ) & ( Sb_sample['v'] > z2v(0.003) ) ]
	Sc_0003_003 = Sc_sample[ ( Sc_sample['v'] < z2v(0.03) ) & ( Sc_sample['v'] > z2v(0.003) ) ]

	if Rext_cut == True:
		Sa_0003_003 = Sa_0003_003[ ( 10**(Sa_0003_003['r_ext']) < 120  ) ]
		Sb_0003_003 = Sb_0003_003[ ( 10**(Sb_0003_003['r_ext']) < 120  ) ]
		Sc_0003_003 = Sc_0003_003[ ( 10**(Sc_0003_003['r_ext']) < 120  ) ]
#	else:
#		continue

	Sa_sample = Sa_0003_003
	Sb_sample = Sb_0003_003
	Sc_sample = Sc_0003_003
	
elif Samp_choice == 0.02:	# --- 0.003 < z < 0.02
	Sa_0003_002 = Sa_sample[ ( Sa_sample['v'] < z2v(0.02) ) & ( Sa_sample['v'] > z2v(0.003) ) ]
	Sb_0003_002 = Sb_sample[ ( Sb_sample['v'] < z2v(0.02) ) & ( Sb_sample['v'] > z2v(0.003) ) ]
	Sc_0003_002 = Sc_sample[ ( Sc_sample['v'] < z2v(0.02) ) & ( Sc_sample['v'] > z2v(0.003) ) ]

	if Rext_cut == True:
		Sa_0003_002 = Sa_0003_002[ ( 10**(Sa_0003_002['r_ext']) < 120  ) ]
		Sb_0003_002 = Sb_0003_002[ ( 10**(Sb_0003_002['r_ext']) < 120  ) ]
		Sc_0003_002 = Sc_0003_002[ ( 10**(Sc_0003_002['r_ext']) < 120  ) ]
#	else:
#		continue

	Sa_sample = Sa_0003_002
	Sb_sample = Sb_0003_002
	Sc_sample = Sc_0003_002

elif Samp_choice == 0.015:	# --- 0.003 < z < 0.015
	Sa_0003_0015 = Sa_sample[ ( Sa_sample['v'] < z2v(0.015) ) & ( Sa_sample['v'] > z2v(0.003) ) ]
	Sb_0003_0015 = Sb_sample[ ( Sb_sample['v'] < z2v(0.015) ) & ( Sb_sample['v'] > z2v(0.003) ) ]
	Sc_0003_0015 = Sc_sample[ ( Sc_sample['v'] < z2v(0.015) ) & ( Sc_sample['v'] > z2v(0.003) ) ]

	if Rext_cut == True:
		Sa_0003_0015 = Sa_0003_0015[ ( 10**(Sa_0003_0015['r_ext']) < 120  ) ]
		Sb_0003_0015 = Sb_0003_0015[ ( 10**(Sb_0003_0015['r_ext']) < 120  ) ]
		Sc_0003_0015 = Sc_0003_0015[ ( 10**(Sc_0003_0015['r_ext']) < 120  ) ]
#	else:
#		continue

	Sa_sample = Sa_0003_0015
	Sb_sample = Sb_0003_0015
	Sc_sample = Sc_0003_0015
	
elif Samp_choice == 0.01:	# --- 0.003 < z < 0.01
	Sa_0003_001 = Sa_sample[ ( Sa_sample['v'] < z2v(0.01) ) & ( Sa_sample['v'] > z2v(0.003) ) ]
	Sb_0003_001 = Sb_sample[ ( Sb_sample['v'] < z2v(0.01) ) & ( Sb_sample['v'] > z2v(0.003) ) ]
	Sc_0003_001 = Sc_sample[ ( Sc_sample['v'] < z2v(0.01) ) & ( Sc_sample['v'] > z2v(0.003) ) ]

	if Rext_cut == True:
		Sa_0003_001 = Sa_0003_001[ ( 10**(Sa_0003_001['r_ext']) < 120  ) ]
		Sb_0003_001 = Sb_0003_001[ ( 10**(Sb_0003_001['r_ext']) < 120  ) ]
		Sc_0003_001 = Sc_0003_001[ ( 10**(Sc_0003_001['r_ext']) < 120  ) ]
#	else:
#		continue

	Sa_sample = Sa_0003_001
	Sb_sample = Sb_0003_001
	Sc_sample = Sc_0003_001

	

# -------------------------- BLOCK OF PROCESS  ------------------------------------:
# ---- We get the position vectors (norm=1) of the gxs --:
Sa_vec = hp.ang2vec( np.pi/2. - Sa_sample['b']*np.pi/180., Sa_sample['l']*np.pi/180., lonlat=False )
Sb_vec = hp.ang2vec( np.pi/2. - Sb_sample['b']*np.pi/180., Sb_sample['l']*np.pi/180., lonlat=False ) 
Sc_vec = hp.ang2vec( np.pi/2. - Sc_sample['b']*np.pi/180., Sc_sample['l']*np.pi/180., lonlat=False ) 
	# ang2vec returns a (numpy) array
RextVec_Sa = [ Sa_sample['r_ext'].to_numpy(), Sa_vec ]
RextVec_Sb = [ Sb_sample['r_ext'].to_numpy(), Sb_vec ]
RextVec_Sc = [ Sc_sample['r_ext'].to_numpy(), Sc_vec ]
	# Let us note that we convert the 1st element to a (numpy) array
	       
S_Sample = [Sa_sample, Sb_sample, Sc_sample]
RextVec_pairs = [RextVec_Sa, RextVec_Sb, RextVec_Sc]

print('###############################################################')
print( 'The checks below should be True')
print('###############################################################')
print('')      
print( 'Sa_sample has the same len than Sa_vec?', len(RextVec_Sa[0]) == len(RextVec_Sa[1]) )
print( 'Sb_sample has the same len than Sb_vec?', len(RextVec_Sb[0]) == len(RextVec_Sb[1]) )
print( 'Sc_sample has the same len than Sc_vec?', len(RextVec_Sc[0]) == len(RextVec_Sc[1]) )
print('')   
print('###############################################################')

# ---- Main loops ---------------------------------------:
t_cpu_start = time.clock()
t_running_start = time.time()
                
for j in range(1,5):
	T_data, T_mask = Rd_map(CMBmap[j-1])
	print('###############################################################')
	print('CMB map in current use', idMap[j-1])
	print('###############################################################')	

    
	for kk in range(len(S_Sample)):
		print('')
		print('###############################################################')
		print('Morphological type of the galaxy sample', id_Sample[kk] )	
		print('###############################################################')

		RextVec = RextVec_pairs[kk]
		Rext_gxs = 10**(RextVec[0]) 		# The array contains r_ext of the gxs in the sample.
		Vec_gxs = RextVec[1]			# The array contains pos+vector for each gxy in the sample.		
		Pix_disks = set()
		N_Disk_in_Mask = 0
		N_in_mask_S = []
	
		for k in range( len(Rext_gxs) ):
			R_disk_ext = (Rext_gxs[k])*arcsec_to_rad		# Radius of galaxy k in radians
#			R_disk_ext = (50.0)*arcsec_to_rad			# Radius of galaxy k in radians
			Disk_k = hp.query_disc( nside=Nside, vec= Vec_gxs[k],\
								   radius= R_disk_ext*x_Rext, inclusive=True, nest=False )
					# We realize that it is better to put inclusive=True.
				
			# We check that all the pixels belong to the CMB maps and no one of them drop into the mask.
			if ( np.all( np.ones(len(Disk_k)) == T_mask[Disk_k] ) ):		
				Pix_disks.update(Disk_k)				# Set collecting the values of T.
			# If not we count the number of galaxies inside the mask.
			else:
				N_Disk_in_Mask = N_Disk_in_Mask + 1			
												
		Mask_S = list(Pix_disks)
		T_mask[Mask_S] = np.zeros(len(Mask_S))
		N_in_mask_S.append(N_Disk_in_Mask)					# We add the number of disks in the mask	
		T_Mask_S = np.where(T_mask, T_data, hp.UNSEEN)
		CMB_MwPlot(T_Mask_S, Tmin, Tmax, idMap[j-1],id_Sample[kk])           
		T_mask[Mask_S] = np.ones(len(Mask_S))
		
t_cpu_end = time.clock()
t_running_end = time.time()


print('###############################################################')
print('CPU elapsed time', t_cpu_end - t_cpu_start )
print('Total running time', t_running_end - t_running_start )		
print('###############################################################')
