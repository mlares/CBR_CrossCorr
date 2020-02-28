	
# --------------------- IMPORT MODULES ------------------------------------ : #

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
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
Samp_choice = 0.04 	# 0.02 , 0.04, etc.		# This parameter chooses the sample
Rext_cut = True			# If True it cut galaxies with 10**(r_ext) > 120.0 arcsec.

ymax = 5.0
ymin = -15.0
xmax = 50
xmin = 5
nbins_plt = 10
bin_sz = (xmax - xmin)/(nbins_plt - 1)
              
print('###############################################################')
print('# USEFUL DEFINITIONS AND PARAMETERS                           #')
print('###############################################################')
print('')
print('Speed of light in km x s^{-1} =', c_luz)
print('Nisde resolution of the CMB maps = ', Nside)
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

# ---- 2MASS galaxy files -------------------------------: 

#fileProject = '/home/ezequiel/Projects/CMB/Anisotropy_signal_around_ScGalaxies/'	#Local PC path
fileProject = '/home/zboero/Projects/CMB/Anisotropy_signal_around_ScGalaxies/'		#Path in cluster
file_TwoMASS = fileProject+'data/2MASS/2mrs_1175_done.dat'

Columns = ['l', 'b', 'r_iso', 'r_ext', 'b/a', 'type', 'v']
# We use Gal Coord (l,b) since it is appropriated at the moment to work with CMB maps.
hd= 9					# Commented lines of header in the 2MASS catalog file


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
    return( Gxs )

def Plot_map(x, Data, id_Sample , ymax, ymin):
	plt.close()
	plt.figure(0)
	plt.plot( x, np.array(Data[0][0])*K_to_muK , label='SMICA', marker='o')
	plt.plot( x, np.array(Data[1][0])*K_to_muK, label='SEVEM', marker='s')
	plt.plot( x, np.array(Data[2][0])*K_to_muK, label='NILC', marker='p')
	plt.plot( x, np.array(Data[3][0])*K_to_muK, label='Commander', marker='*')
	plt.title('$<{\Delta T}>$ around rings in '+str(id_Sample)+' galaxies', loc='center', fontdict={'fontsize': 14} )
	plt.xlabel('Factor of extrapolation radius $A_k$ [dimensionless]')
	plt.ylabel('$\widebar{\Delta T}$ [$\mu$K]')
	plt.ylim(ymin, ymax)
	plt.axhline(y=0.0, color='grey', linewidth=1.0)
	plt.legend()

	plt.savefig(fileProject+'graficos/RingCMBMeanT_on_'+str(id_Sample)+'_z0003_z'+str(Samp_choice)+'.png')
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
id_Sample = ['Sa_sample', 'Sb_sample', 'Sc_sample']
RextVec_pairs = [RextVec_Sa, RextVec_Sb, RextVec_Sc]

print('###############################################################')
print( 'Number of galaxies in each sample')
print('###############################################################')
print('')      
print( 'Nr of gxs Sa = ', len(Sa_sample) )
print( 'Nr of gxs Sb = ', len(Sb_sample) )
print( 'Nr of gxs Sc = ', len(Sc_sample) )
print('')   
print('###############################################################')

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

for kk in range(len(S_Sample)):
	RextVec = RextVec_pairs[kk]
	Rext_gxs = 10**(RextVec[0]) 		# The array contains r_ext of the gxs in the sample.
	Vec_gxs = RextVec[1]			# The array contains pos+vector for each gxy in the sample.
	
	for j in range( len(idMap) ):
		T_data, T_mask = Rd_map(CMBmap[j])
		
		print('###############################################################')
		print('CMB map in current use', idMap[j] )
		print('Nside resolution', Nside	)	
		print('###############################################################')
		print('')	
		print('###############################################################')
		print('Morphological type of the galaxy sample', id_Sample[kk] )	
		print('###############################################################')
		
		T_array = []			# Initialize the array with the CMB mean T computed around the galaxies.
		N_in_mask_array = []
		sizes = np.linspace(xmin, xmax, nbins_plt)		#sizes = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]
		# (This loop run over all the factors of extrapolation radius)
		for ss in sizes:
			x_Rext = ss
			N_Disk_in_Mask = 0			    # Counts the number of disks in the masked CMB map 
			GalDisk_ipix = []			    # List of pixels for the disks
			Tpix_rings = set()				# Set of values with T CMB from the pixels in disks 
			
			print('')
			print('# ------ Factor of the extrapolation radius A_k = ', x_Rext)
			print('# ------ Size of the bin at this A_k = ', bin_sz)
			print('# ------ Nr of disks in mask when begins the run in A_k'\
			      ' (should be 0) = ', N_Disk_in_Mask)	
			print('# ------ List of pixels when begins the run in A_k'\
			      ' (should be []) = ', GalDisk_ipix)	
			print('# ------ CMB T values for pixels when begins the run in A_k'\
			      ' (should be a set {}) = ', Tpix_rings)
			
			# (This loop run over all the galaxies in the sample)
			for k in range( len(Rext_gxs) ):
				R_disk_ext = (Rext_gxs[k])*arcsec_to_rad		# Radius of galaxy k in radians
#				R_disk_ext = (50.0)*arcsec_to_rad			# Radius of galaxy k in radians
				x_Rext0 = x_Rext - bin_sz
				Disk_k0 = hp.query_disc( nside=Nside, vec= Vec_gxs[k],\
							radius= R_disk_ext*x_Rext0, inclusive=True, nest=False)		
				Disk_k = hp.query_disc( nside=Nside, vec= Vec_gxs[k],\
						       radius= R_disk_ext*x_Rext, inclusive=True, nest=False )
			# We realize that it is better to put inclusive=True.
				Ring_k = np.array( list(set(Disk_k) - set(Disk_k0)) )                
			# We check that all the pixels belong to the CMB maps and no one of them drop into the mask.
				if ( np.all( np.ones(len(Disk_k)) == T_mask[Disk_k] ) ):
					Tpix_ring_k = T_data[Ring_k]				# Numpy array with T values in disk k.
					Tpix_rings.update(Tpix_ring_k)				# Set collecting the values of T.
			# If not we count the number of galaxies inside the mask.
				else:
					N_Disk_in_Mask = N_Disk_in_Mask + 1							
			T_val_ss = np.array( list(Tpix_rings) )
			T_mean_ss = np.nanmean(T_val_ss)
			T_array.append(T_mean_ss)						# We add the T CMB
			N_in_mask_array.append(N_Disk_in_Mask)					# We add the number of disks in the mask
		
		print('')
		print('###############################################################')
		print('')
		print('The number of disk in the mask at the "radius" A_k is: ', N_in_mask_array)
		print('The values of T in CMB at "radius" A_ks is: ', T_array)
		print('')
		print('###############################################################')
		
		if j == 0: # The final data for the CMB maps SMICA							
			SMICA_T_array = T_array
			SMICA_N_in_Mask = N_in_mask_array	
		elif j == 1: # The final data for the CMB maps SEVEM
			SEVEM_T_array = T_array
			SEVEM_N_in_Mask = N_in_mask_array	
		elif j == 2: # The final data for the CMB maps NILC
			NILC_T_array = T_array
			NILC_N_in_Mask = N_in_mask_array
		else: # The final data for the CMB maps Commander
			Commander_T_array = T_array
			Commander_N_in_Mask = N_in_mask_array

# We save the dictionaries as pickle files
	if kk == 0: # ---- Sa Sample
		SMICA_Data = [ SMICA_T_array, SMICA_N_in_Mask]
		SEVEM_Data = [ SEVEM_T_array, SEVEM_N_in_Mask]
		NILC_Data = [ NILC_T_array, NILC_N_in_Mask]
		Commander_Data = [ Commander_T_array, Commander_N_in_Mask]
		
		Data_Sa = [ SMICA_Data, SEVEM_Data, NILC_Data, Commander_Data ]
		Plot_map( sizes, Data_Sa, id_Sample[0], ymax, ymin)
		
	# The final data for the CMB maps SEVEM
	elif kk == 1: # ---- Sb Sample
		SMICA_Data = [ SMICA_T_array, SMICA_N_in_Mask]
		SEVEM_Data = [ SEVEM_T_array, SEVEM_N_in_Mask]
		NILC_Data = [ NILC_T_array, NILC_N_in_Mask]
		Commander_Data = [ Commander_T_array, Commander_N_in_Mask]
		
		Data_Sb = [ SMICA_Data, SEVEM_Data, NILC_Data, Commander_Data ]
		Plot_map( sizes, Data_Sb, id_Sample[1], ymax, ymin)
		
	# The final data for the CMB maps Commander
	else:			# ---- Sc Sample
		SMICA_Data = [ SMICA_T_array, SMICA_N_in_Mask]
		SEVEM_Data = [ SEVEM_T_array, SEVEM_N_in_Mask]
		NILC_Data = [ NILC_T_array, NILC_N_in_Mask]
		Commander_Data = [ Commander_T_array, Commander_N_in_Mask]
		
		Data_Sc = [ SMICA_Data, SEVEM_Data, NILC_Data, Commander_Data ]
		Plot_map( sizes, Data_Sc, id_Sample[2], ymax, ymin)
		
	
t_cpu_end = time.clock()
t_running_end = time.time()
	

print('###############################################################')
print('CPU elapsed time', t_cpu_end - t_cpu_start )
print('Total running time', t_running_end - t_running_start )		
print('###############################################################')


