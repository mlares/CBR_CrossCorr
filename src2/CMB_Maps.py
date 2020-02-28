	# --------------------- IMPORT MODULES ------------------------------------ : #

	import numpy as np
	import matplotlib.pyplot as plt
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
	#colormap = 'viridis' 
	#colormap = 'plasma' 
	colormap = 'inferno' 
	#colormap = 'magma'
	#colormap = 'cividis'

	Tmax = 0.0005*K_to_muK
	Tmin = -0.0005*K_to_muK

	print('###############################################################')
	print('# USEFUL DEFINITIONS AND PARAMETERS                           #')
	print('###############################################################')
	print('')
	print('Nisde resolution of the CMB maps = ', Nside)
	print('Color pallete CMB maps = ', colormap)
	print('')
	print('###############################################################')
	print('')
	print('')


	# ---- CMB Maps data files ------------------------------: 

	#fileCMB = '/home/ezequiel/Documentos/Planck+Data+Product+Release_2019/'
	fileProject = '/home/zboero/Projects/CMB/Maps/'
	fileCMB = '/home/zboero/Documentos/Planck+Data+Product+Release_2019/'
	fileMask = '/home/zboero/Documentos/Planck+Data+Product+Release_2019/Masks/'

	SMICA=(fileCMB+'CMB_Maps/COM_CMB_IQU-smica_2048_R3.00_full.fits')
	SEVEM=(fileCMB+'CMB_Maps/COM_CMB_IQU-sevem_2048_R3.01_full.fits')
	NILC=(fileCMB+'CMB_Maps/COM_CMB_IQU-nilc_2048_R3.00_full.fits')
	Commander=(fileCMB+'CMB_Maps/COM_CMB_IQU-commander_2048_R3.00_full.fits')
	common_mask = (fileMask+'COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits')

	CMBmap =[SMICA, SEVEM, NILC, Commander]
	idMap = ['SMICA', 'SEVEM', 'NILC', 'Commander']	# ID's for the maps



	# -------------------------- FUNCTIONS  -------------------------------------: #

	def Rd_map(filename):			# It reads the CMB maps.
		CMB_T, header = hp.read_map( filename, nest=False, h=True, field=(0,3) )
		T_data = CMB_T[0]           	# Temperature (numpy) array 
		T_mask = CMB_T[1]               	# Mask (numpy) array with 0's and 1's.
		return( T_data, T_mask )

	def Rd_map0(filename):			# It reads the CMB maps.
		T_data, header = hp.read_map( filename, nest=False, h=True, field=(0) )
		return(T_data) 


	def CMB_MwPlot(T_data, Tmin, Tmax, idMap):
		plt.close()
		plt.figure().patch.set_facecolor('xkcd:white')    
		hp.mollview( T_data*K_to_muK, coord='G',unit='Temperature [$\mu$K]', xsize=800,\
					title='CMB Temperature Map ('+str(idMap)+')',\
					cbar=True, cmap=colormap, min=Tmin, max=Tmax,\
					badcolor='w' )	
		plt.show()
		plt.savefig(fileProject+'graficos/'+str(idMap)+'_Tmask.png')
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
		plt.savefig(fileProject+'graficos/Hist_'+str(idMap)+'.png')

		plt.close()	


	# -------------------------- BLOCK OF PROCESS  ------------------------------------: 
	t_cpu_start = time.clock()
	t_running_start = time.time()

	print('')
	print('Loading the common mask')
	print('')
	#T_mask, head_Mask = hp.read_map(common_mask, nest=False, h=True)

	print('###############################################################')
	print('Common mask loaded...' )
	print('###############################################################')

	for j in range(1,5):

		print('')
		print('Loading map...')
		print('')

		T_data, T_mask = Rd_map(CMBmap[j-1])

		print('###############################################################')
		print('Current map in use: ', str(idMap[j-1]) )
		print('###############################################################')

		T_conf = np.where(T_mask, T_data, hp.UNSEEN)
		CMB_MwPlot(T_conf, Tmin, Tmax, idMap[j-1])           
		CMB_Hist(T_data, T_mask, idMap[j-1])


	t_cpu_end = time.clock()
	t_running_end = time.time()


	print('###############################################################')
	print('CPU elapsed time', t_cpu_end - t_cpu_start )
	print('Total running time', t_running_end - t_running_start )		
	print('###############################################################')
                
    
