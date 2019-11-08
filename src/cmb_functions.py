def readmap(filename):
    '''
    Reads the CMB map

    Args:
        filename (str): the file name of the map to be read

    Raises:

    Returns:
        readmap: a healpix map, class ?

    '''

    import healpy as hp
    
    hp_data, hp_header = hp.read_map(filename, h=True,field=(0,3))

    hp_data_sel = hp_data[0] # <---- elijo la columna de los datos
    hp_mask = hp_data[1]     # <---- elijo la columna de la máscara

    return(hp_data_sel,hp_mask)


def profdata(nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask):
    """
    Args:
        nside (int): the healpix nside parameter
        fac (number):
        Dsel:
        vec:
        hp_data_sel:
        hp_mask:
    Raises:
        errors?
    Returns:
        profdata:
        Kr ():
        serr ()
    """

    import numpy as np
    import healpy as hp
    from scipy.stats import sem

    Kr=np.zeros(len(rprof))
    serr=np.zeros(len(rprof))

    for k in range(len(Kr)):

        Kd=[]

        for i in range(len(vec)):

            listpix = hp.query_disc(nside, vec[i], ((10.**(Dsel['r_ext'][i]))/206264.99992965)*rprof[k],inclusive=True, fact=4, nest=False) #radianes

            listpix_prev=hp.query_disc(nside, vec[i], ((10.**(Dsel['r_ext'][i]))/206264.99992965)*rprof[k-1],inclusive=True, fact=4, nest=False) 

            if(k==0):

                listpix_prev=[]

            listpix_mask=[]

            for j in range(len(listpix)):

                if(hp_mask[listpix[j]]==1):
    
                    if(listpix[j] not in listpix_prev):             

                        listpix_mask.append(listpix[j])

            meanK_disk=np.nanmean(hp_data_sel[listpix_mask])

            Kd.append(meanK_disk)

        serr[k]=sem(Kd,nan_policy='omit')*10.**fac

        Kr[k]=np.nanmean(Kd)*(10.**fac)

    return(Kr,serr)


def profran(nran,nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask):
    """
    Args:
        nran ():
        nside (int): the healpix nside parameter
        fac (number):
        rprof ():
        Dsel:
        vec:
        hp_data_sel:
        hp_mask:
    Raises:
        errors?
    Returns:
        profran:
    """

    import numpy as np
    import healpy as hp
    from scipy.stats import sem


    Kr_ran=np.zeros((nran,len(rprof)))
    for r in range(nran):
        
        np.random.seed()

        vec_ran=vec

        for rr in range(len(vec)):

            phiran = np.random.uniform(0,2*np.pi,1)
            costheta = np.random.uniform(-1,1,1)
            u = np.random.uniform(0,1)

            thetaran = np.arccos(costheta)
            rran = u**(1/3)

            xran = rran * np.sin(thetaran) * np.cos(phiran)
            yran = rran * np.sin(thetaran) * np.sin(phiran)
            zran = rran * np.cos(thetaran)

            vec_ran[rr,0]=xran
            vec_ran[rr,1]=yran
            vec_ran[rr,2]=zran



        for k in range(len(Kr_ran[r,:])):

            Kd_ran=[]

            for i in range(len(vec_ran)):

                listpix_ran = hp.query_disc(nside, vec_ran[i], ((10.**(Dsel['r_ext'][i]))/206264.99992965)*rprof[k],inclusive=True, fact=4, nest=False) #radianes
                listpix_ran_prev=hp.query_disc(nside, vec_ran[i], ((10.**(Dsel['r_ext'][i]))/206264.99992965)*rprof[k-1],inclusive=True, fact=4, nest=False) 

                if(k==0):
                    listpix_ran_prev=[]

                listpix_mask_ran=[]

                for j in range(len(listpix_ran)):

                    if(hp_mask[listpix_ran[j]]==1):
                    
                        if(listpix_ran[j] not in listpix_ran_prev):             
    
                            listpix_mask_ran.append(listpix_ran[j])

                meanK_disk_ran=np.nanmean(hp_data_sel[listpix_mask_ran])

                Kd_ran.append(meanK_disk_ran)

       
            Kr_ran[r,k]=np.nanmean(Kd_ran)*(10.**fac)

    
    return(Kr_ran)


