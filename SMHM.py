import numpy as np
from scipy import stats
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import illustris_python as il
import h5py

mpl.rcParams['agg.path.chunksize'] = 10000

parttype_stars = 4

if __name__=='__main__':
    
    # information about simulation
    ubications = ['Illustris','IllustrisTNG']
    simulations = ['L75n1820FP','L75n1820TNG']
    simulation_name = ['ILLUSTRIS-1','TNG100-1' ]
    snapnums = [135,99]  # for z = 0
    colors = ['r','b']

    mstarmin = 7.5e8  #Msun
    mstarmax = 3e12 #Msun
    mhalomin = 9e9  #Msun
    mhalomax = 3.5e14 #Msun
    # Asignar mass stellar in bins
    nbins = 20 
#    mstarmin = 2e8  #Msun
#    mstarmax = 2e12 #Msun
    mhalo_bin_edges = np.logspace(np.log10(mhalomin), np.log10(mhalomax), nbins+1)
    mhalo_bin_centeres = np.sqrt(mhalo_bin_edges[:-1] * mhalo_bin_edges[1:]) # geometric mean
    #Nbins = len(mhalo_bin_centeres)
    mstellar_bin_edges = np.logspace(np.log10(mstarmin), np.log10(mstarmax), nbins+1)
    mstellar_bin_centeres = np.sqrt(mstellar_bin_edges[:-1] * mstellar_bin_edges[1:])
    

    # create figure
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)

    #fig = plt.figure(figsize=(8,8))
    #ax2 = fig.add_subplot(111)
    #ax1 = fig.add_subplot(221)

    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5), tight_layout=True)
    #ax1,ax2 = fig.add_subplot(111)
    #fig, axes = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
    
    #fig  = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
    #ax1,ax2 = axes.flatten()
    #fig = plt.figure(figsize=(6,6))
    #fig.add_subplot(111)
    #ax2 = fig.add_subplot(222)
       
    #ax1,ax2 = fig.add_subplot(111)

    for k in range(len(simulations)):
        ubication = ubications[k]
        simulation = simulations[k]
        snapnum = snapnums[k]
        
#        """        
        # selecc the folder with data
        if ubication == 'Illustris':
            #basedir = '/n/hernquistfs1/%s/Runs/%s/output' % (ubication, simulation)
            basedir = '/n/holystore01/LABS/hernquist_lab/Lab/hernquistfs1/%s/Runs/%s/output' % (ubication, simulation)
        elif ubication == 'IllustrisTNG':
            #basedir = '/n/hernquistfs3/%s/Runs/%s/output' % (ubication, simulation)
            basedir = '/n/hernquistfs3/%s/Runs/%s/output' % (ubication, simulation)
        else:
            raise NotImplementedError(ubication)

         
            # select parameters of simulation 
        snap_filename = '%s/snapdir_%03d/snap_%03d.0.hdf5' % (basedir, snapnum, snapnum)
        with h5py.File(snap_filename, 'r') as f_snap:
            header = f_snap['Header']
            h = header.attrs['HubbleParam'] 
            z = header.attrs['Redshift'] 
            box_size = header.attrs['BoxSize'] #Mpc
            # ..
            box_size_Mpc = box_size / 1000.0 /(1.0 + z) / h
            volume_Mpc = box_size_Mpc**3
   
        #####----MASS HALO-----------------#########    
        # load stellar Mass (SubhaloMassInHalfRadType)
        #mstar_in_rad = il.groupcat.loadSubhalos(
            #basedir, snapnum, fields=['SubhaloMassInHalfRadType'])[:, parttype_stars]
        mhalo_200 = il.groupcat.loadHalos(
            basedir, snapnum, fields=['Group_M_Crit200'])
        
        # covert  units to Msun
        mhalo_200_Msun = mhalo_200 * 1e10 /h
        print('num halos',len(mhalo_200_Msun))
        #####----------MASS STELLAR- SUBFIND-------------#########
        
        # load stellar Mass (SubhaloMassType)
        #mstar = il.groupcat.loadSubhalos(
         #   basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
        # load stellar Mass ((SubhaloMassInRadType)
        mstar = il.groupcat.loadSubhalos(
            basedir, snapnum, fields=['SubhaloMassInRadType'])[:, parttype_stars]

        # Valid location
        locs_mstar =  mstar > 0
        
        # covert  units to Msun
        mstar_Msun = mstar * 1e10/h
        print('num star',len(mstar_Msun))
        
        ##===Selecionar el first/primary/most massive Subfind group===##
        
        group_first_sub = il.groupcat.loadHalos(
            basedir, snapnum, fields=['GroupFirstSub'])
        

        
        ####--------MASS OF THE GALXIE PARENT (PRINCIPAL)--------####
        """
        Condicion es para solo tomar la galacia mas importante o padre
        para esto se utiliza la funcion 'SubhaloGrNr', la cual relaciona 
        el index de cada halo con subfind. 
        """
        sub_gr_nr = il.groupcat.loadSubhalos(
            basedir, snapnum, fields=['SubhaloGrNr'])
        #-->Index into the Group table of the FOF host/parent of this Subhalo.
        is_central = group_first_sub[sub_gr_nr] == np.arange(len(sub_gr_nr), dtype=np.int32)
        
        mstar_central_galx_Msun = mstar_Msun[is_central]
        mhalo_200_central_Msun = mhalo_200_Msun[sub_gr_nr[is_central]]
        
        # Valid location
        mhalo200 = mhalo_200[sub_gr_nr]
        locs_mhalo_200 = (mhalo200 > 0) & is_central
        mhalo200_Msun = mhalo200 *1e10 / h
        
        print('Halo Mass max: %e'%max(mhalo_200_central_Msun))
        print('Stellar Mass max: %e'%max(mstar_central_galx_Msun))
        
        #==*********************************************==#
        #PARA QUITARME LAS MASAS CENTRALES IGUALES A CERO
        #==*********************************************==#
        """
        mstar_central_galx_Msun = mstar_central_galx_Msun[mhalo_200_central_Msun > mhalomin]
        mhalo_200_central_Msun = mhalo_200_central_Msun[mhalo_200_central_Msun > mhalomin]
        print('%e'%min(mhalo_200_central_Msun))
        """
        
        """
        #####----------MASS STELLAR- FoF-------------#########
        
        mhalo_stellar = il.groupcat.loadHalos(
            basedir, snapnum, fields=['GroupMassType'])[:,parttype_stars]
        # covert  units to Msun
        mhalo_stellar_Msun = mhalo_stellar * 1e10 /h

        ####--------definir rango de masas --------##
        mhalo_stellar_Msun = mhalo_stellar_Msun[mhalo_200_Msun>mhalomin]
        mhalo_200_Msun = mhalo_200_Msun[mhalo_200_Msun>mhalomin]
        """
               

        #===========  PLOTS  =============

        """                        
        ax1.plot(mhalo_200_Msun,mhalo_stellar_Msun/mhalo_200_Msun, label= simulation, c = colors[k], lw=2.3 )
        ax1.set_xlabel(r'$Halo \ Mass\, [M_{200c}] [M_{Sun}]$', fontsize = 14)
        ax1.set_ylabel(r'$M_{stars}/M_{halo}$', fontsize = 14)
        # Finish plot
        ax1.legend(loc=1, fontsize=14, frameon=False)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax1.set_xlim(mstarmin, mstarmax)
        """        
        #frac_baryon = 0.155 # universal baryon fraction

        #SMHM = mstar_central_galx_Msun/ mhalo_200_central_Msun

        #Phi =  mstar_central_galx_Msun/ mhalo_200_central_Msun/ frac_baryon        
        
        hist_mhalo200, ignored = np.histogram(
            mhalo_200_central_Msun , bins = mhalo_bin_edges) 
        #print(hist_mhalo200)
        #hist_mhalo200, ignored = np.histogram(
         #   mhalo_200_central_Msun, bins = mhalo_bin_edges)
        #hist_mstellar, ignored = np.histogram(
        #    mstar_central_galx_Msun , bins = mhalo_bin_edges)
        #hist_mhalo200, ignored = np.histogram(
        #    Phi, bins = 15)
    
        #x = np.linspace(mhalomin, mhalomax, len(mhalo_200_central_Msun))
        #x = mhalo_200_central_Msun
           
    
        x = mhalo_200_central_Msun
        #x = mhalo200_Msun
        #print(len(x),len(mhalo_200_central_Msun))
        Mstar_halo200_central_median, mhalo_bin_edges_static, binnumber  = stats.binned_statistic(x, mstar_central_galx_Msun, statistic='median', bins = mhalo_bin_edges)
        #Mstar_halo200_central_median, mhalo_bin_edges_static, binnumber  = stats.binned_statistic(x, mstar_Msun, statistic='median', bins = mhalo_bin_edges)
        """
        #y = mstar_Msun
        Mstar_central_static, mstar_bin_edges_static, binnumber  = stats.binned_statistic(x, mstar_central_galx_Msun , statistic='median', bins = mstellar_bin_edges)
        """
        ##=====Funciones para la estadistica======##
        
        #=== Percentil ====
        
        def percentil_bins( p, bin_edges, bin_centeres,locs_mstar, locs_mhalo_200, mh200, mstar,  Nbins):
            """
            a = array for compute percential
            percen = numero percentil to compute
            Nbins = values numbres of bines
            """
            
            percentil = np.zeros(Nbins)
            median = np.zeros(Nbins)
            locs_common = locs_mhalo_200 & locs_mstar
            mhalo_common = mh200[locs_common]
            mstar_common = mstar[locs_common]
            locs_valid = np.zeros(Nbins, dtype=np.bool8)
            #percentil = []
            #loc = np.zeros(Nbins)
            #ones = np.ones(len(a))
            #print(Nbins)
            #dx = (xf-xi)/Nbins
            #print('dx=  %e'%dx)
            
            
            for i in range(Nbins):
                
                #xmin = (xi + i*dx)
                #xmax = (xi + (i+1)*dx)
                xmin = bin_edges[i]
                xmax = bin_edges[i+1]  
                loc = (mhalo_common >= xmin)&(mhalo_common < xmax)
                #aux = mstar[loc]
                if np.sum(loc) > 2:
                    percentil[i] =  np.percentile(mstar_common[loc], p)
                    #median[i] = np.percentile(mstar_common[loc], 50)
                    median[i] = np.median(mstar_common[loc])
                    locs_valid[i] = True
                #print('%e'%xmin)
                #print('%e'%xmax)
                #print(loc)
                #print(aux)
                #print(len(aux))
                """
                if len(aux) == 0:
                    pass
                else:
                    percentil[i] =  np.percentile(aux, p)  
                    #percentil.append( np.percentile(aux, percen) )
                #print(percentil)
                """
                 
                #percentil.append( np.nanpercentile(aux, p) )
                #print('percentile')
                #print(percentil)
            percentil = percentil[locs_valid]
            median = median[locs_valid]
            bin_centeres = bin_centeres[locs_valid]
            return percentil, median, bin_centeres

            

        """
        def median_bins(a, xi, xf, Nbins):
            
            #a = array for compute percential
            #xi = limit initial data
            #xf = limit final data
            #Nbins = values numbres of bines
            
            median = np.zeros(Nbins)
            #loc = np.zeros(Nbins)
            
            for i in range(Nbins):
                dx = (xf-xi)/Nbins
                xmin = (xi + i+dx)
                xmax = (xi + (i+1)*dx)
                loc = (a >= xmin)&(a < xmax)
                aux = a[loc]
                
                median[i] =  np.median(aux)  
                
            return median
        """    
        bin_edges = np.logspace(np.log10(1),np.log10(100),num = nbins+1)
        bin_centeres = np.sqrt(bin_edges[:-1] * bin_edges[1:])
        
        #======== HALO  =============
        Mhalo_percentil_16, median, bin_centeres_m200 = percentil_bins( 16., mhalo_bin_edges, mhalo_bin_centeres,locs_mstar, locs_mhalo_200, mhalo200_Msun, mstar_Msun , nbins)
        
        Mhalo_percentil_84, median, bin_centeres_m200 = percentil_bins( 84., mhalo_bin_edges, mhalo_bin_centeres, locs_mstar, locs_mhalo_200, mhalo200_Msun, mstar_Msun , nbins)
        
        #Mhalo_percentil_84 = percentil_bins,(mhalo_200_central_Msun,84, mhalomin, mhalomax, nbins)
        #Mhalo_median = median_bins(mhalo_200_central_Msun, mhalomin, mhalomax, nbins)
        
        ##======== STELLAR ===========
        
        #Mstellar_percentil_16 = percentil_bins( 16, mstellar_bin_edges, locs_mstar, locs_mhalo_200, mhalo_200_central_Msun, mstar_central_galx_Msun , nbins)
        #np.percentile(mstar_central_galx_Msun,16)
        #Mstellar_percentil_84 = percentil_bins( 84, mstellar_bin_edges, locs_mstar, locs_mhalo_200, mhalo_200_central_Msun, mstar_central_galx_Msun , nbins)
        
        #stellar_percentil_84 = percentil_bins(mstar_central_galx_Msun, 84, mstarmin, mstarmax, nbins)
        #np.percentile(mstar_central_galx_Msun,16)

        #Mstellar_median = median_bins(mstar_central_galx_Msun, mstarmin, mstarmax, nbins)
        
        

        #ax1.plot(mhalo_bin_edges_static[:-1],  Mstar_halo200_central_median, label = simulation_name[k], c=colors[k],lw=1.5 )
        ax1.plot(bin_centeres_m200,  Mstar_halo200_central_median, label = simulation_name[k], c=colors[k],lw=1.5 )

        
        ax1.plot(bin_centeres_m200,  median, '.' , label = 'media', c=colors[k],lw=2 )

        
        #ax1.fill_between(mhalo_bin_edges_static[:-1], Mstellar_percentil_16, Mstellar_percentil_84,  alpha = 0.2, label = 'percentil')

        
        alpha = [0.3, 0.2]
        #ax1.fill_between(mhalo_bin_centeres, Mhalo_percentil_16, Mhalo_percentil_84,  alpha = alpha[k], label = 'percentil-%s'%simulation_name[k])
        ax1.fill_between(bin_centeres_m200, Mhalo_percentil_16, Mhalo_percentil_84,  alpha = alpha[k], label = 'percentil-%s'%simulation_name[k])




        ax1.set_xlabel(r'$M_{200,crit}\, [M_{Sun}]$', fontsize = 14)
        ax1.set_ylabel(r'$M_{stellar}\, [M_{Sun}]$', fontsize = 14)
        ax1.set_title(r'bins=%s'%nbins, fontsize = 15)
        # Finish plot
        ax1.legend(loc=0, fontsize=14, frameon=False)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(mhalomin, mhalomax)

        """
        ax2.plot(mhalo_200_central_Msun,mstar_central_galx_Msun,'.', alpha=0.1)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel(r'$M_{200,crit}\, [M_{Sun}]$', fontsize = 14)
        ax2.set_ylabel(r'$M_{stellar}\, [M_{Sun}]$', fontsize = 14)
        ax2.set_title(r'bins=%s'%nbins, fontsize = 15)
        ax2.set_xlim(mhalomin, mhalomax)
        """

        """
        #===========Mass_in_rad =============
        hist, ignored = np.histogram(mstar_in_rad_Msun,bins =mstar_bin_edges) ## dN = hist
        Phi2 = hist / (d_log_mstar * volume_Mpc)
        ax2.plot(mstar_bin_centeres, Phi2, label = simulation, c=colors[k],lw=1.5 )
        ax2.set_xlabel(r'$Stellar\ half\ mass\ radius\ [M_{Sun}]$', fontsize = 14)
        ax2.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)
        # Finish plot
        ax2.legend(loc=1, fontsize=14, frameon=False)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlim(mstarmin, mstarmax)
#        ax2.text(5e11,5e-3, tex, fontsize=14)
        """

        print('Finished for %s.' %(simulation))



#    ax.set_xlabel(r'$Stellar \ Mass\ [M_{Sun}]$', fontsize = 14)
#    ax.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)

    for ending in ['pdf','png']:
#        fig.suptitle('Mass function galaxies to $z\sim %s$'%(int(z)), fontsize=20)
        fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/SMHM_z%s.%s' %(int(z) , ending))
    

        
