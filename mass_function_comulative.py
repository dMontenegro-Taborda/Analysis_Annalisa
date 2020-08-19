import numpy as np
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
    #only TNG
    ubications = ['IllustrisTNG','IllustrisTNG','IllustrisTNG','IllustrisTNG']
    simulations = ['L205n2500TNG','L75n1820TNG','L205n2500TNG','L75n1820TNG']
    #simulation = ['TNG300-1','TNG100-1' ]
    snapnums = [99,99,50,50]  # for z = 0,1
    colors = ['r','orange', 'b','c']

    name_simulations = ['TNG300-1 z=0','TNG100-1 z=0','TNG300-1 z=1','TNG100-1 z=1']

    """
    ubications = ['Illustris','IllustrisTNG']
    simulations = ['L75n1820FP','L75n1820TNG']
    #simulation = ['ILLUSTRIS-1','TNG100-1' ]
    
    snapnums = [135,99]  # for z = 0
    colors = ['r','b']
    """
    mstarmin = 2e8  #Msun
    mstarmax = 2e12 #Msun
    mhalomin = 1e13  #Msun
    mhalomax = 3e15 #Msun

    
    # Asignar mass stellar in bins
    nbins = 25
#    mstarmin = 2e8  #Msun
#    mstarmax = 2e12 #Msun
    mhalo_bin_edges = np.logspace(np.log10(mhalomin), np.log10(mhalomax), nbins+1)
    mhalo_bin_centeres = np.sqrt(mhalo_bin_edges[:-1] * mhalo_bin_edges[1:]) # geometric mean
#    d_mstar = mhalo_bin_edges[1:]-mhalo_bin_edges[:-1] # bin Delta Mass
#    d_log_m = np.log10(mstar_bin_edges[1:])-np.log10(mstar_bin_edges[:-1]) # bin Delta log Mass
    

    # create figure
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(111)
#    ax1 = fig.add_subplot(221)
    #    ax1,ax2 = fig.add_subplot(111)
    #    fig, axes = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
    
#    fig  = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
#    ax1,ax2 = axes.flatten()
#    fig = plt.figure(figsize=(6,6))
#    fig.add_subplot(111)
#    ax2 = fig.add_subplot(222)
       
#    ax1,ax2 = fig.add_subplot(111)

    
    for k in range(len(simulations)):
        ubication = ubications[k]
        simulation = simulations[k]
        snapnum = snapnums[k]
        name_simulation = name_simulations[k]
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
        mhalo_200_Msun = mhalo_200 * 1e10/h

        #####----------MASS STELLAR- SUBFIND-------------#########
        
        # load stellar Mass (SubhaloMassType)
        mstar = il.groupcat.loadSubhalos(
            basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
        # covert  units to Msun

        mstar_Msun = mstar * 1e10/h

        
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
        is_central = group_first_sub[sub_gr_nr] == np.arange(len(mstar), dtype=np.int32)
        
        mstar_central_galx_Msun = mstar_Msun[is_central]
        mhalo_200_central_Msun = mhalo_200_Msun[sub_gr_nr[is_central]]

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
        
        
        
        #hist_mhalo200, ignored = np.histogram(
         #   mhalo_200_Msun , bins = mhalo_bin_edges) 
        hist_mhalo200, ignored = np.histogram(
            mhalo_200_central_Msun , bins = mhalo_bin_edges)
        hist_mstellar, ignored = np.histogram(
            mstar_central_galx_Msun , bins = mhalo_bin_edges)
        #histogram comulative
        hist_comulative_Mhalo200 = np.sum(hist_mhalo200)-np.cumsum(hist_mhalo200)

        #frac_baryon = 0.155 # universal baryon fraction
        #Phi =  hist_mstellar/ hist_mhalo200 / frac_baryon
        
        """
        Plot comulative histogram 
        """
        
        ax1.plot(mhalo_bin_centeres, hist_comulative_Mhalo200, label = name_simulation, c=colors[k],lw=2 )
        #ax1.plot(mhalo_bin_centeres, hist_comulative_Mhalo200, label = simulation, c=colors[k],lw=2 )
        ax1.set_xlabel(r'$Halo \, Mass \, [M_{200c}][M_{Sun}]$', fontsize = 14)
        ax1.set_ylabel(r'#$Halo(200c)$', fontsize = 14)
        # Finish plot
        ax1.legend(loc=1, fontsize=14, frameon=False)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(mhalomin, mhalomax)
#       ax1.text(5e11,5e-3, tex, fontsize=14)
        



        """
        ax1.plot(mhalo_bin_centeres, Phi, label = simulation, c=colors[k],lw=2 )
        ax1.set_xlabel(r'$Halo \, Mass \, [M_{200c}][M_{Sun}]$', fontsize = 14)
        ax1.set_ylabel(r'$M_{Stellar} / M_{Halo_{200c}} (\Omega_b/\Omega_m)^{-1}$', fontsize = 14)
        # Finish plot
        ax1.legend(loc=1, fontsize=14, frameon=False)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(mhalomin, mhalomax)
#       ax1.text(5e11,5e-3, tex, fontsize=14)
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
        fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/comulative_z.%s' %( ending))
    

        
