import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import illustris_python as il
import h5py

parttype_stars = 4

if __name__=='__main__':
    
    # information about simulation
    ubications = ['Illustris','IllustrisTNG','IllustrisTNG','IllustrisTNG']
    simulations = ['L75n1820FP','L75n1820TNG','L205n2500TNG','L205n1250TNG']
    #simulation=   [Illustris, TNG100-1, TNG300-1, TNG100-2]
    snapnums = [135,99,99,99]  # for z = 0
    colors = ['r','b','m','g']

    # Asignar mass stellar in bins
    nbins = 25
    mstarmin = 2e8  #Msun
    mstarmax = 2e12 #Msun
    mstar_bin_edges = np.logspace(np.log10(mstarmin), np.log10(mstarmax), nbins+1)
    mstar_bin_centeres = np.sqrt(mstar_bin_edges[:-1] * mstar_bin_edges[1:]) # geometric mean
    d_mstar = mstar_bin_edges[1:]-mstar_bin_edges[:-1] # bin Delta Mass
    d_log_mstar = np.log10(mstar_bin_edges[1:])-np.log10(mstar_bin_edges[:-1]) # bin Delta log Mass
    
    # create figure
    fig = plt.figure(figsize=(14,10))
#    fig.add_subplot(111)
    ax1 = fig.add_subplot(221)
    #    ax1,ax2 = fig.add_subplot(111)
    #    fig, axes = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
    
#    fig  = plt.subplots(nrows=1, ncols= 2, figsize=(14,6))
#    ax1,ax2 = axes.flatten()
#    fig = plt.figure(figsize=(6,6))
#    fig.add_subplot(111)
    ax2 = fig.add_subplot(222)
       
#    ax1,ax2 = fig.add_subplot(111)

    for k in range(len(simulations)):
        ubication = ubications[k]
        simulation = simulations[k]
        snapnum = snapnums[k]
        
#        """        
        # selecc the folder with data
        if ubication == 'Illustris':
            basedir = '/n/hernquistfs1/%s/Runs/%s/output' % (ubication, simulation)
        elif ubication == 'IllustrisTNG':
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

        ####==== Condition for save data about mass rTNG300 ====####
    

        if simulation == 'L205n1250TNG':
        
            #### ==== 
            #####-------------------------------------------------#########    
            # load stellar Mass (SubhaloMassInHalfRadType)
            mstar_in_rad = il.groupcat.loadSubhalos(
                basedir, snapnum, fields=['SubhaloMassInHalfRadType'])[:, parttype_stars]
            # covert  units to Msun
            mstar_in_rad_Msun_TNG_2  =mstar_in_rad * 1e10 /h

            ### === Charge mass limulations
            
            mstar_in_rad_Msun_TNG_1 = np.loadtxt('mstar_in_rad_Msun_L75n1820TNG.txt') # Simulation TNG100 -1
            mstar_in_rad_Msun_TNG300 = np.loadtxt('mstar_in_rad_Msun_L205n2500TNG.txt') # simulation TNG300

            Mass_in_rad_rTNG300 = mstar_in_rad_Msun_TNG300 * (mstar_in_rad_Msun_TNG_1/mstar_in_rad_Msun_TNG_2)

            #####-------------------------------------------------#########    
            # load stellar Mass (SubhaloMassType)
            mstar = il.groupcat.loadSubhalos(
                basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
            # covert  units to Msun
            mstar_Msun_TNG_2 = mstar * 1e10 /h
            
            ### === Charge mass limulations
            
            mstar_Msun_TNG_1 = np.loadtxt('mstar_Msun_L75n1820TNG.txt') # Simulation TNG100 -1
            mstar_Msun_TNG300 = np.loadtxt('mstar_Msun_L205n2500TNG.txt') # simulation TNG300            
            
            Mass_rTNG300 = mstar_Msun_TNG300 * (mstar_Msun_TNG_1/mstar_Msun_TNG_2)
            
            #######-----------------------------------------#####
            # Calculate mass funtion (dN/d(Log10(Mass)))
            # N = numers of galaxies 
            tex = r'$z\sim 0$'
            #===========  Mass  =============
            hist, ignored = np.histogram(mstar_Msun,bins =mstar_bin_edges) ## dN = hist
            Phi1 = hist / (d_log_mstar * volume_Mpc)
            ax1.plot(mstar_bin_centeres, Phi1, label = simulation, c=colors[k],lw=2 )
            ax1.set_xlabel(r'$Stellar \ Mass\ [M_{Sun}]$', fontsize = 14)
            ax1.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)
            # Finish plot
            ax1.legend(loc=1, fontsize=14, frameon=False)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
            ax1.set_xlim(mstarmin, mstarmax)
            #       ax1.text(5e11,5e-3, tex, fontsize=14)

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
        
            print('Finished for %s.' %(simulation))



            #    ax.set_xlabel(r'$Stellar \ Mass\ [M_{Sun}]$', fontsize = 14)
            #    ax.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)

            for ending in ['pdf','png']:
                fig.suptitle('Mass function galaxies to $z\sim %s$'%(int(z)), fontsize=20)
                fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/mass_funtion_nbins_%d_z%s.%s' %(nbins, int(z) , ending))
    


        #### ==== Finishid if ===========####
        
        
    else:
        
        #### ==== opening ==========###
        #####-------------------------------------------------#########    
        # load stellar Mass (SubhaloMassInHalfRadType)
        mstar_in_rad = il.groupcat.loadSubhalos(
            basedir, snapnum, fields=['SubhaloMassInHalfRadType'])[:, parttype_stars]
        # covert  units to Msun
        mstar_in_rad_Msun = mstar_in_rad * 1e10 /h

        #### === save mass for cada simulation
        np.savetxt('mstar_in_rad_Msun_%s.txt'%(simulation))

        #####-------------------------------------------------#########    
        # load stellar Mass (SubhaloMassType)
        mstar = il.groupcat.loadSubhalos(
            basedir, snapnum, fields=['SubhaloMassType'])[:, parttype_stars]
        # covert  units to Msun
        mstar_Msun  =mstar * 1e10 /h
        
        #### === save mass for cada simulation
        np.savetxt('mstar_Msun_%s.txt'%(simulation))

        # Calculate mass funtion (dN/d(Log10(Mass)))
        # N = numers of galaxies 
        
        #===========  Mass  =============
        hist, ignored = np.histogram(mstar_Msun,bins =mstar_bin_edges) ## dN = hist
        Phi1 = hist / (d_log_mstar * volume_Mpc)
        ax1.plot(mstar_bin_centeres, Phi1, label = simulation, c=colors[k],lw=2 )
        ax1.set_xlabel(r'$Stellar \ Mass\ [M_{Sun}]$', fontsize = 14)
        ax1.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)
        # Finish plot
        ax1.legend(loc=1, fontsize=14, frameon=False)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(mstarmin, mstarmax)
        #       ax1.text(5e11,5e-3, tex, fontsize=14)

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
        
        print('Finished for %s.' %(simulation))


        
        #    ax.set_xlabel(r'$Stellar \ Mass\ [M_{Sun}]$', fontsize = 14)
        #    ax.set_ylabel(r'$\Phi\ [Mpc^{-3}\ dex^{-1}]$', fontsize = 14)

        for ending in ['pdf','png']:
            fig.suptitle('Mass function galaxies to $z\sim %s$'%(int(z)), fontsize=20)
            fig.savefig('/n/home10/dmontenegro/PythonModules/graficas/mass_funtion_nbins_%d_z%s.%s' %(nbins, int(z) , ending))
    


        #### ==== Finishid else ===========####

