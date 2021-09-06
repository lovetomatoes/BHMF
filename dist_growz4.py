from PYmodule import *

log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {'4':'Akiyama_18','5':'McGreer_18','6':'Matsuoka_18'}


print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

i_bsm = 0
iM = 2

T_tell = 8000
eta = 0.3

T = ascii.read(z6datapre+'LF_z6N3', guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
Nsite = len(T) # N3 <-> 1e7
T['Lbol_z'] = np.zeros(Nsite)
T['M1450_z'] = np.zeros(Nsite)
T['Mstar_z6'] = T['Mstar_z'] # BH mass grown to z=6

# LF
abin_lf = np.linspace(-31,-7,num=100,endpoint=True)  # 100 grids, close to Aki 2018
# abin_lf =np.linspace(-29, -21.75,num=30) # Aki
wid_lf = abs(abin_lf[1:]-abin_lf[:-1])
dmag = abin_lf[1]-abin_lf[0] # print(dmag)

# MF
abin_mf =  np.logspace(2,12,num=40) # default endpoint=True
wid_mf = abin_mf[1:]-abin_mf[:-1]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print(dlog10M)

for flambda in [.11, .12, .13]:
    print('flambda= ', flambda)
    f_duty = 0.5
    z = int(4)
    mu_fit = flambda/f_duty
    for sigma_fit in [.18]:
        print('sigma= ', sigma_fit)
        h0_lf = np.zeros(len(abin_lf)-1); h1_lf = np.zeros(len(abin_lf)-1); h2_lf = np.zeros(len(abin_lf)-1)
        h0_mf = np.zeros(len(abin_mf)-1); h1_mf = np.zeros(len(abin_mf)-1); h2_mf = np.zeros(len(abin_mf)-1)

        # T = T[T['z_col']>z] # all satisfying
        T['Edd_ratio'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=Nsite) # scatter=0.1dex; center=scale
        # T['Edd_ratio'] = np.ones(Nsite)*mu_fit # fix Eddington ratio
        for i in range(Nsite):
            T['Mstar_z'][i]  = T['Mstar_z6'][i]*np.exp( (t_from_z(z)-t_from_z(6)) / t_Edd * f_duty* T['Edd_ratio'][i] )
            # T['Lbol_z'][i] = L_M(T['Mstar_z'][i], T['Edd_ratio'][i]) # 用本身的Edd_ratio
            T['Lbol_z'][i] = L_M(T['Mstar_z'][i], lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit)) # wli 使用rvs随机生成一个lbd
            T['M1450_z'][i] = M1450_Lbol(T['Lbol_z'][i])

    # ------------------------  CUT  ---------------------------
    #  constrain M_BH < 1e10 Msun
        T_ = T[T['Mstar_z']<1e10]
    # ----------------------------------------------------------

        T_H2 = T_[T_['type']==1]
        T_isofail = T_[T_['type']==2]
        T_isoOK =  T_[T_['type']==3]
 
    # LF
        hist0_lf, bin_edges = np.histogram(T_H2['M1450_z'],bins=abin_lf,density=False)
        hist1_lf, bin_edges = np.histogram(T_isofail['M1450_z'],bins=abin_lf,density=False)
        hist2_lf, bin_edges = np.histogram(T_isoOK['M1450_z'],bins=abin_lf,density=False)
        h0_lf = hist0_lf*n_base[iM]/float(Nsite)/dmag*f_duty*1e9
        h1_lf = hist1_lf*n_base[iM]/float(Nsite)/dmag*f_duty*1e9
        h2_lf = hist2_lf*n_base[iM]/float(Nsite)/dmag*f_duty*1e9

    # MF
        hist0_mf, bin_edges = np.histogram(T_H2['Mstar_z'],bins=abin_mf,density=False)
        hist1_mf, bin_edges = np.histogram(T_isofail['Mstar_z'],bins=abin_mf,density=False)
        hist2_mf, bin_edges = np.histogram(T_isoOK['Mstar_z'],bins=abin_mf,density=False)
        h0_mf = hist0_mf*n_base[iM]/float(Nsite)/dlog10M*f_duty
        h1_mf = hist1_mf*n_base[iM]/float(Nsite)/dlog10M*f_duty
        h2_mf = hist2_mf*n_base[iM]/float(Nsite)/dlog10M*f_duty
    

        fig, ax = plt.subplots(2,2,figsize=(20,24),dpi=400)
        ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(z)])
        x = (abin_lf[:-1]+abin_lf[1:])/2.

        # M_1 = -27.2
        # M_2 = -20.7
        # t = (x - M_1) / (M_2 - M_1)
        # h0_lf *= 1/(2*(1-t)+3*t)
        # h1_lf *= 1/(2*(1-t)+3*t)
        # h2_lf *= 1/(2*(1-t)+3*t)

        ax[0,0].bar(x, h0_lf,width=wid_lf,color='C'+str(0),alpha=0.5,label=typenames[0])
        ax[0,0].bar(x, h1_lf,width=wid_lf,bottom=h0_lf,color='C'+str(1),alpha=0.5,label=typenames[1])
        ax[0,0].bar(x, h2_lf,width=wid_lf,bottom=(h0_lf+h1_lf),color='C'+str(2),alpha=0.5,label=typenames[2])
        ax[0,0].tick_params(labelsize=fstick)
        ax[0,0].set_xlabel(r'$\mathrm{mag}$',fontsize=fslabel)
        ax[0,0].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Gpc^{-3}mag^{-1}]}$',fontsize=fslabel)
        ax[0,0].set_yscale('log')
        ax[0,0].set_title('z='+str(z),fontsize=fslabel)
        ax[0,0].set_xlim(-22,-29);  ax[0,0].set_ylim(1.e-2,1e4)
        ax[0,0].grid(True)
        ax[0,0].legend(fontsize=fslegend,loc='best')
        # whole M1450 range 
        ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(z)])
        ax[0,1].bar(x, h0_lf,width=wid_lf,color='C'+str(0),alpha=0.5,label=typenames[0])
        ax[0,1].bar(x, h1_lf,width=wid_lf,bottom=h0_lf,color='C'+str(1),alpha=0.5,label=typenames[1])
        ax[0,1].bar(x, h2_lf,width=wid_lf,bottom=(h0_lf+h1_lf),color='C'+str(2),alpha=0.5,label=typenames[2])
        ax[0,1].tick_params(labelsize=fstick)
        ax[0,1].set_xlabel(r'$\mathrm{mag}$',fontsize=fslabel)
        ax[0,1].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Gpc^{-3}mag^{-1}]}$',fontsize=fslabel)
        ax[0,1].set_yscale('log')
        ax[0,1].set_title('z='+str(z),fontsize=fslabel)
        ax[0,1].set_xlim(abin_lf[-1],abin_lf[0]); ax[0,1].set_ylim()
        # ax[0,1].set_xlim(-22,-29); ax[0,1].set_ylim() # wli: just for check
        ax[0,1].grid(True)
        ax[0,1].legend(fontsize=fslegend,loc='best')

        x = (abin_mf[:-1]+abin_mf[1:])/2.
        if z==4 or z==6:
            ax[1,0].plot(x, MF(x,z),label='Willott 2010 or extrapolate')
        ax[1,0].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
        ax[1,0].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
        ax[1,0].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
        ax[1,0].tick_params(labelsize=fstick)
        ax[1,0].set_xlabel(r'$\mathrm{M_{\bullet}}$',fontsize=fslabel)
        ax[1,0].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
        ax[1,0].set_xscale('log'); ax[1,0].set_yscale('log')
        ax[1,0].set_title('z='+str(z),fontsize=fslabel)
        ax[1,0].set_xlim(1e7,10**10.5);  ax[1,0].set_ylim(1.e-9,1e-4)
        ax[1,0].grid(True)
        ax[1,0].legend(fontsize=fslegend,loc='best')
        if z==4 or z==6:
            ax[1,1].plot(x, MF(x,z),label='Willott 2010 or extrapolate')
        ax[1,1].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
        ax[1,1].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
        ax[1,1].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
        ax[1,1].tick_params(labelsize=fstick)
        ax[1,1].set_xlabel(r'$\mathrm{M_{\bullet}}$',fontsize=fslabel)
        ax[1,1].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
        ax[1,1].set_xscale('log'); ax[1,1].set_yscale('log')
        ax[1,1].set_title('z='+str(z),fontsize=fslabel)
        ax[1,1].set_xlim(abin_mf[0],abin_mf[-1]); ax[1,1].set_ylim(bottom=1e-10)
        ax[1,1].grid(True)
        ax[1,1].legend(fontsize=fslegend,loc='best')

        plt.savefig(z4figpre+'z'+str(z)+'f'+str(int(flambda*100))+'s'+str(int(sigma_fit*100))+'N'+str(int(np.log10(Nsite)-4))+'.png')
        ascii.write(T, z4datapre+'selfEddLF_z'+str(z)+'f'+str(int(flambda*100))+'s'+str(int(sigma_fit*100))+'N'+str(int(np.log10(Nsite)-4)),overwrite=True)
