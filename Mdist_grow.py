from PYmodule import *


log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']

pres = ['./data/1e9','./data/1e10','./data/1e11','./data/1e12','./data/1e13']


print('z=6 to 4, t in Myr: ', (t_from_z(4)-t_from_z(6))/Myr)

M_grow_ratio = 5.
f_lambda = np.log(M_grow_ratio) * t_Edd / (t_from_z(4)-t_from_z(6))
print(f_lambda)

exit(0)

i_bsm = 0
iM = 2

T_tell = 8000
eta = 0.3

T=ascii.read(pres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
T['Mstar0'] = np.zeros(len(T))
T['Mstar_z'] = np.zeros(len(T))
T['Lbol_z'] = np.zeros(len(T))
T['M1450_z'] = np.zeros(len(T))
for i in range(len(T)):
    T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
T['on'] = np.ones(len(T))
T['Edd_ratio'] = np.ones(len(T))
# print("max T['Edd_ratio']=",np.max(T['Edd_ratio']))


abin_mf =  np.logspace(2,12,num=40) # default endpoint=True
wid_mf = abin_mf[1:]-abin_mf[:-1]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print(dlog10M)

for z in [4,6]:
    f_duty = 0.5
    mu_fit = .3
    sigma_fit = .15
    
    N_concatenate = int(1e0)
    # if z==6:
    #     N_concatenate = int(1e4)
    # else:
    #     N_concatenate = int(1e0)
    h0_mf = np.zeros(len(abin_mf)-1); h1_mf = np.zeros(len(abin_mf)-1); h2_mf = np.zeros(len(abin_mf)-1)

    for i_concatenate in range(N_concatenate):
        # T = T[T['z_col']>z]
        T['Edd_ratio'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=len(T)) # scatter=0.1dex; center=scale
        for i in range(len(T)):
            # T['Edd_ratio'][i] = .3
            T['Mstar_z'][i] = T['Mstar0'][i] * np.exp( (t_from_z(z)-t_from_z(T['z_col'][i])) / t_Edd * f_duty* T['Edd_ratio'][i] )
        # print(np.argmax(T['Mstar_z'])," max Mstar:", np.max(T['Mstar_z']), "corresponding Edding_ratio",T['Edd_ratio'][np.argmax(T['Mstar_z'])], "t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Mstar_z'])]))/Myr)

        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]

        hist0_mf, bin_edges = np.histogram(T_H2['Mstar_z'],bins=abin_mf,density=False)
        hist1_mf, bin_edges = np.histogram(T_isofail['Mstar_z'],bins=abin_mf,density=False)
        hist2_mf, bin_edges = np.histogram(T_isoOK['Mstar_z'],bins=abin_mf,density=False)
        h0_mf += hist0_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
        h1_mf += hist1_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
        h2_mf += hist2_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
    

    fig, ax = plt.subplots(1,2,figsize=(20,10),dpi=400)
    ax[0].plot(abin_mf, MF(abin_mf,z),label='Willott 2010 or extrapolate')
    x = (abin_mf[:-1]+abin_mf[1:])/2.   
    ax[0].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    ax[0].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
    ax[0].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
    ax[0].tick_params(labelsize=fstick)
    ax[0].set_xlabel(r'$\mathrm{M_{\bullet}}$',fontsize=fslabel)
    ax[0].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    ax[0].set_xscale('log'); ax[0].set_yscale('log')
    ax[0].set_title('z='+str(int(z)),fontsize=fslabel)
    ax[0].set_xlim(1e7,10**10.5);  ax[0].set_ylim(1.e-9,1e-4)
    ax[0].grid(True)
    ax[0].legend(fontsize=fslegend,loc='best')
    ax[1].plot(abin_mf, MF(abin_mf,z),label='Willott 2010 or extrapolate')
    ax[1].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    ax[1].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
    ax[1].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
    ax[1].tick_params(labelsize=fstick)
    ax[1].set_xlabel(r'$\mathrm{M_{\bullet}}$',fontsize=fslabel)
    ax[1].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    ax[1].set_xscale('log'); ax[1].set_yscale('log')
    ax[1].set_title('z='+str(int(z)),fontsize=fslabel)
    ax[1].set_xlim(abin_mf[0],abin_mf[1]); ax[1].set_ylim()
    ax[1].grid(True)
    ax[1].legend(fontsize=fslegend,loc='best')
    plt.savefig(figpre+'MF_z'+str(int(z))+'eta'+str(int(10*eta))+'f'+str(f_duty)+'mu'+str(mu_fit)+'sigma'+str(sigma_fit)+'N'+str(int(np.log10(N_concatenate)))+'.png')
    ascii.write(Table([(abin_mf[:-1]+abin_mf[1:])/2.,h0_mf,h1_mf,h2_mf, ],names=['bin_center','hist_0','hist_1','hist_2']), datapre+'histBHmass_z'+str(int(z))+'eta'+str(int(10*eta))+'f'+str(f_duty)+'mu'+str(mu_fit)+'sigma'+str(sigma_fit)+'N'+str(int(np.log10(N_concatenate))),overwrite=True)