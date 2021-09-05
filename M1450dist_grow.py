from PYmodule import *

log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']

pres = ['./data/1e9','./data/1e10','./data/1e11','./data/1e12','./data/1e13']

print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

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


abin_lf = np.linspace(-31,-7,num=30,endpoint=True)  # endpoint=True
wid_lf = abs(abin_lf[1:]-abin_lf[:-1])
dmag = abin_lf[1]-abin_lf[0] # print(dmag)


for z in [6,5,4]:
    f_duty = 0.5
    mu_fit = .3
    sigma_fit = .15

    N_concatenate = int(1e0)
    # N_concatenate = (z==6)?int(1e4):int(1e2)

    h0_lf = np.zeros(len(abin_lf)-1); h1_lf = np.zeros(len(abin_lf)-1); h2_lf = np.zeros(len(abin_lf)-1)

    for i_concatenate in range(N_concatenate):
        # T = T[T['z_col']>z]
        T['Edd_ratio'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=len(T)) # scatter=0.1dex; center=scale
        for i in range(len(T)):
            # T['Edd_ratio'][i] = .3
            T['Mstar_z'][i] = T['Mstar0'][i] * np.exp( (t_from_z(z)-t_from_z(T['z_col'][i])) / t_Edd * f_duty* T['Edd_ratio'][i] )
            # T['Lbol_z'][i] = L_M(T['Mstar_z'][i], T['Edd_ratio'][i]) # 或使用rvs随机生成一个lbd
            T['Lbol_z'][i] = L_M(T['Mstar_z'][i], lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit)) # 使用rvs随机生成一个lbd 未用本身的Edd_ratio
            T['M1450_z'][i] = M1450_Lbol(T['Lbol_z'][i])

        # print(np.argmax(T['Edd_ratio'])," max Eddingratio:", np.max(T['Edd_ratio']), "corresponding Mstar",T['Mstar_z'][np.argmax(T['Edd_ratio'])],"t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Edd_ratio'])]))/Myr)
        # print(np.argmax(T['Mstar_z'])," max Mstar:", np.max(T['Mstar_z']), "corresponding Edding_ratio",T['Edd_ratio'][np.argmax(T['Mstar_z'])], "t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Mstar_z'])]))/Myr)
        # print(np.max(T['Lbol_z']))
        # print(np.min(T['M1450_z']))
        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
        hist0_lf, bin_edges = np.histogram(T_H2['M1450_z'],bins=abin_lf,density=False)
        hist1_lf, bin_edges = np.histogram(T_isofail['M1450_z'],bins=abin_lf,density=False)
        hist2_lf, bin_edges = np.histogram(T_isoOK['M1450_z'],bins=abin_lf,density=False)

        h0_lf += hist0_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty
        h1_lf += hist1_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty
        h2_lf += hist2_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty

    # print(h0_lf)
    # print(h1_lf)
    # print(h2_lf)

    fig, ax = plt.subplots(2,2,figsize=(20,24),dpi=400)
    ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9)
    x = (abin_lf[:-1]+abin_lf[1:])/2.
    M_1 = -27.2
    M_2 = -20.7
    t = (x - M_1) / (M_2 - M_1)
    h0_lf *= 1e9/(2*(1-t)+3*t)
    h1_lf *= 1e9/(2*(1-t)+3*t)
    h2_lf *= 1e9/(2*(1-t)+3*t)
    ax[0,0].bar(x, h0_lf,width=wid_lf,color='C'+str(0),alpha=0.5,label=typenames[0])
    ax[0,0].bar(x, h1_lf,width=wid_lf,bottom=h0_lf,color='C'+str(1),alpha=0.5,label=typenames[1])
    ax[0,0].bar(x, h2_lf,width=wid_lf,bottom=(h0_lf+h1_lf),color='C'+str(2),alpha=0.5,label=typenames[2])
    ax[0,0].tick_params(labelsize=fstick)
    ax[0,0].set_xlabel(r'$\mathrm{mag}$',fontsize=fslabel)
    ax[0,0].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Gpc^{-3}mag^{-1}]}$',fontsize=fslabel)
    ax[0,0].set_yscale('log')
    ax[0,0].set_title('z='+str(int(z)),fontsize=fslabel)
    ax[0,0].set_xlim(-22,-29);  ax[0,0].set_ylim(1.e-2,1e3)
    ax[0,0].grid(True)
    ax[0,0].legend(fontsize=fslegend,loc='best')
    # whole M1450 range 
    ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9)
    ax[0,1].bar(x, h0_lf,width=wid_lf,color='C'+str(0),alpha=0.5,label=typenames[0])
    ax[0,1].bar(x, h1_lf,width=wid_lf,bottom=h0_lf,color='C'+str(1),alpha=0.5,label=typenames[1])
    ax[0,1].bar(x, h2_lf,width=wid_lf,bottom=(h0_lf+h1_lf),color='C'+str(2),alpha=0.5,label=typenames[2])
    ax[0,1].tick_params(labelsize=fstick)
    ax[0,1].set_xlabel(r'$\mathrm{mag}$',fontsize=fslabel)
    ax[0,1].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Gpc^{-3}mag^{-1}]}$',fontsize=fslabel)
    ax[0,1].set_yscale('log')
    ax[0,1].set_title('z='+str(int(z)),fontsize=fslabel)
    ax[0,1].set_xlim(abin_lf[-1],abin_lf[0]); ax[0,1].set_ylim()
    ax[0,1].grid(True)
    ax[0,1].legend(fontsize=fslegend,loc='best')
    plt.savefig(figpre+'LF_M1450_z'+str(int(z))+'eta'+str(int(10*eta))+'f'+str(f_duty)+'mu'+str(mu_fit)+'sigma'+str(sigma_fit)+'N'+str(int(np.log10(N_concatenate)))+'.png')
    ascii.write(Table([x,h0_lf,h1_lf,h2_lf],names=['bin_center','hist_0','hist_1','hist_2']), datapre+'histM1450_z'+str(int(z))+'eta'+str(int(10*eta))+'f'+str(f_duty)+'mu'+str(mu_fit)+'sigma'+str(sigma_fit)+'N'+str(int(np.log10(N_concatenate))),overwrite=True)
