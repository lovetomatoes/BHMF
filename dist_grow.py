from PYmodule import *

log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {'4':'Akiyama_18','5':'McGreer_18','6':'Matsuoka_18'}
pres = [datapre+'1e9',datapre+'1e10',datapre+'1e11',datapre+'1e12',datapre+'1e13']
figprefix = './figs/'

print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

i_bsm = int(1)
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
T['type'] = np.ones(len(T))
T['Edd_ratio'] = np.ones(len(T))

# LF
abin_lf = np.linspace(-31,-7,num=30,endpoint=True)  # endpoint=True
wid_lf = abs(abin_lf[1:]-abin_lf[:-1])
dmag = abin_lf[1]-abin_lf[0] # print(dmag)
# MF
abin_mf =  np.logspace(2,12,num=40) # default endpoint=True
wid_mf = abin_mf[1:]-abin_mf[:-1]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print(dlog10M)
flambda = .19
print('flambda= ', flambda)

#for f_duty in np.linspace(0.4,0.7,num=2):
for f_duty in [0.5]:
    for z in [6]:
        mu_fit = flambda/f_duty
        sigma_fit = .12
        N_concatenate = int(1e2)
        
        T_long = Table(
            [[], [], [], [], [], []],
            names=('flambda','M1450_z','Mstar_z','type','z_col','Mstar0')
        )

        h0_lf = np.zeros(len(abin_lf)-1); h1_lf = np.zeros(len(abin_lf)-1); h2_lf = np.zeros(len(abin_lf)-1)
        h0_mf = np.zeros(len(abin_mf)-1); h1_mf = np.zeros(len(abin_mf)-1); h2_mf = np.zeros(len(abin_mf)-1)
        print(time.ctime())
        for i_concatenate in range(N_concatenate):
            print(i_concatenate)
            # T = T[T['z_col']>z]
            T['Edd_ratio'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=len(T)) # scatter=0.1dex; center=scale
            #T['Edd_ratio'] = np.ones(len(T))*.38
            for i in range(len(T)):
               # T['Edd_ratio'][i] = max(1.,T['Edd_ratio'][i])
                T['Mstar_z'][i] = T['Mstar0'][i] * np.exp( (t_from_z(z)-t_from_z(T['z_col'][i])) / t_Edd * f_duty* T['Edd_ratio'][i] )
                # T['Lbol_z'][i] = L_M(T['Mstar_z'][i], T['Edd_ratio'][i]) # 或使用rvs随机生成一个lbd
                T['Lbol_z'][i] = L_M(T['Mstar_z'][i], lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit)) # 使用rvs随机生成一个lbd 未用本身的Edd_ratio
                T['M1450_z'][i] = M1450_Lbol(T['Lbol_z'][i])
                if T['Tg_max'][i] <= T_tell:
                    T['type'][i] = 1
                elif T['iso_col'][i] == 0:
                    T['type'][i] = 2
                else:
                    T['type'][i] = 3
            T_long = vstack( [T_long,Table([ np.ones(len(T))*flambda,T['M1450_z'],T['Mstar_z'],T['type'],T['z_col'],T['Mstar0']],names=('flambda','M1450_z','Mstar_z','type','z_col','Mstar0')) ])
            # print(T_long)
            # print(np.argmax(T['Edd_ratio'])," max Eddingratio:", np.max(T['Edd_ratio']), "corresponding Mstar",T['Mstar_z'][np.argmax(T['Edd_ratio'])],"t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Edd_ratio'])]))/Myr)
            # print(np.argmax(T['Mstar_z'])," max Mstar:", np.max(T['Mstar_z']), "corresponding Edding_ratio",T['Edd_ratio'][np.argmax(T['Mstar_z'])], "t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Mstar_z'])]))/Myr)
            # print(np.max(T['Lbol_z']))
            # print(np.min(T['M1450_z']))
            T_H2 = T[T['Tg_max']<=T_tell] 
            T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
            T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
            
        # LF
            hist0_lf, bin_edges = np.histogram(T_H2['M1450_z'],bins=abin_lf,density=False)
            hist1_lf, bin_edges = np.histogram(T_isofail['M1450_z'],bins=abin_lf,density=False)
            hist2_lf, bin_edges = np.histogram(T_isoOK['M1450_z'],bins=abin_lf,density=False)
            h0_lf += hist0_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty
            h1_lf += hist1_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty
            h2_lf += hist2_lf*n_base[iM]/(1e4*N_concatenate)/dmag*f_duty
        # MF
            hist0_mf, bin_edges = np.histogram(T_H2['Mstar_z'],bins=abin_mf,density=False)
            hist1_mf, bin_edges = np.histogram(T_isofail['Mstar_z'],bins=abin_mf,density=False)
            hist2_mf, bin_edges = np.histogram(T_isoOK['Mstar_z'],bins=abin_mf,density=False)
            h0_mf += hist0_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
            h1_mf += hist1_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
            h2_mf += hist2_mf*n_base[iM]/(1e4*N_concatenate)/dlog10M*f_duty
    
        print("Nconcatenate",N_concatenate,"time",time.ctime())

        fig, ax = plt.subplots(2,2,figsize=(20,24),dpi=400)
        ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(int(z))])
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
        ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(int(z))])
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
        ax[1,0].set_title('z='+str(int(z)),fontsize=fslabel)
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
        ax[1,1].set_title('z='+str(int(z)),fontsize=fslabel)
        ax[1,1].set_xlim(abin_mf[0],abin_mf[-1]); ax[1,1].set_ylim(bottom=1e-10)
        ax[1,1].grid(True)
        ax[1,1].legend(fontsize=fslegend,loc='best')

        plt.savefig(z6figpre+'/f'+str(int(flambda*100))+'s'+str(int(sigma_fit*100))+'N'+str(int(np.log10(N_concatenate)))+'alpha1.png')
        ascii.write(T_long, z6datapre+'LF_z'+str(int(z))+'N'+str(int(np.log10(N_concatenate)))+'bsm'+str(i_bsm)+'alpha1',overwrite=True)