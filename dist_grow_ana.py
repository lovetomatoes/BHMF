from PYmodule import *
# wli try LF ana form; may from dist_grow_MBH
log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {'4':'Akiyama_18','5':'McGreer_18','6':'Matsuoka_18'}
pres = [datapre+'1e9',datapre+'1e10',datapre+'1e11',datapre+'1e12',datapre+'1e13']
figprefix = './figs/'

print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

iM = 2

T_tell = 8000
eta = 0.3

Ts = [] # bsm=0,1 two files
for i_bsm in range(2):
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
    Ts.append(T)

# LF
abin_lf = np.linspace(-31,-7,num=30,endpoint=True)  # endpoint=True
wid_lf = abs(abin_lf[1:]-abin_lf[:-1])
dmag = abin_lf[1]-abin_lf[0] # print(dmag)

# luminosity bin same as Matsuoka18 setting
bin_cen = [-29,-27.5,-26.75,-26.25,-25.75,-25.25,-24.75,-24.25,-23.75,-23.25,-22.75,-22]
bin_wid = [2, 1, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1]
Nbin = len(bin_cen)
bin_edg = np.zeros(Nbin+1)
for i in range(Nbin):
    bin_edg[i] = bin_cen[i] - bin_wid[i]/2.
bin_edg[-1] =  bin_cen[i] + bin_wid[i]/2.
# print(bin_edg)

flambda = .19 # .18 .20
print('flambda= ', flambda)
z = 6
k_M = 1.1

for f_duty in [0.5]: # .6 .4 
    for sigma_fit in [.12]: # .10  .14
        mu_fit = flambda/f_duty
        dP_0 = 1e-9*np.ones(Nbin)
        dP = np.zeros(Nbin)
        Phi = np.zeros(Nbin)

        for ibin in range(1): # Nbin
            for i_bsm in range(1): # 2 
                jcount = 0
                while (abs(dP[ibin]/dP_0[ibin]-1.)>1.0e-4):
                    jcount +=1 
                    T = Ts[i_bsm]
                    # T = T[T['z_col']>z]
                    dP_0[ibin] = dP[ibin]
                    for i in range(1): # len(T):
                        M_BH = T['Mstar0'][i]*np.exp((t_from_z(z)-t_from_z(T['z_col'][i]))/t_Edd*f_duty*mu_fit)
                        print('%5.1e'%(T['Mstar0'][i]*np.exp((t_from_z(z)-t_from_z(T['z_col'][i]))/t_Edd*f_duty*mu_fit)))
                        print('z_col = %3.1f'%(T['z_col'][i]))
            
                        x0 = kernel_MBH(M_BH/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit*np.log(10))
                        x1 = kernel_MBH(k_M*M_BH/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit*np.log(10))
                        dP_MBH = math.erf(x1) - math.erf(x0)

                        x0 = kernel_M1450(bin_edg[ibin+1], M_BH*math.sqrt(k_M), mu_fit, sigma_fit*np.log(10))
                        x1 = kernel_M1450(bin_edg[ibin], M_BH*math.sqrt(k_M), mu_fit, sigma_fit*np.log(10))
                        dP_M1450 = math.erf(x1) - math.erf(x0)
                        dP_M1450 = np.exp(-5*jcount)
                        dP_MBH = 1

                        dP[ibin] += dP_MBH*dP_M1450
                    print('dP[ibin]=%5.3e'%dP[ibin])
                    M_BH *= k_M
                    # print('%5.1e'%x0,'%5.1e'%x1,'%5.1e'%dP_M1450)
                    print('M_BH=%5.1e'%M_BH)
                # hist Phi compared w/ LF
                Phi[ibin] += dP[ibin]*n_base[iM]*f_bsm[i_bsm]/1e4*f_duty*1e9
                print('%5.1e'%x0,'%5.1e'%x1,'%5.1e'%dP_M1450)
                print('jcount=',jcount, 'Phi[ibin]=%5.3e'%Phi[ibin])
                # t3 = time.time()
                # print("t3-t2=",(t3-t2)/3600, "hrs")
        exit(0)
        fig, ax = plt.subplots(2,2,figsize=(20,24),dpi=400)
        ax[0,0].bar(bin_cen, Phi,width=bin_wid,color='grey',alpha=0.5,label='Model')
        ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(int(z))])
        ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9/(1-f_obsc_const),label='CO=5 corrected')
        ax[0,0].plot(abin_lf, LF_M1450(abin_lf,z)*1e9*corr_U14D20(abin_lf),label='DO corrected')
        ax[0,0].tick_params(labelsize=fstick)
        ax[0,0].set_xlabel(r'$\mathrm{mag}$',fontsize=fslabel)
        ax[0,0].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Gpc^{-3}mag^{-1}]}$',fontsize=fslabel)
        ax[0,0].set_yscale('log')
        ax[0,0].set_title('z='+str(int(z)),fontsize=fslabel)
        ax[0,0].set_xlim(bin_edg[-1],bin_edg[0]);  ax[0,0].set_ylim(1.e-2,1e3)
        ax[0,0].grid(True)
        ax[0,0].legend(fontsize=fslegend,loc='best')

        # whole M1450 range 
        x = (abin_lf[:-1]+abin_lf[1:])/2.
        ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9,label=lfnames[str(int(z))])
        ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9/(1-f_obsc_const),label='CO=5 corrected')
        ax[0,1].plot(abin_lf, LF_M1450(abin_lf,z)*1e9*corr_U14D20(abin_lf),label='DO corrected')
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
        # ax[1,0].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
        # ax[1,0].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
        # ax[1,0].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
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
        # ax[1,1].bar(x, h0_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
        # ax[1,1].bar(x, h1_mf,width=wid_mf,bottom=h0_mf,color='C'+str(1),alpha=0.5,label=typenames[1])
        # ax[1,1].bar(x, h2_mf,width=wid_mf,bottom=(h0_mf+h1_mf),color='C'+str(2),alpha=0.5,label=typenames[2])
        ax[1,1].tick_params(labelsize=fstick)
        ax[1,1].set_xlabel(r'$\mathrm{M_{\bullet}}$',fontsize=fslabel)
        ax[1,1].set_ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
        ax[1,1].set_xscale('log'); ax[1,1].set_yscale('log')
        ax[1,1].set_title('z='+str(int(z)),fontsize=fslabel)
        ax[1,1].set_xlim(abin_mf[0],abin_mf[-1]); ax[1,1].set_ylim(bottom=1e-10)
        ax[1,1].grid(True)
        ax[1,1].legend(fontsize=fslegend,loc='best')

        plt.savefig(z6figpre+'/anafl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm01alpha1.png')
        T = Table(
            [bin_cen, Phi],
            names=('bin_cen','Phi')
        )
        ascii.write(T, z6datapre+'anaPhi_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm01alpha1',formats={'bin_cen':'%6.2f','Phi':'4.2e'},overwrite=True)
        # ascii.write(T_long, z6datapre+'LF_z'+str(int(z))+'N'+str(int(np.log10(N_concatenate)))+'bsm'+str(i_bsm)+'alpha1',overwrite=True)