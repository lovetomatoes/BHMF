from PYmodule import *

typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']

pres = ['./data/1e9','./data/1e10','./data/1e11','./data/1e12','./data/1e13']

z = 6
i_bsm = 0
iM = 2

T_tell = 8000
eta = 0.3

T=ascii.read(pres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
T['Mstar0'] = np.zeros(len(T))
T['Mstar_z'] = np.zeros(len(T))
T['Lbol_z'] = np.zeros(len(T))
for i in range(len(T)):
    T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
T['on'] = np.ones(len(T))
T['Edd_ratio'] = np.ones(len(T))
# print("max T['Edd_ratio']=",np.max(T['Edd_ratio']))
f_duty = 0.5

print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

abin = np.logspace(38,50,num=20)  # endpoint=True
wid = abin[1:]-abin[:-1]

dlogL = np.log10(abin[1]/abin[0]) # print(dlogL)

h_0 = np.zeros(len(abin)-1); h_1 = np.zeros(len(abin)-1); h_2 = np.zeros(len(abin)-1)
h00 = np.zeros(len(abin)-1); h11 = np.zeros(len(abin)-1); h22 = np.zeros(len(abin)-1)

N_concatenate = int(1e0)
for i_concatenate in range(N_concatenate):
    # T = T[T['z_col']>z]
    mu_fit = .3
    sigma_fit = .15
    T['Edd_ratio'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=len(T)) # scatter=0.1dex; center=scale
    for i in range(len(T)):
        # T['Edd_ratio'][i] = .3
        T['Mstar_z'][i] = T['Mstar0'][i] * np.exp( (t_from_z(z)-t_from_z(T['z_col'][i])) / t_Edd * f_duty* T['Edd_ratio'][i] )
        # T['Lbol_z'][i] = L_M(T['Mstar_z'][i], T['Edd_ratio'][i]) # 或使用rvs随机生成一个lbd
        T['Lbol_z'][i] = L_M(T['Mstar_z'][i], lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit)) # 或使用rvs随机生成一个lbd
    # print(np.argmax(T['Edd_ratio'])," max Eddingratio:", np.max(T['Edd_ratio']), "corresponding Mstar",T['Mstar_z'][np.argmax(T['Edd_ratio'])],"t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Edd_ratio'])]))/Myr)
    # print(np.argmax(T['Mstar_z'])," max Mstar:", np.max(T['Mstar_z']), "corresponding Edding_ratio",T['Edd_ratio'][np.argmax(T['Mstar_z'])], "t_grow Myr",(t_from_z(6)-t_from_z(T['z_col'][np.argmax(T['Mstar_z'])]))/Myr)

    T_H2 = T[T['Tg_max']<=T_tell] 
    T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
    T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
    hist0, bin_edges = np.histogram(T_H2['Lbol_z'],bins=abin,density=False)
    hist1, bin_edges = np.histogram(T_isofail['Lbol_z'],bins=abin,density=False)
    hist2, bin_edges = np.histogram(T_isoOK['Lbol_z'],bins=abin,density=False)

    h_0 += hist0*n_base[iM]/(1e4*N_concatenate)/dlogL
    h_1 += hist1*n_base[iM]/(1e4*N_concatenate)/dlogL
    h_2 += hist2*n_base[iM]/(1e4*N_concatenate)/dlogL

# print(h_0)
# continue
plt.figure(figsize=(10,8),dpi=400)
plt.plot(abin, LF(abin))
plt.bar( (abin[:-1]+abin[1:])/2.,h_0*f_duty,width=wid,color='C'+str(0),alpha=0.5,label=typenames[0])
plt.bar( (abin[:-1]+abin[1:])/2.,h_1*f_duty,width=wid,bottom=h_0*f_duty,color='C'+str(1),alpha=0.5,label=typenames[1])
plt.bar( (abin[:-1]+abin[1:])/2.,h_2*f_duty,width=wid,bottom=(h_0+h_1)*f_duty,color='C'+str(2),alpha=0.5,label=typenames[2])

plt.tick_params(labelsize=fstick)
s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
plt.xlabel(r'$\mathrm{M_\bullet}$'+r' $(M_{\odot})$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{dn/d\logL}$'+' '+r'$\mathrm{[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
plt.xscale('log'); plt.yscale('log')
plt.title('z='+str(int(z)),fontsize=fslabel)
plt.xlim(abin[0],abin[-1]); #plt.ylim(1.e-10,1e-6)
plt.grid(True)
plt.legend(fontsize=fslegend,loc='best')
plt.savefig(figpre+'wholeLF_Lbol_z'+str(int(z))+'eta'+str(int(10*eta))+'f'+str(f_duty)+'mu'+str(mu_fit)+'sigma'+str(sigma_fit)+'N'+str(int(np.log10(N_concatenate)))+'.png')