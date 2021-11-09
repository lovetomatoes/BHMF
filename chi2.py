from PYmodule import *
# intrinsic LF @z=6; similar w/ chi2z6, using best para as input; plot w/ finer bins

Chi2_min = {'no_corr':1e10, 'CO':1e10, 'DO':1e10}
Chi2_min = 1e10
z = int(6)

bin_edg = np.arange(-29,-22,0.1)
bin_wid = bin_edg[1:] - bin_edg[:-1]
bin_cen = bin_edg[:-1] + bin_wid/2.
N_lf = len(bin_cen)
# print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)
M1450_min = -30.; M1450_max = -22.

alpha = 1.

f_duty = .5; mu_fit = .38; sigma_fit = .12

# MF
N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
Ts = [] # bsm=0,1 two files
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        T['Mstar_z'] = np.zeros(len(T))
        T['Lbol_z'] = np.zeros(len(T))
        T['M1450_z'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['type'] = np.ones(len(T))
        T['Edd_ratio'] = np.ones(len(T))
        TM.append(T)
    Ts.append(TM)
abin_mf =  np.logspace(6,12,num=100) # default endpoint=True
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1
dn_MBH = np.zeros(N_mf)
for ibin in range(N_mf): # N_mf
    for iM in range(N_Mh):
        for i_bsm in range(Nbsm):
            T = Ts[iM][i_bsm]
            dP = 0
            x0 = kernel_MBH(abin_mf[ibin]/T['Mstar0'],t_from_z(z)-t_from_z(T['z_col']),f_duty, mu_fit, sigma_fit)
            x1 = kernel_MBH(abin_mf[ibin+1]/T['Mstar0'],t_from_z(z)-t_from_z(T['z_col']),f_duty, mu_fit, sigma_fit)
            x0 = ma.masked_invalid(x0); x1 = ma.masked_invalid(x1) 
            dP_MBH = .5*(special.erfc(x0) - special.erfc(x1))
            dP = np.sum(dP_MBH)/Ntr
            dn_MBH[ibin] += dP*n_base[iM]*f_bsm[i_bsm]
T = Table(
    [abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0]), dn_MBH],
    names=('M_BH','dn_MBH')
)

T  = T[np.logical_and(T['M_BH']>1e6,T['M_BH']<2e10)] # select M_BH range

Phi = np.zeros(N_lf)
for ibin in range(N_lf): # N_lf
    M_BH = T['M_BH']
    kernel = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
    kernel = ma.masked_outside(kernel, -10., 10.)
    x0 = kernel_M1450(bin_edg[ibin+1], M_BH, mu_fit, sigma_fit)
    x1 = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
    dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))
    dPhi = np.sum(T['dn_MBH']*dP_M1450)
    Phi[ibin] += dPhi/bin_wid[ibin]*f_duty

T = Table(
    [bin_cen, Phi*1e9, Phi*1e9*(1.-f_obsc_const), Phi*1e9/corr_U14D20(bin_cen)],
    names=('bin_cen','Phi','Phi_CO','Phi_DO')
)

plt.figure(figsize=(10,10),dpi=200)
# plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
plt.plot(bin_cen, LF_M1450(bin_cen,z)*1e9, label=lfnames[str(z)], c='C0')
plt.scatter(bin_cen, T['Phi_DO'], c='C1', label='histogram')
plt.text(-26,1e3,'f=%3.2f'%f_duty+'\n'+r'$\mu=$'+'%3.2f'%mu_fit+'\n'+r'$\sigma=$'+'%3.2f'%sigma_fit,fontsize = fstxt)
plt.xlim(bin_edg[-1],bin_edg[0]); plt.ylim(5e-3,1e4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.savefig(z6figpre+'_chi2_Phi_DO_'+str(z)+'MBH2e10'+
                'M'%+
                'f%3.2f'%f_duty+
                'm%3.2f'%mu_fit+
                's%3.2f'%sigma_fit+
                'alpha%.1f'%alpha+'.png')