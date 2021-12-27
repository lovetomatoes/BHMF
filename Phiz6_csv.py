from PYmodule import *

# print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
# print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
# print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
# print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
print('z:20 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(32))/Myr)
# print('t_Edd         : %3.2f Myr', t_Edd/Myr)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
alpha = 1.


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
        # print(np.min(T['z_col']),np.min(T['Mstar0']),np.max(T['Mstar0']))
    Ts.append(TM)

print(M1M0(np.array([1e2,1e5]),t_from_z(z)-t_from_z(20),1,.1,.1,0.001))
# MF
abin_mf =  np.logspace(2,15,num=200) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]

# bin_edg = np.linspace(-29.,-22.)
# bin_wid = bin_edg[1:] - bin_edg[:-1]
# bin_cen = (bin_edg[1:] + bin_edg[:-1])/2.

N_lf = len(bin_cen)

# for M1450 in bin_cen:
#     print('M1450',M1450,'Eddington accretion Mbh = %3.2e'%M_L(Lbol_M1450(M1450),.1))


d_range = [0., .001, .01, .1]
e_range = [0., .05, .1, .15, .2]
f_range = np.arange(.5, 1., .1)
m_range = np.arange(.2, .5, .01)
s_range = np.arange(.01, 0.2, .01)

d_range = [.0]
e_range = [.0]
f_range = [.7]
m_range = [.21]
s_range = [.15]

i = 0
for delta_fit in d_range: # .001 in [0., .001, .01, .1]
    for eta8 in e_range:# .1 in [0., .05, .1, .15, .2]:
        if (eta8 == 0. and delta_fit == 0.):
            pass
        elif eta8*delta_fit == 0.:
            continue
        for f_duty in f_range: # .7 in np.arange(.5, 1., .1):
            for mu_fit in m_range: # .21 in np.arange(.2, .5, .01):
                for sigma_fit in s_range: # .15 in np.arange(.01, 0.2, .01):
                    i = i+1
                    # continue
                    dn_MBH = np.zeros(N_mf)
                    
                    for iM in range(N_Mh):
                        for i_bsm in range(Nbsm):
                            T = Ts[iM][i_bsm]
                            for ibin in range(N_mf): # N_mf
                                dt = t_from_z(z)-t_from_z(T['z_col'])
                                if (eta8 == 0. and delta_fit == 0.): # eta=0.1 constant; e,d not into fitting
                                    x0 = kernel_MBH1(abin_mf[ibin]/T['Mstar0'],dt,f_duty, mu_fit, sigma_fit)
                                    x1 = kernel_MBH1(abin_mf[ibin+1]/T['Mstar0'],dt,f_duty, mu_fit, sigma_fit)
                                else :
                                    x0 = kernel_MBH2(abin_mf[ibin],T['Mstar0'],dt,f_duty,mu_fit,sigma_fit,eta8,delta_fit)
                                    x1 = kernel_MBH2(abin_mf[ibin+1],T['Mstar0'],dt,f_duty,mu_fit,sigma_fit,eta8,delta_fit)
                                # x0 = ma.masked_invalid(x0); x1 = ma.masked_invalid(x1) # 必删! 否则not conserved
                                dP_MBH = .5*(special.erfc(x0) - special.erfc(x1))
                                # print(np.sum(~np.isnan(dP_MBH)))
                                dP = np.nansum(dP_MBH)/Ntr
                                dn_MBH[ibin] += dP*n_base[iM]*f_bsm[i_bsm]
                    T = Table(
                        [M_BH, dn_MBH],
                        names=('M_BH','dn_MBH')
                    )

                    ascii.write( Table([T['M_BH'], T['dn_MBH']/dlog10M], names=['M_BH','dn_dlog10M']),
                                   z6datapre+
                                   'MF'+
                                   'f%3.2f'%f_duty+
                                   'm%3.2f'%mu_fit+
                                   's%3.2f'%sigma_fit+
                                   'e%.3f'%eta8+
                                   'd%.3f'%delta_fit+
                                   'alpha%.1f'%alpha,
                                formats={'M_BH':'4.2e','dn_dlog10M':'4.2e'},
                                overwrite=True)
                    # print('np.sum(dn_MBH)',np.nansum(dn_MBH),'np.sum(n_base)',np.sum(n_base))
                    print('conserved fraction',np.nansum(dn_MBH)/np.sum(n_base))
                    exit(0)
                    T  = T[np.logical_and(True,T['M_BH']<2e10)] # select M_BH range

                    Phi = np.zeros(N_lf)
                    for ibin in range(N_lf): # N_lf
                        M_BH = T['M_BH']
                        kernel = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
                        kernel = ma.masked_outside(kernel, -10., 10.)
                        x0 = kernel_M1450(bin_edg[ibin+1], M_BH, mu_fit, sigma_fit)
                        x1 = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
                        dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))
                        dPhi = np.nansum(T['dn_MBH']*dP_M1450)
                        Phi[ibin] += dPhi/bin_wid[ibin]*f_duty

                    T = Table(
                        [bin_cen, Phi*1e9, Phi*1e9*(1.-f_obsc_const), Phi*1e9/corr_U14D20(bin_cen), Phi_obs],
                        names=('bin_cen','Phi','Phi_CO','Phi_DO','Phi_obs')
                    )
                    ascii.write(T, z6datapre+
                                'LF2e10_'+
                                'f%3.2f'%f_duty+
                                'm%3.2f'%mu_fit+
                                's%3.2f'%sigma_fit+
                                'e%.3f'%eta8+
                                'd%.3f'%delta_fit+
                                'alpha%.1f'%alpha,
                                formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_CO':'4.2e','Phi_DO':'4.2e','Phi_obs':'4.2e'},
                                overwrite=True)
                    # exit(0)
print(i)