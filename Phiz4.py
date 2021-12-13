from PYmodule import *

# print('z:45, t_Hubble: %3.2f Myr'%(t_from_z(45)/Myr))
# print('z:30, t_Hubble: %3.2f Myr'%(t_from_z(30)/Myr))
# print('z:10, t_Hubble: %3.2f Myr'%(t_from_z(10)/Myr))
# print('z:5,  t_Hubble: %3.2f Myr'%(t_from_z(5)/Myr))
# print('z:32 to 6     : %3.2f Myr'%((t_from_z(6)-t_from_z(32))/Myr))
# print('t_Edd         : %3.2f Myr'%(t_Edd/Myr))

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(4)

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
    Ts.append(TM)


# LF bins same w/ Aki18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)

# for M1450 in bin_cen:
#     print('M1450',M1450,'Eddington accretion Mbh = %3.2e'%M_L(Lbol_M1450(M1450),.1))
# print('M1450 of 1e6Msun: 10Edd',M1450_Lbol(L_M(1e6,10)))

# MF
abin_mf =  np.logspace(2,13,num=200) # default endpoint=True
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
# M_BH: mass bin center; dn_MBH: # density in mass bins
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
N_mf = len(M_BH)
dn_MBH = np.zeros(N_mf)

# ## --------------------------   z=5   ----------------------------
# # best fit from z=5 LF (M1450_min = ?; M1450_max = ?)
# # currently first fitting best
f_duty = .6;  mu_fit = .1; sigma_fit = .4
eta8 = .05; delta_fit = .4

f_best_5 = z5datapre+'MF_best_5'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
T = ascii.read(f_best_5, guess=False, delimiter=' ') 
# print(T)

d_range = [.05, .1, .15, .2, .25, .3, .35, .4, .45, .5]
e_range = [.05, .08, .1, .12, .15]
f_range = np.arange(.1, 1.1, .1)
m_range = np.logspace(-2, 0, num=5)
s_range = np.arange(.1, 0.5, .1)

d_range = [.1]
e_range = [.05]
f_range = [1.]
m_range =[.1]
s_range = [.1]

d_range = [.3]
e_range = [.1]
f_range = [.4]
m_range = [.4]
s_range = [.2]

d_range = [.35]
e_range = [.05]
f_range = [.1]
m_range = [.8]
s_range = [.01]

i = 0
## --------------------------   z=z   ----------------------------
for delta_fit in d_range: # .25 in np.arange(.1, .4, .05)
    for eta8 in e_range: # .07 in np.arange(.05, .13, .01)
        if (eta8 == 0. and delta_fit == 0.):
            pass
        elif eta8*delta_fit == 0.:
            continue

        for f_duty in f_range: # .6 in np.arange(.2, 1., .1)
            for mu_fit in m_range: # .1 in np.logspace(-2, 0, num=5)
                for sigma_fit in s_range: # .29 in np.arange(.2, 0.3, .01)
                    i = i+1
                    # continue
                    fname = z4datapre+'LF_'+'z%d'%z+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    # # check if file exists
                    # if os.path.isfile(fname):
                    #     continue

                    # converge check timestep
                    dt = t_from_z(z)-t_from_z(5.)
                    # M_BHz: mass bin center at z
                    M_BHz = M1M0(M_BH,dt,f_duty,mu_fit,eta8,delta_fit)
                    dn_MBH = np.zeros(N_mf)
                    for ibin in range(N_mf):
                        x0 = kernel_MBH2(M_BHz[ibin],T['bin_right'],dt,f_duty,mu_fit,sigma_fit,eta8,delta_fit)
                        x1 = kernel_MBH2(M_BHz[ibin],T['bin_left'],dt,f_duty,mu_fit,sigma_fit,eta8,delta_fit)
                        # x0 = ma.masked_invalid(x0); x1 = ma.masked_invalid(x1) # no use
                        dP_MBH = .5*(special.erfc(x0) - special.erfc(x1)) * T['dn_MBH']
                        dn_MBH[ibin] = np.nansum(dP_MBH)
                    # print('i=%d'%i,'z=%.1f:'%z,'percent of bins',np.nansum(dn_MBH)/N_bin_ini)
                    # print('min M_BHz @ z=%.1e:'%np.nanmin(M_BHz),'max M_BHz @ z=%.1e:'%np.nanmax(M_BHz))

                    Tz = Table(
                        [M_BHz, dn_MBH],
                        names=('M_BH','dn_MBH')
                    )
                    if (eta8 == 0. and delta_fit == 0.):
                        dlog10Mz = dlog10M
                    else:
                        dlog10Mz = np.log10(M1M0(abin_mf[1:],dt,f_duty,mu_fit,eta8,delta_fit)/M1M0(abin_mf[:-1],dt,f_duty,mu_fit,eta8,delta_fit))
                    #!!!!!!!!!!!!! dlog10Mz[0] instead of dlog10Mz to make MF@z4 smooth...
                    ascii.write(Table([Tz['M_BH'], Tz['dn_MBH']/dlog10Mz[0]], names=['M_BH','dn_dlog10M']),
                                   z4datapre+
                                   'MF_'+
                                   # 'piece'
                                   'z%d'%z+
                                   'f%3.2f'%f_duty+
                                   'm%3.2f'%mu_fit+
                                   's%3.2f'%sigma_fit+
                                   'e%.3f'%eta8+
                                   'd%.3f'%delta_fit+
                                   'alpha%.1f'%alpha,
                                   formats={'M_BH':'4.2e','dn_dlog10M':'4.2e'},
                                   overwrite=True)                    
                    # exit(0)

                    Phi = np.zeros(N_lf)
                    for ibin in range(N_lf): # N_lf
                        kernel = kernel_M1450(bin_edg[ibin], M_BHz, mu_fit, sigma_fit)
                        # kernel = ma.masked_outside(kernel, -10., 10.)
                        x0 = kernel_M1450(bin_edg[ibin+1], M_BHz, mu_fit, sigma_fit)
                        x1 = kernel_M1450(bin_edg[ibin], M_BHz, mu_fit, sigma_fit)
                        dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))
                        dPhi = np.nansum(Tz['dn_MBH']*dP_M1450)
                        Phi[ibin] += dPhi/bin_wid[ibin]*f_duty

                    Tlf = Table(
                        [bin_cen, Phi*1e9, Phi*1e9*(1.-f_obsc_const), Phi*1e9/corr_U14D20(bin_cen), Phi_obs],
                        names=('bin_cen','Phi','Phi_CO','Phi_DO','Phi_obs')
                    )
                    ascii.write(Tlf, z4datapre+
                                    'LF_'+
                                    'z%d'%z+
                                    'f%3.2f'%f_duty+
                                    'm%3.2f'%mu_fit+
                                    's%3.2f'%sigma_fit+
                                    'e%.3f'%eta8+
                                    'd%.3f'%delta_fit+
                                    'alpha%.1f'%alpha,
                                    formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_CO':'4.2e','Phi_DO':'4.2e','Phi_obs':'4.2e'},
                                    overwrite=True)
                    # if i==2:
                    #     exit(0)
print(i)