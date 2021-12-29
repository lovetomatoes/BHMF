from PYmodule import *
# Phiz6 paras: t_life, f_duty, lognorm(mu,sigma) or Schechter(l_cut,a)

# print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
# print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
# print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
# print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
# print('z:17.58 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(17.58))/Myr)
# print('t_Edd         : %3.2f Myr', t_Edd/Myr)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.


Ts = [] # bsm=0,1 two files
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        T['Lbol_z'] = np.zeros(len(T))
        T['M1450_z'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['t_col'] = t_from_z(T['z_col'])
        TM.append(T)
        # print(np.min(T['z_col']),np.max(T['z_col']),np.min(T['Mstar0']),np.max(T['Mstar0']))
    Ts.append(TM)

# print(M1M0(np.array([1e2,1e5]),tz-t_from_z(20),1,.1,.1,0.001))
# print(M1M0(np.array([1e2,1e5]),50*Myr,1,.1,.1,0.001))

# MF
abin_mf =  np.logspace(2,12,num=400) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
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


d_range = np.array([0., .001, .01, .1])
e_range = np.array([0., .05, .1, .15, .2])
f_range = np.arange(.5, 1., .1)
m_range = np.arange(.2, .5, .01)
s_range = np.arange(.01, 0.2, .01)
t_range = np.array([100, 200, 500, 1000])*Myr

d_range = [.0]
e_range = [.0]
f_range = [.7]
m_range = [.21]
s_range = [.15]
a_range = [1] # a>0 total P convergent
l_range = [.1] # l_cut

# for iM in range(N_Mh):
#     for i_bsm in range(Nbsm):
#         T = Ts[iM][i_bsm]
#         dt = tz-t_from_z(T['z_col'])
#         print('Mh=',Mhpres[iM][-2:],'ibsm=%d'%i_bsm,'Nmax=%.2f'%np.max(dt/t_life),'Nmin=%.2f'%np.min(dt/t_life))

# f_duty=.7;mu_fit=.21*30;sigma_fit=.15
# M0 = 1e4; dt = 500*Myr
# print('M1=%.1e'%(M0*np.exp(mu_fit*f_duty*dt/t_Edd)), 'N_mf=%d'%N_mf)
# x0 = kernel_MBH1(bin_left /M0, dt, f_duty, mu_fit, sigma_fit)
# x1 = kernel_MBH1(bin_right/M0, dt, f_duty, mu_fit, sigma_fit)
# dP_seed = -.5*(special.erfc(x1) - special.erfc(x0))
# # print(x0, x1)
# print('log-normal: P_tot=',np.nansum(dP_seed))

# x0 = kernelS_MBH(bin_left /M0, dt, f_duty, l_cut)
# x1 = kernelS_MBH(bin_right/M0, dt, f_duty, l_cut)
# x0[x0<0] = 0.; x1[x1<0] = 0.
# dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)
# # print(x0, x1)
# x0 = ma.masked_less(x0,0.); x1 = ma.masked_less(x1,0.) # 必删! 否则not conserved
# print('Schechter: P_tot=',np.nansum(dP_seed))
# print('Schechter: P_tot_ana=',special.gammainc(a,np.inf))
# exit(0)

i = 0
t_life = 100.*Myr
l_cut = .21
a = 1.

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
                    for iM in range(N_Mh): #[0]:#
                        for i_bsm in range(Nbsm): #[0]:#
                            T = Ts[iM][i_bsm]
                            Nt = np.max((tz-T['t_col'])//t_life)
                            Nmax = Nt
                            dP_MBH = np.zeros(N_mf)
                            while Nt>=0:
                                t_point = tz - Nt*t_life
                                T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
                                dt_seed = t_point - T_seed['t_col']
                                dP_MBH_prev = dP_MBH.copy()
                                # M_BHt = M1M0(M_BH,dt,f_duty,mu_fit,eta8,delta_fit) # if not exp growth, may needed
                                for ibin in range(N_mf): #[100]:#
                                    # new seeds 
                                    if len(T_seed):
                                        # #----------- log-norm lbd -----------
                                        x0 = kernel_MBH1(bin_left[ibin] /T_seed['Mstar0'],dt_seed,f_duty, mu_fit, sigma_fit)
                                        x1 = kernel_MBH1(bin_right[ibin]/T_seed['Mstar0'],dt_seed,f_duty, mu_fit, sigma_fit)
                                        # x0 = ma.masked_invalid(x0); x1 = ma.masked_invalid(x1) # 必删! or not conserved
                                        dP_seed = -.5*(special.erfc(x1) - special.erfc(x0)) # P(<x) = 1-.5*erfc(x)
                                        #----------- Schechter lbd -----------
                                        x0 = kernelS_MBH(bin_left[ibin] /T_seed['Mstar0'], dt_seed, f_duty, l_cut)
                                        x1 = kernelS_MBH(bin_right[ibin]/T_seed['Mstar0'], dt_seed, f_duty, l_cut)
                                        x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                        dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)

                                        dP_seed = np.nansum(dP_seed)/len(T) 
                                    else:
                                        dP_seed = 0.
                                    # prev BHMF
                                    # #----------- log-norm lbd -----------
                                    x0 = kernel_MBH1(M_BH[ibin]/bin_right,t_life,f_duty, mu_fit, sigma_fit)
                                    x1 = kernel_MBH1(M_BH[ibin]/bin_left, t_life,f_duty, mu_fit, sigma_fit)
                                    dP_MBH[ibin] = np.nansum(.5*(special.erfc(x0) - special.erfc(x1)) * dP_MBH_prev) + dP_seed
                                    #----------- Schechter lbd -----------
                                    x0 = kernelS_MBH(M_BH[ibin]/bin_right, t_life, f_duty, l_cut)
                                    x1 = kernelS_MBH(M_BH[ibin]/bin_left,  t_life, f_duty, l_cut)
                                    x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                    dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
                                Nt -= 1
                            dn_MBH += dP_MBH*n_base[iM]*f_bsm[i_bsm]

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
                                   'alpha%.1f'%alpha+
                                   'N%d'%Nmax,
                                formats={'M_BH':'4.2e','dn_dlog10M':'4.2e'},
                                overwrite=True)
                    consv_ratio = np.nansum(dn_MBH)/np.sum(n_base)
                    print('conserved fraction=%.10f'%consv_ratio)
                    if consv_ratio<.9:
                        print('conserved fraction=%.10f'%consv_ratio)
                    exit(0)
                    T  = T[np.logical_and(True,T['M_BH']<2e10)] # select M_BH range

                    Phi = np.zeros(N_lf)
                    for ibin in range(N_lf): # N_lf
                        M_BH = T['M_BH']
                        # ----------- log-norm lbd -----------
                        kernel = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
                        kernel = ma.masked_outside(kernel, -10., 10.)
                        x0 = kernel_M1450(bin_edg[ibin+1], M_BH, mu_fit, sigma_fit)
                        x1 = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
                        dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))
                        # ----------- Schechter lbd -----------

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
                                'alpha%.1f'%alpha+
                                'N%d'%Nmax,
                                formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_CO':'4.2e','Phi_DO':'4.2e','Phi_obs':'4.2e'},
                                overwrite=True)
                    # exit(0)
print(i)