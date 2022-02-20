from PYmodule import *
# Phiz6 paras: f_0, t_life, lbd ~Schechter(l_cut,a)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.

# print(M1M0(np.array([1e2,1e5]),tz-t_from_z(20),1,.1,.1,0.001))
# print(M1M0(np.array([1e2,1e5]),50*Myr,1,.1,.1,0.001))

# f_duty=.7;mu_fit=.21*30;sigma_fit=.15
# M0 = 1e4; dt = 500*Myr
# print('M1=%.1e'%(M0*np.exp(mu_fit*f_duty*dt/t_Edd)), 'N_mf=%d'%N_mf)

# MF
abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
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
N_lf = len(bin_cen)

t_range = np.array([100, 200, 500, 1000])*Myr
f_range = np.arange(.1, 1., .3)
l_range = np.logspace(-2, 1., num=4) # l_cut
a_range = np.array([.1, 1., 2., 3.]) # a>0 total P convergent
# print(len(t_range)*len(f_range)*len(l_range)*len(a_range) )
# print(a_range); exit(0)

t_range = np.arange(100,1000,100)*Myr
f_range = np.arange(.1, 1., .1)
l_range = np.append( [.01,.05], np.arange(.1,2.,.1))
a_range = np.arange(.1, 3., .1) # a>0 total P convergent
# print(len(t_range)*len(f_range)*len(l_range)*len(a_range) )
# exit(0)

t_range = [120.*Myr]
f_range = [1.]
d_range = [.25]
l_range = [.9]
a_range = [.1]

t_range = [6.93021461e+02*Myr]
f_range = [2.54469738e+00]
d_range = [9.42247174e-01]
l_range = [1.78416875e+00]
a_range = [6.13968336e-03]



# Chi2_M = 1.7; Chi2_L = 0.68 quite small.
# off_M = 1.6; off_L = 2.2

T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

i = 0
Chi2_min = 1e10; find_min = False
for t_life in t_range:
    for f_0 in f_range:
        for d_fit in d_range:
            for l_cut in l_range:
                for a in a_range:
                    i = i+1
                    # continue
                ## --------- Mass Function ---------
                    dn_MBH = np.zeros(N_mf)
                    Nt = np.max((tz-T['t_col'])//t_life)
                    Nmax = Nt
                    dP_MBH = np.zeros(N_mf)
                    # t_tot = np.zeros(len(T))
                    while Nt>=0:
                        t_point = tz - Nt*t_life
                        T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
                        dt_seed = t_point - T_seed['t_col']
                        dP_MBH_prev = dP_MBH.copy()
                        for ibin in range(N_mf):
                            # new seeds
                            if len(T_seed):
                                x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                                x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                                x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                for i in range(len(x0)):
                                    assert x0[i]<= x1[i]
                                dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)

                                dP_seed = np.nansum(dP_seed)/len(T)
                            else:
                                dP_seed = 0.
                            x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, f_0, l_cut, d_fit)
                            x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, f_0, l_cut, d_fit)
                            x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                            for i in range(len(x0)):
                                assert x0[i]<= x1[i]
                            dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
                        Nt -= 1
                    dn_MBH = dP_MBH*n_base*f_bsm

                    consv_ratio = np.nansum(dn_MBH)/n_base
                    print('MF conserved fraction=%.10f'%consv_ratio)
                    # if consv_ratio<.9:
                    #     print('conserved fraction=%.10f'%consv_ratio)

                    T = Table(
                        [M_BH, dn_MBH/dlog10M, MF(M_BH)],
                        names=('M_BH','Phi','W10_MF')
                    )
                    MFname = z6datapre+'MF_SC'+'t%.1e'%(t_life/Myr)+ \
                            'f%.1f'%f_0+ \
                            'd%.1f'%d_fit+ \
                            'l%.1e'%l_cut+ \
                            'a%.3f'%a+ \
                            'alpha%.1f'%alpha
                    MFname = z6datapre+'MFIMF'
                    ascii.write( Table([np.log10(T['M_BH']), T['Phi'], T['W10_MF']],
                                names=['M_BH','Phi','W10_MF']),
                                MFname,
                                formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e'},
                                overwrite=True)
                    # exit(0)
                    # T  = T[np.logical_and(True,T['M_BH']<2e10)] # select M_BH range
                    index = np.logical_and(T['M_BH']>1e7, T['M_BH']<1e10)
                    Chi2_M = np.sum(pow(np.log(T['Phi']/T['W10_MF'])[index], 2))/(np.sum(index)-1)
                    off_M = np.max(abs(np.log(T['Phi']/T['W10_MF'])[index]))

                # # --------- Luminosity Function ---------
                    Phi = np.zeros(N_lf)
                    Phi_csv = 0.
                    T['dn_MBH'] = T['Phi']*dlog10M
                    for ibin in range(N_lf):
                        #----------- Schechter lbd -----------
                        x0 = kernelS_M1450(bin_edg[ibin+1], T['M_BH'], l_cut)
                        x1 = kernelS_M1450(bin_edg[ibin],   T['M_BH'], l_cut)
                        dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

                        dPhi = np.nansum(T['dn_MBH']*dP_M1450)
                        Phi_csv += dPhi
                        Phi[ibin] = dPhi/bin_wid[ibin]
                    # print('consv of dP_M1450:',Phi_csv/np.sum(n_base))

                    Phi *= 1e9
                    Phi_DO = Phi/corr_U14D20(bin_cen)
                    Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
                    off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

                    T = Table(
                        [bin_cen,Phi_obs,Phi_DO,Phi,Chi2*np.ones(N_lf)],
                        names=('bin_cen','Phi_obs','Phi_DO','Phi','Chi2')
                    )
                    LFname = z6datapre+'LF_SC'+'t%.1e'%(t_life/Myr)+ \
                            'f%.1f'%f_0+ \
                            'd%.1f'%d_fit+ \
                            'l%.1e'%l_cut+ \
                            'a%.3f'%a+ \
                            'alpha%.1f'%alpha
                    LFname = z6datapre+'LFIMF'
                    ascii.write(T, LFname,
                                formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
                                overwrite=True)

                    if np.nanmin([Chi2, Chi2_min]) == Chi2:
                        find_min = True
                        Chi2_min = Chi2
                        T_min = T
                        LFname_min = LFname
                    # exit(0)
# print(i)
if find_min:
    print('log_prob=',-.5*Chi2_min*(len(Phi_obs)-1))
    #, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)