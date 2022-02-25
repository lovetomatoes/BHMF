from PYmodule import *
# same as LFmcz6, using only 1 tree file (Mh=1e11, ibsm=0 fbsm=1); tighter prior constraints
  
# pre-determined in PYmodule
# eta = 0.3
# alpha = 1.

t_range = [120.]
f_range = [1.]
d_range = [.25]
l_range = [.9]
a_range = [.1]
# Chi2_M = 1.7; Chi2_L = 0.68 quite small.
# off_M = 1.6; off_L = 2.2

# LF bins same w/ Matsu18
z = int(6)
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)


T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

def lnlike(theta):
    t_life, f_0, d_fit, l_cut, a = theta
    t_life *= Myr
## --------- Mass Function ---------
    tz = t_from_z(z)
    dn_MBH = np.zeros(N_mf)
    Nt = np.max((tz-T['t_col'])//t_life)
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
                # #----------- Schechter lbd -----------
                x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)
                dP_seed = np.nansum(dP_seed)/len(T)
            else:
                dP_seed = 0.
            # prev BHMF
            # #----------- Schechter lbd -----------
            x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, f_0, l_cut, d_fit)
            x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, f_0, l_cut, d_fit)
            x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
            dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
        Nt -= 1
    dn_MBH = dP_MBH*n_base*f_bsm

    consv_ratio = np.nansum(dn_MBH)/n_base
    # print('MF consv_ratio',consv_ratio)
    # wli: too loose constraint!
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta)
        print('consv_ratio: ',consv_ratio)
        # assert 0


# # --------- Luminosity Function ---------
    Phi = np.zeros(N_lf)
    for ibin in range(N_lf):
        #----------- Schechter lbd -----------
        x0 = kernelS_M1450(bin_edg[ibin+1], M_BH, l_cut)
        x1 = kernelS_M1450(bin_edg[ibin],   M_BH, l_cut)
        dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

        dPhi = np.nansum(dn_MBH*dP_M1450)
        Phi[ibin] = dPhi/bin_wid[ibin]

    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(bin_cen)
    # Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
    # print('Chi2=%.2e'%Chi2)
    ys = np.log(Phi_obs)
    y_model = np.log(Phi_DO)
    y_err = np.log(Phi_err)
    Chi2 = np.sum( pow((ys - y_model)/y_err, 2))
    if not np.isfinite(Chi2):
        print('theta=',theta)
        print('inf or nan? Chi2=',Chi2)
        return -np.inf
    return -.5*Chi2

# # 7even* ... & posPara
# def lnprior(theta):
#     t_life, f_0, d_fit, l_cut, a = theta
#     if t_life>0. and f_0>0. and d_fit>0. and l_cut>0. and a>0.:
#         return 0.0
#     else:
#         return -np.inf

# # 4* ...
# def lnprior(theta):
#     t_life, f_0, d_fit, l_cut, a = theta
#     if 1e1<t_life<1e3 and d_fit>0. and f_0>0. and 0.1<f_0*l_cut<10. and a>0.:
#         return 0.0
#     else:
#         return -np.inf

# range: 1e1<t_life<1e3 and .1<f_0<10. and d_fit>0. and .1<l_cut<10. and a>0.:
# range1: 1e1<t_life<1e3 and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range2: 1e1<t_life and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range3: 1e1<t_life<200 and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range4: 1e1<t_life<200 and .1<f_0<2. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range5: 1e1<t_life<200 and .1<f_0<2. and 0<d_fit<.5 and .1<l_cut<10. and a>0.:

def lnprior(theta):
    t_life, f_0, d_fit, l_cut, a = theta
    if 1e1<t_life<200. and .1<f_0<2. and 0<d_fit<.5 and .1<l_cut<10. and a>0.:
        return 0.0
    else:
        return -np.inf

def lnprobab(theta):
    lp = lnprior(theta)
    # print('theta',theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)