from PYmodule import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.

logMs = np.linspace(7,10,num=4)

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
    t_life, d_fit, logM0 = theta
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

        # new seeds (using 2d meshgrids)
        if len(T_seed):
            z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit, logM0)
            z_mesh[z_mesh<0] = 0.
            dP_seed = special.gammainc(a,z_mesh[1:,:]) - special.gammainc(a,z_mesh[:-1,:])
            dP_seed = np.nansum(dP_seed, axis=1)/len(T)
        else:
            dP_seed = 0.
        # prev BHMF
        z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit, logM0)
        z_mesh[z_mesh<0] = 0.
        dP_MBH = np.nansum((special.gammainc(a,z_mesh[:,:-1])-special.gammainc(a,z_mesh[:,1:]))*dP_MBH_prev, axis=1) + dP_seed

        Nt -= 1
    dn_MBH = dP_MBH*n_base*f_bsm

    consv_ratio = np.nansum(dn_MBH)/n_base
    # print('MF consv_ratio',consv_ratio)
    # wli: too loose constraint!
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta)
        print('consv_ratio: ',consv_ratio)
        return -np.inf

    # 30 N_M in 1e7-1e10 range, plus 12 N_L
    index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    xs = M_BH[index]
    ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
    y_model = np.log10( (dn_MBH/dlog10M) [index] )
    y_err = 1.
    Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
    Phi = np.zeros(N_lf)
    
    z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
    P_mesh = special.gammainc(a,z_mesh[:-1,:])-special.gammainc(a,z_mesh[1:,:])
    dPhi_mesh = np.nansum(P_mesh*dn_MBH,axis=1)
    Phi = dPhi_mesh/bin_wid

    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(bin_cen)
    # Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
    # print('Chi2=%.2e'%Chi2)
    ys = np.log(Phi_obs)
    y_model = np.log(Phi_DO)
    y_err = np.log(Phi_err)
    Chi2_L = np.sum( pow((ys - y_model)/y_err, 2))
    if not np.isfinite(Chi2_L):
        print('theta=',theta)
        print('inf or nan? Chi2_L=',Chi2_L)
        return -np.inf
    return -.5*(Chi2_M + Chi2_L)


# range: 1e1<t_life<1e3 and .1<f_0<10. and d_fit>0. and .1<l_cut<10. and a>0.:
# range1: 1e1<t_life<1e3 and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range2: 1e1<t_life and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range3: 1e1<t_life<200 and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range4: 1e1<t_life<200 and .1<f_0<2. and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range5: 1e1<t_life<200 and .1<f_0<2. and 0<d_fit<.5 and .1<l_cut<10. and a>0.:
# range6: 1e1<t_life<200 and .1<l_cut<10. and 0.<a<1.:
# range7: 1e1<t_life<200 and .1<l_cut<10. and 0.1<a<0.5:
# 2prange1: 1e1<t_life<200. and 0.1<d_fit<0.5:
# 4prange1: 1e1<t_life<200. and 0.1<d_fit<0.5:
# 3prange1: 1e1<t_life<200. and 0.1<d_fit<0.5 and 5.<logM0<8.:

def lnprior(theta):
    t_life, d_fit, logM0 = theta
    if 1e1<t_life<200. and 0.1<d_fit<0.5 and 5.<logM0<8.:
        return 0.0
    else:
        return -np.inf

def lnprobab(theta):
    lp = lnprior(theta)
    # print('theta',theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)