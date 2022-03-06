from PYmodule import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.
# N_Mh = 1 now

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

f_0 = 1.
d_fit = 0.
l_cut = .7 # l_cut=2., l_cut' = l_cut/2; M=M0=1e7 grow as Eddington
a = .3

def model(theta, z = int(6)):
    # t_life, f_0, d_fit, l_cut, a = theta
    t_life = theta
    t_life = t_life * Myr # wli: not *=!!!!!! theta changed by this
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
        print('theta: ',theta, 'consv_ratio: ',consv_ratio)
        # assert 0

    # 30 N_M in 1e7-1e10 range, plus 12 N_L
    index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    xs = M_BH[index]
    ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
    y_model = np.log10( (dn_MBH/dlog10M) [index] )
    y_err = 1.
    Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

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
    ys = np.log(Phi_obs)
    y_model = np.log(Phi_DO)
    y_err = np.log(Phi_err)
    Chi2_L = np.sum( pow((ys - y_model)/y_err, 2))
    # print('in model -.5*(Chi2_L+Chi2_M)',-.5*(Chi2_L+Chi2_M))
    return {'M_BH':M_BH, 'MF':dn_MBH/dlog10M, 'MF_data':MF(M_BH), 'Chi2_M':Chi2_M,
            'M1450':bin_cen, 'LF':Phi_DO, 'LF_data':Phi_obs, 'Chi2_L':Chi2_L}