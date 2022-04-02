from PYmodule import *
from PYmodule.l_intg import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.
# N_Mh = 1 now

z = int(6)
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]


def model(theta, z = int(6), f_0=f_0, d_fit=d_fit, logM0=logM0, l_cut= l_cut, a=a):
    # t_life, f_0, d_fit, l_cut, a = theta
    if len(theta) == 1:
        t_life = theta
    elif len(theta) == 2:
        t_life, d_fit = theta
    elif len(theta) == 3:
        # t_life, l_cut, a = theta
        t_life, d_fit, logM0 = theta
    elif len(theta) == 4:
        t_life, d_fit, l_cut, a = theta
    else:
        assert 0
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

        # new seeds (using 2d meshgrids)
        if len(T_seed):
            z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
            z_mesh[z_mesh<x0] = x0
            Ps = integral(a,z_mesh)/I_toinf
            dP_seed = Ps[1:,:] - Ps[:-1,:]
            dP_seed = np.nansum(dP_seed, axis=1)/len(T)
        else:
            dP_seed = 0.
        # prev BHMF
        z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh)/I_toinf
        dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed

        Nt -= 1

    dn_MBH = dP_MBH*n_base*f_bsm

    consv_ratio = np.nansum(dn_MBH)/n_base
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta, 'consv_ratio: ',consv_ratio)
        # assert 0

    # 30 N_M in 1e7-1e10 range, plus 12 N_L
    index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    xs = M_BH[index]
    ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
    y_model = np.log10( (dn_MBH/dlog10M) [index] )
    y_err = 1.
    y_err = (np.log10(xs)-8.5)**2/3 + .5 # from 0.5 to 1.2 
    Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
    z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh)/I_toinf
    dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)

    Phi = dPhi_mesh/bin_wid

    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(bin_cen)
    ys = np.log10(Phi_obs)
    y_model = np.log10(Phi_DO)
    y_err = np.log10(Phi_err)
    Chi2_L = np.sum( pow((ys - y_model)/y_err, 2))
    # print('in model -.5*(Chi2_L+Chi2_M)',-.5*(Chi2_L+Chi2_M))
    return {'M_BH':M_BH, 'MF':dn_MBH/dlog10M, 'MF_data':MF(M_BH), 'Chi2_M':Chi2_M,
            'M1450':bin_cen, 'LF':Phi_DO, 'LF_data':Phi_obs, 'Chi2_L':Chi2_L}