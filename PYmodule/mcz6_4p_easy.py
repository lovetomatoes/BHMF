from PYmodule import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.

t_range = [120.]
d_range = [.25]
l_range = [.9]
a_range = [.1]
# Chi2_M = 1.7; Chi2_L = 0.68 quite small.
# off_M = 1.6; off_L = 2.2

logMs = np.linspace(7,10,num=4)


T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

def lnlike(theta, logMs = logMs, z = int(6)):
    t_life, d_fit, l_cut, a = theta
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
                x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
                x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
                x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)
                dP_seed = np.nansum(dP_seed)/len(T)
            else:
                dP_seed = 0.
            # prev BHMF
            # #----------- Schechter lbd -----------
            x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, 1., l_cut, d_fit)
            x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, 1., l_cut, d_fit)
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
        return -np.inf

    T_MF = Table([M_BH, dn_MBH/dlog10M],names=('M_BH','Phi'))

    T_log = np.log10(M_BH)
    T_Phi = T_MF['Phi']
    ys = []
    y_model =  np.log10(MF(pow(10.,logMs)))
    # print(y_model)
    y_err = np.array([2., 1., .4, 2.])
    for logM0 in logMs:
        i = np.max(np.where(T_log<logM0))
        t = (logM0 - T_log[i])/ (T_log[i+1]-T_log[i])

        y1 = np.log10(T_Phi[i]*(1-t) + T_Phi[i+1]*t)
        if not np.isfinite(y1):
            print('theta=',theta)
            print('inf or nan? y1=',y1)
            return -np.inf
        ys.append(y1)
    ys = np.array(ys)
    # print(pow(10.,ys), -.5*np.sum( pow((ys - y_model)/y_err, 2)))
    return -.5*np.sum( pow((ys - y_model)/y_err, 2))


# range: 1e1<t_life<1e3 and d_fit>0. and .1<l_cut<10. and a>0.:
# range1: 1e1<t_life<1e3 and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range2: 1e1<t_life and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range3: 1e1<t_life<200 and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range4: 1e1<t_life<200 and 0<d_fit<1. and .1<l_cut<10. and a>0.:
# range5: 1e1<t_life<200 and 0<d_fit<.5 and .1<l_cut<10. and a>0.:

def lnprior(theta):
    t_life, d_fit, l_cut, a = theta
    if 1e1<t_life<200. and 0<d_fit<.5 and .1<l_cut<10. and a>0.:
        return 0.0
    else:
        return -np.inf

def lnprobab(theta):
    lp = lnprior(theta)
    # print('theta',theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)