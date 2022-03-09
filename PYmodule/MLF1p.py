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

d_fit = 0.
l_cut = .7 # l_cut=2., l_cut' = l_cut/2; M=M0=1e7 grow as Eddington
a = .3

def lnlike(theta):
    t_life = theta
    t_life *= Myr
## --------- Mass Function ---------
    tz = t_from_z(z)
    dn_MBH = np.zeros(N_mf)
    Nt = np.max((tz-T['t_col'])//t_life)
    print('Nt=',Nt, 'N latest seed to z6', np.min((tz-T['t_col'])//t_life))
    dP_MBH = np.zeros(N_mf)
    # t_tot = np.zeros(len(T))
    t_discrete = 0; t_mf = 0 # recording computation time
    while Nt>=0:
        t_point = tz - Nt*t_life
        T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
        dt_seed = t_point - T_seed['t_col']
        dP_MBH_prev = dP_MBH.copy()
        # # #------------------
        # # each bin
        # for ibin in range(N_mf):
        #     t0 = time.time()
        #     # new seeds
        #     if len(T_seed):
        #         # #----------- Schechter lbd -----------
        #         x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
        #         x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
        #         x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
        #         dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)
        #         dP_seed = np.nansum(dP_seed)/len(T)
        #     else:
        #         dP_seed = 0.
        #     t_discrete += time.time() - t0
        #     t0 = time.time()
        #     # prev BHMF
        #     # #----------- Schechter lbd -----------
        #     x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, 1., l_cut, d_fit)
        #     x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, 1., l_cut, d_fit)
        #     x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
        #     dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
        #     t_mf += time.time() - t0
        
        # using 2d meshgrids
        # new seeds
        if len(T_seed):
            z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
            z_mesh[z_mesh<0] = 0.
            dP_seed = special.gammainc(a,z_mesh[1:,:]) - special.gammainc(a,z_mesh[:-1,:])
            dP_seed = np.nansum(dP_seed, axis=1)/len(T)
        else:
            dP_seed = 0.
        # prev BHMF
        z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)
        z_mesh[z_mesh<0] = 0.
        # _, dP_prev_mesh = np.meshgrid(M_BH,bin_left)
        #  *dP_prev_mesh same as * dP_MBH_prev
        dP_MBH = np.nansum((special.gammainc(a,z_mesh[:,:-1])-special.gammainc(a,z_mesh[:,1:]))*dP_MBH_prev, axis=1) + dP_seed
        #-----------------------

        Nt -= 1
    dn_MBH = dP_MBH*n_base*f_bsm
    print('t_discrete=%.1f'%t_discrete, 't_mf=%.1f'%t_mf)
    consv_ratio = np.nansum(dn_MBH)/n_base
    print('MF consv_ratio',consv_ratio)
    # wli: too loose constraint!
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta)
        print('consv_ratio: ',consv_ratio)
        return -np.inf

    # # prev 4 points:
    # T_MF = Table([M_BH, dn_MBH/dlog10M],names=('M_BH','Phi'))
    # T_log = np.log10(M_BH)
    # T_Phi = T_MF['Phi']
    # ys = []
    # y_model =  np.log10(MF(pow(10.,logMs)))
    # # print(y_model)
    # y_err = np.array([2., 1., .4, 2.])
    # for logM0 in logMs:
    #     i = np.max(np.where(T_log<logM0))
    #     t = (logM0 - T_log[i])/ (T_log[i+1]-T_log[i])

    #     y1 = np.log10(T_Phi[i]*(1-t) + T_Phi[i+1]*t)
    #     if not np.isfinite(y1):
    #         print('theta=',theta)
    #         print('inf or nan? y1=',y1)
    #         return -np.inf
    #     ys.append(y1)
    # ys = np.array(ys)

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
    print(Phi)
    z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
    P_mesh = special.gammainc(a,z_mesh[:-1,:])-special.gammainc(a,z_mesh[1:,:])
    n_mesh, _ = np.meshgrid(dn_MBH,bin_wid)
    # print('len(bin_wid),len(M_BH)',len(bin_wid),len(M_BH),z_mesh.shape,P_mesh.shape,n_mesh.shape)
    dPhi_mesh = np.nansum(dn_MBH*P_mesh,axis=1)
    Phi = dPhi_mesh/bin_wid
    print(Phi)

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
# range7: 1e1<t_life<200 and .1<l_cut<3. and 0.<a<1.:

def lnprior(theta):
    t_life = theta
    if 1e1<t_life<200.:
        return 0.0
    else:
        return -np.inf

def lnprobab(theta):
    lp = lnprior(theta)
    # print('theta',theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)