from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.

# LF bins same w/ Matsu18
N_lf = len(bin_cen)

# # initial: lambda_0=0.01, logM0 = 8.
# t_life, d_fit, l_cut, a = 20, .01, 1., 0.1 # f_seed = .01, log_prob= -9.89
# t_life, d_fit, l_cut, a = 25, .01, 1.2, -0.2 # f_seed = .1, log_prob= -36.82
# t_life, d_fit, l_cut, a = 30, .01, 1., -.2 # f_seed = 1., log_prob= -13.88
# # # best:
# # t_life, d_fit, l_cut, a = 19.8, 1.2e-3, 1.1557, -1.8e-01 # f_seed = 1.

# new_nbase initial: lambda_0=0.01, logM0 = 8.
t_life, d_fit, l_cut, a = 30, .01, 1., 0.1 # f_seed = .01, log_prob= -9.11
t_life, d_fit, l_cut, a = 35, .01, 1.2, -0.2 # f_seed = .1, log_prob= -16.37
t_life, d_fit, l_cut, a = 40, .01, .9, -.2 # f_seed = 1., log_prob= -11.45

t_life, d_fit, l_cut, a = 21.4, 0., .89, .15 # f_seed = 0.1, M1M0_e
t_life, d_fit, l_cut, a = 21.4, 0.001, .89, .15 # f_seed = 0.1, M1M0_d

# after more abin_mf bins, after calibration w/ direct sampling; 
# easycali initial: 
t_life, logd_fit, l_cut, a = 21.8, -1, .88, .19 # f_seed = 0.01
t_life, logd_fit, l_cut, a = 21.4, -3, .89, .15 # f_seed = 0.1
t_life, logd_fit, l_cut, a = 22.2, -2.98, .99, -.04 # f_seed = 1
# easycali best:
t_life, logd_fit, l_cut, a = 19.9, -1.08, .87, .17 # f_seed = 0.01, easycali
# t_life, logd_fit, l_cut, a = 19.6, -2.96, .87, .12 # f_seed = 0.1


x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

d_fit = pow(10.,logd_fit)
print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

# print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

i = 0
Chi2_min = 1e10; find_min = False

## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-T['t_col'])//t_life)
Nmax = Nt
dP_MBH = np.zeros(N_mf)
# t_tot = np.zeros(len(T))
while Nt>=0:
    t_point = tz - Nt*t_life
    T_seed = T[np.logical_and(t_point-t_life<T['t_col'],T['t_col']<=t_point)]
    dt_seed = t_point - T_seed['t_col']
    dP_MBH_prev = dP_MBH.copy()
    # new seeds (using 2d meshgrids)
    if len(T_seed):
        # z_mesh = kernelS_MBHmesh(abin_mf, T_seed['Mstar0'], dt_seed, l_cut)
        z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh,x0)/I_toinf
        dP_seed = Ps[1:,:] - Ps[:-1,:]
        dP_seed = np.nansum(dP_seed, axis=1)/len(T)
    else:
        dP_seed = np.zeros(N_mf)
    # prev BHMF
    dP_MBH_prev = np.exp(np.interp(np.log(M0s),np.log(M_BH),np.log(dP_MBH))) # M0 grow, consv_ratio=0
    # print(dP_MBH_prev)
    for iM1 in range(N_mf):
        # kernelS_MBH_M(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
        l1 = kernelS_MBH_M(M_BH[iM1],         M0s,t_life,1.,l_cut,d_fit)
        l2 = kernelS_MBH_M(M_BH[iM1]*(1.+eps),M0s,t_life,1.,l_cut,d_fit)
        dlnldlogM1 = np.log(l2/l1)/np.log10(1.+eps)
        klbd = dlnldlogM1 * pow(l1,a)*np.exp(-l1)/I_toinf
        klbd[l1<x0] = 0
        # print( 'sum of Plambda',np.nansum(klbd/dlnldlogM1),dlog10M )
        dP_MBH[iM1] = np.nansum(klbd*dP_MBH_prev*dlog10M0) + dP_seed[iM1]/dlog10M
    # print('each cycle: consv_ratio =',np.nansum(dP_MBH*dlog10M))

    Nt -= 1

dn_MBH = dP_MBH*n_base*f_bsm*f_seed*dlog10M

consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
print('in Phi_easyconv: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'Phi_easyconvMF'
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

# 10 N_M in 1e7-1e10 range, plus 12 N_L
index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
xs = M_BH[index][::len(index[0])//10]
ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
y_model = np.log10( (dn_MBH/dlog10M)[index][::len(index[0])//10] )
y_err = pow(np.log10(xs)-8.5,2)/3. + .2 # from 0.2 to 0.95
Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
z_mesh[z_mesh<x0] = x0
Ps = integral(a,z_mesh,x0)/I_toinf
dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)

# using abin_lf: test consv
# consv_ratio = np.nansum(dPhi_mesh)/(n_base*f_seed)
# print('LF consv_ratio:',consv_ratio)

Phi = dPhi_mesh/bin_wid

Phi *= 1e9
Phi_DO = Phi/corr_U14D20(bin_cen)
Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

T = Table(
    [bin_cen,Phi_obs,Phi_DO,Phi,Chi2*np.ones(N_lf)],
    names=('bin_cen','Phi_obs','Phi_DO','Phi','Chi2')
)

LFname = z6datapre+'Phi_easyconvLF'
ascii.write(T, LFname,
            formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            overwrite=True)

if np.nanmin([Chi2, Chi2_min]) == Chi2:
    find_min = True
    Chi2_min = Chi2
    T_min = T
    LFname_min = LFname

print('Chi2_M=',Chi2_M,'Chi2_L=',Chi2)
print('log_prob=',-.5*(Chi2_min*(len(Phi_obs)-1)+Chi2_M))
#, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)

print('time=',time.time()-t1)