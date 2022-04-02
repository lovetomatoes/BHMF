from PYmodule import *
from PYmodule.l_intg import *
# Phiz6 paras: f_0, t_life, lbd ~Schechter(l_cut,a)

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.

# LF bins same w/ Matsu18
N_lf = len(bin_cen)

t_life, d_fit, logM0, l_cut, a = 80 ,  0.3 ,  8 ,  .5 ,  .5
t_life, d_fit, logM0, l_cut, a = 50 ,  0.3 ,  8 ,  .5 ,  .3
t_life, d_fit, logM0, l_cut, a = 10 ,  1.,  8 ,  1,  .5; #almost impossible
t_life, d_fit, logM0, l_cut, a = 150 ,  0.,  8 ,  .5,  .3
t_life, d_fit, logM0, l_cut, a = 120 ,  0.3,  8 ,  1.,  -.7

t_life, d_fit, logM0, l_cut, a, x0 = 150 ,  0.163,  8 ,  1.,  -.7, 0.1
t_life, d_fit, logM0, l_cut, a, x0 = 66.43062875,  0.06772173,  8 ,  1.,  -.7, 0.1
t_life, d_fit, logM0, l_cut, a, x0 = 77,  0.01,  8 ,  1.,  -.7, 0.1
# t_life, d_fit = 1.80438151e+02, 1.00071295e-02

t_life, d_fit, logM0, l_cut, a, x0 = 50 ,  0.2,  8 ,  1.,  -.1,  0.01
t_life, d_fit, logM0, l_cut, a, x0 = 100 ,  0.2,  8 ,  1.,  -.5,  0.01
t_life, d_fit, logM0, l_cut, a, x0 = 84 ,  0.27,  8 ,  .9,  .1,  0.01

# positive a, cut at x0
t_life, d_fit, logM0, l_cut, a, x0 = 80,   0.25,  8,   .9,  .1, 0.01
# best fit (t_life, d_fit)
t_life, d_fit, logM0, l_cut, a, x0 = 62,   0.01,  8,   .9,  .1, 0.01
# t_life, d_fit, logM0, l_cut, a, x0 = 74,   0.04,  8,   .9,  .1, 0.001

I_toinf =  integral_toinf(a)

print('t_life, d_fit, logM0, l_cut, a: ', t_life,', ',d_fit,', ',logM0,', ', l_cut,', ',a)

t_life *= Myr
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

i = 0
Chi2_min = 1e10; find_min = False

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
    # new seeds (using 2d meshgrids)
    if len(T_seed):
        z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit, logM0)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh)/I_toinf
        dP_seed = Ps[1:,:] - Ps[:-1,:]
        dP_seed = np.nansum(dP_seed, axis=1)/len(T)
    else:
        dP_seed = 0.
    # prev BHMF
    z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit, logM0)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh)/I_toinf
    dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed

    Nt -= 1

dn_MBH = dP_MBH*n_base*f_bsm

consv_ratio = np.nansum(dn_MBH)/n_base
print('in Phi_easy: MF conserved fraction=%.10f'%consv_ratio)
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

# 30 N_M in 1e7-1e10 range, plus 12 N_L
index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
xs = M_BH[index]
ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
y_model = np.log10( (dn_MBH/dlog10M) [index] )
y_err = 1.
y_err = (np.log10(xs)-8.5)**2/3 + .5 # from 0.5 to 1.2 
Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
Phi = np.zeros(N_lf)

T['dn_MBH'] = T['Phi']*dlog10M
z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
z_mesh[z_mesh<x0] = x0
Ps = integral(a,z_mesh)/I_toinf
dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)

Phi = dPhi_mesh/bin_wid

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


if find_min:
    print('log_prob=',-.5*(Chi2_min*(len(Phi_obs)-1)+Chi2_M))
    #, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)

print('time: ',time.time()-t1)