from PYmodule import *
# Phiz6 paras: f_0, t_life, lbd ~Schechter(l_cut,a)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)


t_life, d_fit, logM0 = 198.08025684,   0.2509987,    5.96

print('theta: ',t_life,d_fit,logM0, l_cut,a)

t_life *= Myr
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

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

index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
xs = M_BH[index]
ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
y_model = np.log10( (dn_MBH/dlog10M) [index] )
y_err = 1.
Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
Phi = np.zeros(N_lf)

T['dn_MBH'] = T['Phi']*dlog10M
z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
P_mesh = special.gammainc(a,z_mesh[:-1,:])-special.gammainc(a,z_mesh[1:,:])
dPhi_mesh = np.nansum(P_mesh*dn_MBH,axis=1)
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