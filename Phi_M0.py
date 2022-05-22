from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
tz = 300*Myr

d_fit = 1e-3
l_cut = 1.
a = .5
t_life = 5


x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr
T = Ts[0][0]
T['t_col'] = 20*Myr*np.ones(len(T))
T['Mstar0'] = 1000
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
DT = 0
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
        DT+=dt_seed
    else:
        dP_seed = 0.
    # prev BHMF
    # z_mesh = kernelS_MBHmesh(M_BH, abin_mf, t_life, l_cut)
    DT+= t_life
    z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh,x0)/I_toinf
    dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed

    z_mesh_left = kernelS_MBH_M_mesh(bin_left, abin_mf, t_life, 1., l_cut, d_fit)
    z_mesh_left[z_mesh_left<x0] = x0
    Ps = integral(a,z_mesh_left,x0)/I_toinf
    dP_MBH_left = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed
    z_mesh_right = kernelS_MBH_M_mesh(bin_right, abin_mf, t_life, 1., l_cut, d_fit)
    z_mesh_right[z_mesh_right<x0] = x0
    Ps = integral(a,z_mesh_right,x0)/I_toinf
    dP_MBH_right = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed
    # 3 points averaging; tried 5 points -> dP_MBH, no use
    dP_MBH = (dP_MBH+dP_MBH_left+dP_MBH_right)/3.

    Nt -= 1
print('DT=',np.mean((DT-t_life)/Myr))
dn_MBH = dP_MBH*n_base*f_bsm*f_seed

consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
print('in Phi_easy: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'Phi_M0MF'
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

lbin = np.arange(-3,3,0.1)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf
with open(z6datapre+"Phi_M0ERDFz6.txt",'w') as f:
    np.savetxt(f, np.array([lbin[:-1],Pana[1:]-Pana[:-1]]).transpose(), fmt='%10.3e')

print('time=',time.time()-t1)
exit(0)

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
    [bin_cen,Phi_obs,Phi_DO,Phi],
    names=('bin_cen','Phi_obs','Phi_DO','Phi')
)

LFname = z6datapre+'LFIMF_easy_newnbase'
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

# print('time=',time.time()-t1)