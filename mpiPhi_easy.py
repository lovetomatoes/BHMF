from PYmodule import *
import sys 
sys.path.append('/gpfs/share/home/1801110214')
from mpi4py import MPI

# Phiz6 paras: f_0, t_life, lbd ~Schechter(l_cut,a)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.


T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)

t_range = [120.*Myr]*3
f_range = [1.]
d_range = [.25]
l_range = [.9]
a_range = [.1]

# range3: 1e1<t_life<1e3 and .1<f_0<10. and 0<d_fit<1. and .1<l_cut<10. and a>0.:

t_range = np.append(np.arange(10., 100., 10.),np.arange(100., 1000.,100.))*Myr
# f_range = np.arange(.1, 10.,1.)
f_range = np.arange(.1, 2.3,.05)
d_range = np.arange(.1, 1., .1)
l_range = np.arange(.1, 10.,1.)
a_range = np.logspace(-3,1,num=10) # a>0 total P convergent


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
assert len(l_range)*len(a_range) == size
# 100 cores

i = 0
Chi2_min = 1e10; find_min = False; para_min =''; para_str = ''

a = a_range[ rank//len(l_range)]
l_cut = l_range[ rank%len(l_range) ]


for t_life in t_range:
    for f_0 in f_range:
        for d_fit in d_range:
            i = i+1
            # continue
        ## --------- Mass Function ---------
            dn_MBH = np.zeros(N_mf)
            Nt = np.max((tz-T['t_col'])//t_life)
            Nmax = Nt
            dP_MBH = np.zeros(N_mf)
            while Nt>=0:
                t_point = tz - Nt*t_life
                T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
                dt_seed = t_point - T_seed['t_col']
                dP_MBH_prev = dP_MBH.copy()
                for ibin in range(N_mf):
                    # new seeds
                    if len(T_seed):
                        x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                        x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                        x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                        dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)

                        dP_seed = np.nansum(dP_seed)/len(T)
                    else:
                        dP_seed = 0.
                    # prev BHMF
                    x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, f_0, l_cut, d_fit)
                    x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, f_0, l_cut, d_fit)
                    x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                    dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
                Nt -= 1
            dn_MBH = dP_MBH*n_base*f_bsm

            consv_ratio = np.nansum(dn_MBH)/n_base
            # print('MF conserved fraction=%.10f'%consv_ratio)
            T_MF = Table(
                [M_BH, dn_MBH/dlog10M, MF(M_BH)],
                names=('M_BH','Phi','W10_MF')
            )
            index = np.logical_and(T_MF['M_BH']>1e7, T_MF['M_BH']<1e10)
            Chi2_M = np.sum(pow(np.log(T_MF['Phi']/T_MF['W10_MF'])[index], 2))/(np.sum(index)-1)
            off_M = np.max(abs(np.log(T_MF['Phi']/T_MF['W10_MF'])[index]))

        # # --------- Luminosity Function ---------
            Phi = np.zeros(N_lf)
            Phi_csv = 0.
            T_MF['dn_MBH'] = T_MF['Phi']*dlog10M
            for ibin in range(N_lf):
                #----------- Schechter lbd -----------
                x0 = kernelS_M1450(bin_edg[ibin+1], T_MF['M_BH'], l_cut)
                x1 = kernelS_M1450(bin_edg[ibin],   T_MF['M_BH'], l_cut)
                dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

                dPhi = np.nansum(T_MF['dn_MBH']*dP_M1450)
                Phi_csv += dPhi
                Phi[ibin] = dPhi/bin_wid[ibin]
            # print('consv of dP_M1450:',Phi_csv/np.sum(n_base))

            Phi *= 1e9
            Phi_DO = Phi/corr_U14D20(bin_cen)
            Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
            off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

            if np.nanmin([Chi2, Chi2_min]) == Chi2:
                find_min = True
                Chi2_min = Chi2
                para_min = 't=%10.1e'%(t_life/Myr)+ \
                    ' f=%10.1f'%f_0+ \
                    ' d=%10.2f'%d_fit+ \
                    ' l=%10.1e'%l_cut+ \
                    ' a=%10.3f'%a+ \
                    ' al=%10.1f'%alpha+ \
                    ' CL=%10.1e'%Chi2+ \
                    ' oL=%10.1e'%off_L+ \
                    ' CM=%10.1e'%Chi2_M+ \
                    ' oM=%10.1e'%off_M

            para_str += '%10.1e'%(t_life/Myr)+ \
                    '%10.1f'%f_0+ \
                    '%10.2f'%d_fit+ \
                    '%10.1e'%l_cut+ \
                    '%10.3f'%a+ \
                    '%10.1f'%alpha+ \
                    '%10.1e'%Chi2+ \
                    '%10.1e'%off_L+ \
                    '%10.1e'%Chi2_M+ \
                    '%10.1e\n'%off_M

# all parameter and results from all cores
paras = np.array(comm.gather(para_str, root=0))
Chi2_mins = np.array(comm.gather(Chi2_min, root=0))
# parameter giving min Chi2_L and results from all cores
para_mins = np.array(comm.gather(para_min, root=0))

if rank == 0:
    with open("paras_range4.dat","w") as f:
        f.write('%10s'%'t_life'+'%10s'%'f_0'+'%10s'%'d_fit'+'%10s'%'l_cut'+ \
                '%10s'%'a'+'%10s'%'alpha'+ \
                '%10s'%'Chi2_L'+'%10s'%'off_L'+\
                '%10s'%'Chi2_M'+'%10s\n'%'off_M')
        for i in range(len(paras)):
            f.write(paras[i])

    print('Chi2_mins',Chi2_mins)
    print('Chi2_min=\n',np.nanmin(Chi2_mins), 'file:\n',para_mins[np.nanargmin(Chi2_mins)])

print(rank)
