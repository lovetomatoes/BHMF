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


Ts = [] # bsm=0,1 two files
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        T['Lbol_z'] = np.zeros(len(T))
        T['M1450_z'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['t_col'] = t_from_z(T['z_col'])
        TM.append(T)
        # print(np.min(T['z_col']),np.max(T['z_col']),np.min(T['Mstar0']),np.max(T['Mstar0']))
    Ts.append(TM)

# MF
abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)

t_range = np.array([100, 200, 500, 1000])*Myr
f_range = np.arange(.1, 1., .3)
l_range = np.logspace(-2, 1., num=4) # l_cut
a_range = np.array([.1, 1., 2., 3.]) # a>0 total P convergent
# print(len(t_range)*len(f_range)*len(l_range)*len(a_range) )
# print(a_range); exit(0)

t_range = np.arange(100,1000,100)*Myr
f_range = np.arange(.1, 1., .1)
l_range = np.append( [.01,.05], np.arange(.1,2.,.1))
a_range = np.arange(.1, 3., .1) # a>0 total P convergent
# print(len(t_range)*len(f_range)*len(l_range)*len(a_range) )
# exit(0)

t_range = [120.*Myr]*3
f_range = [1.]
d_range = [.25]
l_range = [.9]
a_range = [.1]
# Chi2_M = 1.7; Chi2_L = 0.68 quite small.
# off_M = 1.6; off_L = 2.2

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
assert len(t_range)*len(f_range) == size
# 24 cores

i = 0
Chi2_min = 1e10; find_min = False; LFname_min = ''; para_str = ''

t_life = t_range[ rank//len(f_range)]
f_0 = f_range[ rank%len(f_range) ]
for d_fit in d_range:
    for l_cut in l_range:
        for a in a_range:
            i = i+1
            # continue
        ## --------- Mass Function ---------
            dn_MBH = np.zeros(N_mf)
            for iM in range(N_Mh):
                for i_bsm in range(Nbsm):
                    T = Ts[iM][i_bsm]
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
                                # #----------- Schechter lbd -----------
                                # x0 = kernelS_MBH(bin_left[ibin] /T_seed['Mstar0'], dt_seed, f_0, l_cut)
                                # x1 = kernelS_MBH(bin_right[ibin]/T_seed['Mstar0'], dt_seed, f_0, l_cut)
                                x0 = kernelS_MBH_M(bin_left[ibin],  T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                                x1 = kernelS_MBH_M(bin_right[ibin], T_seed['Mstar0'], dt_seed, f_0, l_cut, d_fit)
                                x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)

                                dP_seed = np.nansum(dP_seed)/len(T)
                            else:
                                dP_seed = 0.
                            # prev BHMF
                            # #----------- Schechter lbd -----------
                            # x0 = kernelS_MBH(M_BH[ibin]/bin_right, t_life, f_0, l_cut)
                            # x1 = kernelS_MBH(M_BH[ibin]/bin_left,  t_life, f_0, l_cut)
                            x0 = kernelS_MBH_M(M_BH[ibin], bin_right, t_life, f_0, l_cut, d_fit)
                            x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  t_life, f_0, l_cut, d_fit)
                            x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                            dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
                        Nt -= 1
                    dn_MBH += dP_MBH*n_base[iM]*f_bsm[i_bsm]

            consv_ratio = np.nansum(dn_MBH)/np.sum(n_base)
            # print('conserved fraction=%.10f'%consv_ratio)
            # if consv_ratio<.9:
            #     print('conserved fraction=%.10f'%consv_ratio)
            T = Table(
                [M_BH, dn_MBH/dlog10M, MF(M_BH)],
                names=('M_BH','Phi','W10_MF')
            )
            index = np.logical_and(T['M_BH']>1e7, T['M_BH']<1e10)
            Chi2_M = np.sum(pow(np.log(T['Phi']/T['W10_MF'])[index], 2))/(np.sum(index)-1)
            off_M = np.max(abs(np.log(T['Phi']/T['W10_MF'])[index]))

        # # --------- Luminosity Function ---------
            Phi = np.zeros(N_lf)
            Phi_csv = 0.
            T['dn_MBH'] = T['Phi']*dlog10M
            for ibin in range(N_lf):
                #----------- Schechter lbd -----------
                x0 = kernelS_M1450(bin_edg[ibin+1], T['M_BH'], l_cut)
                x1 = kernelS_M1450(bin_edg[ibin],   T['M_BH'], l_cut)
                dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

                dPhi = np.nansum(T['dn_MBH']*dP_M1450)
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
    with open("paras.dat","w") as f:
        f.write('%10s'%'t_life'+'%10s'%'f_0'+'%10s'%'d_fit'+'%10s'%'l_cut'+ \
                '%10s'%'a'+'%10s'%'alpha'+ \
                '%10s'%'Chi2_L'+'%10s'%'off_L'+\
                '%10s'%'Chi2_M'+'%10s\n'%'off_M')
        for i in range(len(paras)):
            f.write(paras[i])

    print('Chi2_mins',Chi2_mins)
    print('Chi2_min=\n',np.nanmin(Chi2_mins), 'file:\n',para_mins[np.nanargmin(Chi2_mins)])