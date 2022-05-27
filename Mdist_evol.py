from PYmodule import *
from PYmodule.l_intg import *
from datetime import datetime
# inherited from sk1:BHMF_prev/plotBHMF.py, dist_grow.py;
# plot histogram version of seed BHMF, after grow, Willott2010 curve.
t1 = time.time()

z1 = 6
t_end = t_from_z(z1)

T = Ts[0][0]
Nsite = len(T)

f_bsm = 1.
n_base = n_base[0]
N_concatenate = int(1e1) # Nsite * N_concatenate samples
N = Nsite*N_concatenate

# 不同f_seed 不影响BH growth 只是参数不同
# nbase1e-03,4pr8 nsteps=2e4 的 best fits
# 用f_seed=0.01 相当于1e2*N_BH sample
t_life, logd_fit, l_cut, a = 21.8, -3., 0.88, 0.19 # f_seed = 0.01
t_life, logd_fit, l_cut, a = 21.4, -3, .89, .15 # f_seed = 0.1
t_life, logd_fit, l_cut, a = 22.2, -2.98, .99, -.04 # f_seed = 1

# easycali bests; f_seed=0.01 still error, because of logd_fit=-1 not exp growth; debug??
t_life, logd_fit, l_cut, a = 19.9, -1.08, .87, .17; f_seed = 0.01
# t_life, logd_fit, l_cut, a = 19.6, -2.96, .87, .12; f_seed = 0.1

d_fit = pow(10.,logd_fit)
t_life *= Myr

# table stores the cdf of lambda
x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x,x0)/I_toinf

Nt = int(np.max((t_end-T['t_col'])//t_life)+1)

# # N_bri: count active bright quasar (L>L_bright) numbers from all sample at time t
# # potentially comparable with luminous quasar density Wang+2019b
# N_bri = [0]; ts = [0]; L_bright = 1e47

# print('Mstar0',T['Mstar0'])
# print("Nsite",Nsite, 't:',t/Myr,'t_end:',t_end/Myr, 'Nt:',Nt)

prex = '../4pevol/distf{0:d}N{1:d}_'.format(int(abs(np.log10(f_seed))),int(np.log10(N)))+datetime.now().strftime('%m%d%H%M_')
Mfname = prex+'Mevol.txt'
lfname = prex+'levol.txt'
z6fname = prex+'BHatz6.txt'

with open(Mfname,'w') as fM, open(z6fname, 'w') as fz6, open(lfname, 'w') as fl:
    fz6.write('{0:10s}{1:10s}{2:10s}{3:10s}\n'.format('M1','ls','L1','M1450_1'))
    for i_concatenate in range(N_concatenate):
        M0 = T['Mstar0']; t = T['t_col']; 
        M1s = [M0]; ts =  [t/Myr]; l1s = [.01*np.ones(len(M0))]
        for i_ in range(Nt):
            if np.any(t + t_life > t_end): 
                dt = t_life * np.ones(Nsite)
                dt[ np.where(t + t_life > t_end) ]  = t_end - t[ np.where(t + t_life > t_end)] 
                t = t + dt
            else:
                t = t + t_life
                dt = t_life
            uniform_a = np.random.uniform(size=Nsite)
            ls = np.zeros(Nsite)
            for i in range(Nsite):
                ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
            ls = ls*l_cut

            M1 = M1M0_d(M0,ls,dt,d_fit)
            # M1 = M1M0_e(M0,dt,ls) #试了 MF还是对不上

            M0 = M1
            L1 = L_M(M1,ls)
            M1450_1 = M1450_Lbol(L1)
            M1s.append(M1); ts.append(t/Myr); l1s.append(ls)
            # N_bri.append(len(np.where(L1>=L_bright)[0]))
        # print(np.allclose(t_end,t))

        M1s = np.array(M1s); l1s = np.array(l1s)
        M_low = 1e8
        index = np.where(M1>M_low)
        np.savetxt(fM, M1s.transpose()[index], fmt='%10.3e')
        np.savetxt(fl, l1s.transpose()[index], fmt='%10.3e')
        np.savetxt(fz6, np.array([M1,ls,L1,M1450_1]).transpose(), fmt='%10.3e')

pyname = sys.argv[0]
print(pyname,' time: ',time.time()-t1)
print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))

T_z6 = ascii.read(prex+'BHatz6.txt', guess=False, delimiter=' ')
M1, ls, L1, M1450_1 = T_z6['M1'], T_z6['ls'], T_z6['L1'], T_z6['M1450_1']

# MF histogram comparing w/ Phi_easy 
hist, bin_edges = np.histogram(M1,bins=abin_mf,density=False)
hist = hist*n_base*f_bsm*f_seed/N/dlog10M
ascii.write(Table([np.log10(M_BH),hist,MF(M_BH)]), 
z6datapre+'Mdist_evolMFz6.txt',
names=['M_BH','hist','W10'],
formats={'M_BH':'10.2e','hist':'10.2e','W10':'10.2e'},
overwrite=True)

# MF histogram comparing w/ Phi_easy 
hist, bin_edges = np.histogram(M1450_1,bins=bin_edg,density=False)
hist = hist*n_base*f_bsm*f_seed/N/bin_wid*1e9
hist_DO = hist/corr_U14D20(bin_cen)
ascii.write(Table([bin_cen,Phi_obs,hist_DO,hist]), 
z6datapre+'Mdist_evolLFz6.txt',
names=['M1450','Phi_obs','hist_DO','hist'],
formats={'M1450':'10.2f','hist':'10.2e','hist_DO':'10.2e','Phi_obs':'10.2e'},
overwrite=True)

exit(0)
# write lambda histogram: lbin, hist_tot, hist_L46, Pana
lbin = np.linspace(-2,1.2,num=20)
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf

L_limit = 1e46
index = np.where(L_limit<L1)
M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
print('len(ls) after selection:',len(ls_),' / ',Nsite*N_concatenate)
print('len of all L1',len(L1))
hist_, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)

strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
# lambda histogram file
ascii.write(Table([bin_edges[:-1],hist_/len(ls),hist/len(ls),Pana[1:]-Pana[:-1]]), 
prex+strhist_+'.txt',
names=['log_l',strhist_,'hist_tot','ana'],
formats={'log_l':'10.2f',strhist_:'10.2e','hist_tot':'10.2e','ana':'10.2e'},
overwrite=True)