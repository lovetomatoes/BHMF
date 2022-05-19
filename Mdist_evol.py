from PYmodule import *
from PYmodule.l_intg import *
from datetime import datetime
# inherited from sk1:BHMF_prev/plotBHMF.py, dist_grow.py;
# plot histogram version of seed BHMF, after grow, Willott2010 curve.
t1 = time.time()

z1 = 6
t_end = t_from_z(z1)

T = Ts[0][0][:5]
Nsite = len(T)

f_bsm = 1.
n_base = n_base[0]
N_concatenate = int(2e0) # Nsite * N_concatenate samples

# h_seed, bin_edges = np.histogram(M0,bins=abin_mf,density=False)
# Phi_seed = h_seed*n_base/float(Nsite)/dlog10M

# 不同f_seed 不影响BH growth 只是参数不同
# nbase1e-03,4pr8 nsteps=2e4 的 best fits
# 用f_seed=0.01 相当于1e2*N_BH sample
f_seed,t_life, logd_fit, l_cut, a = 0.01, 21.8, -1., 0.88, 0.19
t_life = 500

d_fit = pow(10.,logd_fit)
t_life *= Myr

# table stores the cdf of lambda
x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x,x0)/I_toinf

Nt = int(np.max((t_end-T['t_col'])//t_life)+1)

# #暂时不用 暂时只写z=6 Mdist Ldist
# #应该是 M0 * N_concatenate之类的数组 
# M1s = [M0*np.ones(N_BH)]
# zs = [z0]
# L1s = [L_M(M0,.01)*np.ones(N_BH)]

# # N_bri: count active bright quasar (L>L_bright) numbers from all sample at time t
# # potentially comparable with luminous quasar density Wang+2019b
# N_bri = [0]; ts = [0]; L_bright = 1e47

# print('Mstar0',T['Mstar0'])
# print("Nsite",Nsite, )

# print(t/Myr, 't_end:', t_end/Myr, 'Nt:',Nt)

prex = '../4p/dist'+datetime.now().strftime('%Y-%m-%d_%H-%M_')
Mfname = prex+'Mevol.dat'
Mfname = 'Mevol.txt'
z6fname = 'z6.txt'

with open(Mfname,'w') as fM, open(z6fname, 'w') as fz6:
    for i_concatenate in range(N_concatenate):
        M0 = T['Mstar0']; t = t = T['t_col']; 
        M1s = [M0]; ts =  [t/Myr]; L1s = [L_M(M0,.01)]
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
            M0 = M1
            L1 = L_M(M1,ls)
            M1s.append(M1); L1s.append(L1); ts.append(t/Myr)
            # N_bri.append(len(np.where(L1>=L_bright)[0]))
        
        M1s = np.array(M1s); L1s = np.array(L1s)
        M_low = 1e2
        index = np.where(M1>M_low)
        np.savetxt(fM, M1s.transpose()[index], fmt='%10.3e')
        np.savetxt(fz6, np.array([M1,ls,L1]).transpose(), fmt='%10.3e')

print('time after evol:',time.time()-t1, 'Nt=%d'%Nt)

print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))
# print('Mstar0',T['Mstar0'])

exit(0)
prex = z6datapre+'/f{0:d}N{1:d}'.format(abs(int(np.log10(f_seed))),int(np.log10(N_BH)))
prex = '../dist'+datetime.now().strftime('%Y-%m-%d_%H-%M_')
# prex = '../f{0:d}N{1:d}'.format(abs(int(np.log10(f_seed))),int(np.log10(N_BH)))

# z=6 BH mass, λ, L_bol
ascii.write(Table([M1, ls, L1]),prex+'BHatz6.dat',
names=['M1','ls','L1'],formats={'M1':'10.2e','ls':'10.2e','L1':'10.2e'},overwrite=True)

# all time samples
M1s = np.array(M1s)
L1s = np.array(L1s)
ascii.write(Table([zs,ts,N_bri]),prex+'tN_evol.dat',names=['zs','ts','N_bri'],
formats={'zs':'10.2f','ts':'10.2e','N_bri':'10.2e'},overwrite=True)

# BH mass evol. M growth tracks above M_low
M_low = 1e8
index = np.where(M1>M_low)
M1s = M1s.transpose()[index]
L1s = L1s.transpose()[index]
print('M1s.shape:',M1s.shape)
np.savetxt(prex+'Mevol.dat', M1s,fmt='%10.3e')

# M1450=-26, Lbol=1e47; M1450=-22, Lbol=2e45
L_limit = 1e46
print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))

lbin = np.linspace(-2,1.2,num=20)
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf

index = np.where(L_limit<L1)
M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
print('len(ls) after selection:',len(ls_),' / ',N_BH)
hist_, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)

strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
# lambda histogram file
ascii.write(Table([bin_edges[:-1],hist_/len(ls),hist/len(ls),Pana[1:]-Pana[:-1]]), 
prex+strhist_+'.dat',
names=['log_l',strhist_,'hist_tot','ana'],
formats={'log_l':'10.2f',strhist_:'10.2e','hist_tot':'10.2e','ana':'10.2e'},
overwrite=True)