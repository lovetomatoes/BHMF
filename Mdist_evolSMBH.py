from PYmodule import *
from PYmodule.l_intg import *
from datetime import datetime
# inherited from Mdist_evol.py; select growth tracks of SMBH matching luminous quasar at z=7.642 
# ref: Wang+2021 Mbh = 1.6e9Msun, lambda = 0.67 

t1 = time.time()

z1 = 7.642
t_end = t_from_z(z1)

T = Ts[0][0]
Nsite = len(T)

f_bsm = 1.
n_base = n_base[0]
N_concatenate = int(1e2) # Nsite * N_concatenate samples
N = Nsite*N_concatenate

# 不同f_seed 不影响BH growth 只是参数不同; 影响BHMF hist normalization
# 用f_seed=0.01 相当于1e2*N_BH sample

# M0r8 bests:
# f1 [20.07157851 -2.98140382  0.89453609  0.12195823] -4.137038049964094
t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
# f2 [18.7555167  -1.2574505   0.87372563  0.20389703] -3.1054538991409824
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01

d_fit = pow(10.,logd_fit)
t_life *= Myr

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
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

prex = '../4p/distf{0:d}N{1:d}_'.format(int(abs(np.log10(f_seed))),int(np.log10(N)))+datetime.now().strftime('%m%d%H%M_')
prex = '../temp/distf{0:d}N{1:d}_'.format(int(abs(np.log10(f_seed))),int(np.log10(N)))
Mfname = prex+'Mevol.txt'
lfname = prex+'levol.txt'
# z6fname = prex+'BHatz6.txt'

with open(Mfname,'a') as fM, open(lfname, 'a') as fl:
    # fz6.write('{0:10s}{1:10s}{2:10s}{3:10s}\n'.format('M1','ls','L1','M1450_1'))
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
        # # select close to M = 1.6e9, l=0.67
        # # not able to do: t varying w/ seeds
        #     if t>=t_from_z(7.642) and t-dt<t_from_z(7.642):
        #         print("max M_BH at z=7.642 is %.1e"%np.max(M1))
        #         index = np.logical_and(np.logical_and(1e9<M1, M1<2e9), np.logical_and(ls>0.6,ls<0.7))
        #         if np.sum(index) ==0:
        #             continue
            M1s.append(M1); ts.append(t/Myr); l1s.append(ls)
            # N_bri.append(len(np.where(L1>=L_bright)[0]))
        # print(np.allclose(t_end,t))

        M1s = np.array(M1s); l1s = np.array(l1s)
        # select close to M = 1.6e9, l=0.67
        index = np.logical_and(np.logical_and(1e9<M1, M1<2e9), np.logical_and(ls>0.6,ls<0.7))
        np.savetxt(fM, M1s.transpose()[index], fmt='%10.3e')
        np.savetxt(fl, l1s.transpose()[index], fmt='%10.3e')
        # np.savetxt(fz6, np.array([M1,ls,L1,M1450_1]).transpose(), fmt='%10.3e')

pyname = sys.argv[0]
print(pyname,' time: ',time.time()-t1)
print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))