from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

z0 = 35.
t0 = t_from_z(z0)
z1 = 6.
t_end = t_from_z(z1)

Dt = t_end - t0
M0 = 1e3

t_life, d_fit, l_cut, a = 59 ,  0., 1., 0.03
t_life *= Myr

# table stores the cdf of lambda 
I_toinf = integral_toinf(a)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x)/I_toinf

N_BH = int(1e7)
N_BH = int(1e5)

Nt = int(Dt/t_life+1)

M1s = [M0*np.ones(N_BH)]
zs = [z0]
L1s = [L_M(M0,.01)*np.ones(N_BH)]

# N_act: count active bright quasar (L>L_limit) numbers from all sample at time t
# potentially comparable with luminous quasar density Wang+2019b
N_act = [0]; ts = [0]; L_limit = 1e47

t = t0
L_low = 1e43
for i in range(Nt):
    if t + t_life > t_end:
        dt = t_end - t
        t = t_end
    else:
        t = t + t_life
        dt = t_life
    uniform_a = np.random.uniform(size=N_BH)
    ls = np.zeros(N_BH)
    for i in range(N_BH):
        ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
    ls = ls*l_cut/2.
    # xx,yy = np.meshgrid(uniform_a,Pa_x)
    # # argmax: the first yy>xx
    # ls = x[(yy>xx).argmax(axis=0)] / 2.

    M1 = M1M0(M0,dt,ls)
    M0 = M1
    L1 = L_M(M1,ls)
    M1s.append(M1)
    L1s.append(L1)
    zs.append(z_tH(t/Myr))
    ts.append(t/Myr)
    N_act.append(len(np.where(L1>=1e47)[0]))

# z=6 BH mass, λ, L_bol
ascii.write(Table([M1, ls, L1]),'../BHatz6.dat',names=['M1','ls','L1'],formats={'M1':'10.2e','ls':'10.2e','L1':'10.2e'},overwrite=True)
print('time after evol:',time.time()-t1)

# all time samples
M1s = np.array(M1s)
L1s = np.array(L1s)

# print('N_act at zs', N_act)
ascii.write(Table([zs,ts,N_act]),'../Nact_evol.dat',names=['zs','ts','N_act'],formats={'zs':'10.2f','ts':'10.2e','N_act':'10.2e'},overwrite=True)

# all BHs, all time: lambda distribution
abin = np.log10(x)
hist, bin_edges = np.histogram(np.log10(ls),bins=abin,density=False)
plt.figure(figsize=(10,8),dpi=400)
plt.scatter( bin_edges[:-1],hist/len(ls))
# print(np.sum(hist)/len(ls))
# print(np.sum((Pa_x[1:]-Pa_x[:-1])))
plt.plot(np.log10(x[:-1])/2.,(Pa_x[1:]-Pa_x[:-1]),c='C1')
plt.yscale('log')
plt.savefig('../Plambda_all_MBHevol.png')

# BH mass evol.
index = np.logical_and(1e9<M0,M0<1e10)
index = np.nonzero(index)[0]
print('index',index, len(index))

prex = '../BH'
plt.figure(figsize=(10,10),dpi=400)
# print('M1s.transpose()[index]',M1s.transpose()[index].shape)
# print('M1s.transpose()[0]',M1s.transpose()[0].shape)
for i in index:
    plt.plot(zs,M1s.transpose()[i])
index = np.argmax(M0)
plt.plot(zs,M1s.transpose()[index])
plt.yscale('log')
plt.xlim(10,z1)
plt.ylim(1e2,1e10)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_M.png')

# L_bol evol.
index = np.logical_and(1e9<M0,M0<1e10)
index = np.nonzero(index)[0]
# print(len(index))

plt.figure(figsize=(10,10),dpi=400)
for i in index:
    plt.plot(zs,L1s.transpose()[i])
index = np.argmax(M0)
plt.plot(zs,L1s.transpose()[index])
plt.yscale('log')
plt.xlim(10,z1)
plt.ylim(1e41,1e48)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{L_{bol}}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_L.png')