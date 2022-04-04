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
x = np.logspace(np.log10(x0),1.2,num=1000)
Pa_x = integral(a,x)/I_toinf

N_BH = int(1e5)
N_BH = int(1e4)

Nt = int(Dt/t_life+1)

M1s = [M0*np.ones(N_BH)]
zs = [z0]
L1s = [L_M(M0,.01)*np.ones(N_BH)]

t = t0

for i in range(Nt):
    if t + t_life > t_end:
        dt = t_end - t
        t = t_end
    else:
        t = t + t_life
        dt = t_life
    aaa = np.random.uniform(size=N_BH)
    xx,yy = np.meshgrid(aaa,Pa_x)
    # the first yy>xx
    ls = x[(yy>xx).argmax(axis=0)]

    M1 = M1M0(M0,dt,ls)
    M0 = M1
    L1 = L_M(M1,ls)
    M1s.append(M1)
    L1s.append(L1)
    zs.append(z_tH(t/Myr))

M1s = np.array(M1s)
L1s = np.array(L1s)

print('time after evol:',time.time()-t1)

index = np.logical_and(1e9<M0,M0<1e10)
index = np.nonzero(index)[0]
print('index',index, len(index))

prex = '../BH'
plt.figure(figsize=(10,10),dpi=400)
print('M1s.transpose()[index]',M1s.transpose()[index].shape)
print('M1s.transpose()[0]',M1s.transpose()[0].shape)
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
plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_M.png')


index = np.logical_and(1e9<M0,M0<1e10)
index = np.nonzero(index)[0]
print(len(index))

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
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_L.png')