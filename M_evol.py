from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

z0 = 35.
t0 = t_from_z(z0)
z1 = 6.
t_end = t_from_z(z1)

Dt = t_end - t0
M0 = 1e3

f_seed = 0.01
t_life, d_fit, l_cut, a = 1.68814272e+01, 7.31406767e-03, 1.02157385e+00, 1.46099557e-01 # f_seed = .01
t_life *= Myr

# table stores the cdf of lambda 
I_toinf = integral_toinf(a)
x = np.logspace(np.log10(lambda_0),1.2,num=200)
Pa_x = integral(a,x)/I_toinf

N_BH = int(1e5)
# N_BH = int(1e1)

Nt = int(Dt/t_life+1)

M1s = [M0*np.ones(N_BH)]
zs = [z0]
L1s = [L_M(M0,.01)*np.ones(N_BH)]

# N_bri: count active bright quasar (L>L_limit) numbers from all sample at time t
# potentially comparable with luminous quasar density Wang+2019b
N_bri = [0]; ts = [0]; L_bright = 1e47

t = t0
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
    ls = ls*l_cut
    # xx,yy = np.meshgrid(uniform_a,Pa_x)
    # # argmax: the first yy>xx
    # ls = x[(yy>xx).argmax(axis=0)]

    M1 = M1M0_d(M0,ls,dt,d_fit)
    M0 = M1
    L1 = L_M(M1,ls)
    M1s.append(M1)
    L1s.append(L1)
    zs.append(z_tH(t/Myr))
    ts.append(t/Myr)
    N_bri.append(len(np.where(L1>=L_bright)[0]))

print('time after evol:',time.time()-t1, 'Nt=%d'%Nt)

print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))

prex = z6datapre+'/f{0:.0e}N{1:d}'.format(f_seed,int(np.log10(N_BH)))
prex = '../'
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
x = np.logspace(np.log10(lambda_0),1.2,num=len(lbin)) # for Pana
Pana = integral(a,x)/I_toinf

index = np.where(L_limit<L1)
M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
print('len(ls) after selection:',len(ls_),' / ',N_BH)
hist_, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)

strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
# histogram file
ascii.write(Table([bin_edges[:-1],hist_/len(ls),hist/len(ls),Pana[1:]-Pana[:-1]]), 
prex+strhist_+'.dat',
names=['log_l',strhist_,'hist_tot','ana'],
formats={'log_l':'10.2f',strhist_:'10.2e','hist_tot':'10.2e','ana':'10.2e'},
overwrite=True)