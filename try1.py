from time import sleep
from PYmodule import *
from PYmodule.l_intg import *
from PYmodule.MLF2p import *
from scipy.stats import norm, uniform
# N1 = 1; N2 = 2; N3 = 3; N4 = 4; N5 = 5
# a = np.ones((N1,N2,N3,N4,N5))

# print('most Nt:',(t_from_z(6)-t_from_z(20))/(200.*Myr))
# print((t_from_z(6)-t_from_z(5))/Myr)
# print((t_from_z(5)-t_from_z(4))/Myr)
# print(M1450_Lbol(1e13*Lsun))
# print(t_Edd/Myr)

t1 = time.time()
l = np.logspace(-2,1,num=1000)
a = 1.
# print(P_left_norm(a,l))
P_left_norm(a,l)
t2 = time.time()
print('t2-t1',t2-t1)
# print((special.gammainc(a,l)-special.gammainc(a,0.01))/(special.gammainc(a,10)-special.gammainc(a,0.01)))
(special.gammainc(a,l)-special.gammainc(a,0.01))/(special.gammainc(a,10)-special.gammainc(a,0.01))
t3 = time.time()
print('t3-t2',t3-t2)

l = [.01,.1,1,10]
a = 1.
print(P_left_norm(a,l))
print((special.gammainc(a,l)-special.gammainc(a,0.01))/(special.gammainc(a,10)-special.gammainc(a,0.01)))

x = np.arange(logx_min,logx_max,dlogx)
print(len(x))

exit(0)

t_range = [120.]
f_range = [1.]
d_range = [.25]
l_range = [.9]
a_range = [.1]


t1 =  time.time()
t_life = t_range[0]
f_0  = f_range[0]
d_fit = d_range[0]
l_cut = l_range[0]
a = a_range[0]

x = (t_life, f_0, d_fit, l_cut, a)
# print(lnlike(x))
print('time=',time.time()-t1)
exit(0)

print(np.sum(np.NINF+2))
for i in range(1):
    # try:
    #     a = np.log10(i)
    # except ZeroDivisionError:
    #     print('i is 0')
    a = np.log10(i)
    if np.isinf(a):
        print('isinf')
    print('print a=',a)
print(np.nanmax([np.inf,-np.inf, np.nan, 100]))
print(np.random.randn(2))
exit(0)

a = np.arange(1,10,1)
index = np.logical_and(3<a, a<5)
b = np.arange(2,11,1)
print(np.append(b[index],.1))
index_ = b<8
print(a, index,np.sum(index),index_,index*index_)

for L1450 in [1e40, 1e41, 1e42, 1e43]:
    print(M1450_Lbol(L1450*fbol_1450))
    print(-20.1-2.5*np.log10(L1450/1e44))
    # print(34.1-2.5*np.log10(L1450/(3e18/1450*1e7)))

print(t_Edd/Myr)
exit(0)
fig, ax = plt.subplots(1, 1, dpi=400)
x = np.linspace(uniform.ppf(0.01),
                uniform.ppf(0.99), 100)
y = uniform.pdf(x)
ax.plot(x, y,
       'r-', lw=5, alpha=0.6, label='uniform pdf')
rv = uniform()
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
r = uniform.rvs(size=10000)
ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.savefig('../unifrom.png')
hist, bin_edges = np.histogram(r,bins=np.append(x,1.1),density=False)
z = hist/len(r)/(x[1]-x[0])
ascii.write(Table([x,y,z]),'aaa',names={'x','y','z'},formats={'x':'10.5f','y':'10.5f','z':'10.5f'},overwrite=True)
ascii.write(Table([r]),'bbb',names={'r'},formats={'r':'10.5f'},overwrite=True)

fig, ax = plt.subplots(1, 1, dpi=400)
x = np.linspace(norm.ppf(0.01),
                norm.ppf(0.99), 100)
y = norm.pdf(x)
ax.plot(x, y,
       'r-', lw=5, alpha=0.6, label='norm pdf')
rv = norm()
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
r = norm.rvs(size=10000)
ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.savefig('../norm.png')
hist, bin_edges = np.histogram(r,bins=np.append(x,x[-1]**2/x[-2]),density=False)
z = hist/len(r)/(x[1]-x[0])
ascii.write(Table([x,y,z]),'normaaa',names={'x','y','z'},formats={'x':'10.5f','y':'10.5f','z':'10.5f'},overwrite=True)
ascii.write(Table([r]),'normbbb',names={'r'},formats={'r':'10.5f'},overwrite=True)

# Ls = np.logspace(41,48)
# # for L in Ls:
# M1 = M1450_Lbol(Ls)
# M2 = -25.-2.5*(np.log10(Ls/Lsun)-13.15364)
# ascii.write(Table([Ls,M1,M2]),'aaa',overwrite=True)

def filter(x,y):
    z = []
    for ax in x:
        z.append(math.floor(ax/y)*y)
    return z

xs = np.linspace(1,20,num=10000)
plt.plot(xs,filter(xs,3.0))
plt.savefig('../filter.png')
exit(0)

T = ascii.read('fort.10')
ls = T['ls']
print(np.max(ls),np.min(ls))
print('len(ls):',len(T))
bin_edge_l = np.linspace(np.min(ls),np.max(ls))
bin_cen_l = (bin_edge_l[1:]+bin_edge_l[:-1]) / 2.
bin_wid_l = bin_edge_l[1:]-bin_edge_l[:-1]
hist, bin_edges = np.histogram(ls,bins=bin_edge_l,density=True)
alpha = .5; beta = .4
plt.figure(figsize=(10,10),dpi=400)
plt.bar( bin_cen_l,hist*bin_cen_l,width=bin_wid_l,color='C'+str(1))
plt.plot(bin_cen_l, (bin_cen_l/beta)**alpha*np.exp(-bin_cen_l/beta)/special.gamma(alpha) )
plt.xscale('log')
plt.savefig('../P_lbd.png')
exit(0)

x0 = kernel_M1450(10,1e8,.6,.3)
x1 = kernel_M1450(-20,1e8,.6,.3)
print(special.erfc(x0)-special.erfc(x1))
print('M1450: L=1e41 ', M1450_Lbol(1e41))
a = .5
x1 = np.inf; x0 = 0
print('P_tot:',special.gammainc(a,x1) - special.gammainc(a,x0))
exit(0)

x = np.arange(-21,-29,-1)
y = corr_U14D20(x)
print(x,y)
plt.figure(figsize=(10,10),dpi=400)
z = 4 + 2.5*np.tanh(.5*(x+21))
plt.plot(x,y)
plt.plot(x,z)
plt.ylim(1,7)
plt.savefig('../corr.png')
exit(0)

print('z: 50 to 17.58    : %3.2f Myr', (t_from_z(17.58)-t_from_z(50))/Myr)
print('z: 17.58 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(17.58))/Myr)
print('z: 6 to 4    : %3.2f Myr', (t_from_z(4)-t_from_z(6))/Myr)

for M in [1e4,1e5,1e6,1e7,1e8,1e9,1e10]:
    print('M1450=%.1f'%M1450_Lbol(L_M(M,10)))
for i in range(100):
    print(i//10)
t1=time.time()
sleep(2)
t2=time.time()
print(t2-t1)

exit(0)
# ERDF: 和Schulze2015 作比
def Schechter(x,a):
    return pow(x,a)*np.exp(-x)

plt.figure(figsize=(10,10),dpi=400)
x = np.logspace(-2,0)
# SC: 
a = .5; l_cut = .1
y_SC = Schechter(x/l_cut,a)/special.gamma(a)/np.log10(math.e)
plt.plot(x,y_SC,label='Schechter')
# LN:
mu_fit = .05; sigma_fit = .3
y_LN = 1/np.sqrt(2*math.pi)/sigma_fit*np.exp(- np.log10(x/mu_fit)**2 /(2*sigma_fit**2))
plt.plot(x,y_LN,label='Log-norm')

# SC_M: 
a = .1; l_cut = .8/2.
y_SC = Schechter(x/l_cut,a)/special.gamma(a)/np.log10(math.e)
plt.plot(x,y_SC,label='Schechter + f(M)')

# plt.xlim(1e-2,1)
plt.xscale('log'); plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.xlabel(r'$\mathrm{\lambda}$',fontsize=fslabel)
# plt.ylabel(r'$\mathrm{Mpc^{-3} dex^{-1}}$',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig('../lbd_dist.png')

Nt = 3; Nf = 5
ts = .1*np.arange(Nt)
fs = np.arange(Nf)
for i in range(Nt*Nf):
    print(ts[i//Nf], fs[i%Nf])
exit(0)

# i = 0
# b = []
# for i1 in range(N1):
#     for i2 in range(N2):
#         for i3 in range(N3):
#             for i4 in range(N4):
#                 for i5 in range(N5):
#                     i += 1
#                     b.append(i-1)
#                     # print('i%(N2*N3);',(i-1)%(N3*N2), 'i2=',i2)
#                     # print('i%(N5*N4*N2*N3);',(i-1)%(N5*N4*N3*N2*N1), 'i1=',i1)
#                     # print((i-1)%(N5), 'i5=',i5)
#                     # print((i-1)//N5%N4, 'i4=',i4)
#                     # print((i-1)//N5//N4%N3, 'i3=',i3)
#                     # print((i-1)//N5//N4//N3%N2, 'i2=',i2)
#                     print((i-1)//N5//N4//N3//N2%N1, 'i1=',i1)


N1 = 2; N2 = 3; N3 = 4
a = np.ones((N1,N2,N3))

i = 0
b = []
N_tot = N1*N2*N3
p1 = np.zeros(N_tot); p2 = np.zeros(N_tot); p3 = np.zeros(N_tot)
for i1 in range(N1):
    for i2 in range(N2):
        for i3 in range(N3):
            p1[i] = i1
            p2[i] = i2
            p3[i] = i3
            b.append(i)
            i += 1

print(a.shape)
T = Table([p1,p2,p3,b],names=['p1','p2','p3','b'])
a = np.ones(4)
b[:] = a
# b = a
# b = a.copy()
a += 1
print(b)
print('1e9, 0.1 Eddington ratio: M1450=%.1f'%M1450_Lbol(L_M(1e9,.1)))
print('M1450=-25, edd=0.1, mass:%.1e'%M_L(Lbol_M1450(-25.),.1))
print('M1450=-29, edd=1., mass:%.1e'%M_L(Lbol_M1450(-29.),1.))
print('M1450=-29, edd=10, mass:%.1e'%M_L(Lbol_M1450(-29.),10.))

a = special.gamma([0.5, 1, 5])
a = special.gammainc(1,1)
print(a,1-1./np.exp(1))
# ascii.write(T,'../p1p2p3b',formats={'p1':'5.1e','p2':'5.1e','p3':'5.1e','b':'5.1e'},overwrite=True)