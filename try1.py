from time import sleep
from PYmodule import *
# N1 = 1; N2 = 2; N3 = 3; N4 = 4; N5 = 5
# a = np.ones((N1,N2,N3,N4,N5))

print((t_from_z(6)-t_from_z(15))/Myr)
print((t_from_z(6)-t_from_z(5))/Myr)
print((t_from_z(5)-t_from_z(4))/Myr)
print(M1450_Lbol(1e13*Lsun))
print(Lsun)

print('z: 50 to 17.58    : %3.2f Myr', (t_from_z(17.58)-t_from_z(50))/Myr)
print('z: 17.58 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(17.58))/Myr)
print('z: 6 to 4    : %3.2f Myr', (t_from_z(4)-t_from_z(6))/Myr)

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