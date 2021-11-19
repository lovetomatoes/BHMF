from numpy.lib.function_base import msort
from PYmodule import *
from scipy.optimize import curve_fit

global T
T = []
def read():
    name = "hmf_Pk.dat"
    global T
    if len(T)==0:
        print("reading")
        T = ascii.read(name, guess=False)


a = np.array([True,False,True])
print(isinstance(a, float))

a = np.linspace(-1,1,num=3) # a = [-1, 0, 1]
import numpy.ma as ma

b = np.ones(len(a))
b = ma.masked_where(a<0,b)
print(b,np.sum(b+1))
############### nanmin: ignore nan in array #############
print(np.nanmin(np.array([np.nan, 1, 3])))
print(np.nanmin([1e10,np.nan]))

a = np.linspace(0,5,6)
b = np.ones(len(a))
print(a,np.logical_or(a<=1., a>=4.))
c = ma.masked_where(np.logical_or(a<1, a>4), b)
print('ma.masked_where:', c)

# 我要的 6080个
f_duty = np.arange(.2, 1., .1) # .5
mu_fit = np.arange(.1, .5, .01) # .38
sigma_fit = np.arange(.01, 0.2, .01) # .12
# 暂时try 180 个
# f_duty = np.arange(.2, 1., 1) # .5
# mu_fit = np.arange(.1, .5, .02) # .38
# sigma_fit = np.arange(.01, 0.1, .01) # .12
print('len',len(f_duty)*len(mu_fit)*len(sigma_fit))

# a = [ [[[] for i in range(4)] for j in range(2)] for k in range(5)]
# a[1][0][2].append(1)
# print(a)
print((t_from_z(6)-t_from_z(4))/Myr)

a = np.array([0,1,2])
b = np.array([2,3,4])
print(a+b)
print(t_Edd/Myr)

def main():

    Js = []
    lam1 = .2
    lam2 = .3
    N = 100000
    Nh1s = poisson.rvs(lam1, size=N)
    Nh2s = poisson.rvs(lam2, size=N)
    Nhs = Nh1s + Nh2s
    # plt.hist(Nhs,density=True)
    # plt.show()

    sigma = .4
    sigma1 = .2
    sigma2 = .1
    log10Jmean = 1
    log10Jmean1 = 1.
    log10Jmean2 = 1.1

    for i in range(N):
        if Nh1s[i] ==1:
            rs = lognorm.rvs(sigma1*np.log(10), scale=10**log10Jmean1, size=1)
            Js.extend(rs)
        elif Nh2s[i] ==1:
            rs = lognorm.rvs(sigma2*np.log(10), scale=10**log10Jmean2, size=1)
            Js.extend(rs)
        else:
            Js.append(0)

    xs = np.logspace(-2,np.log10(max(Js)),num=50)
    ys = 1./np.sqrt(2*pi)/sigma * np.exp(- (np.log10(xs)-log10Jmean)**2/2./sigma**2) / (xs*np.log(10))
    y1s = 1./np.sqrt(2*pi)/sigma1 * np.exp(- (np.log10(xs)-log10Jmean1)**2/2./sigma1**2) / (xs*np.log(10))
    y2s = 1./np.sqrt(2*pi)/sigma2 * np.exp(- (np.log10(xs)-log10Jmean2)**2/2./sigma2**2) / (xs*np.log(10))

    plt.hist(Js,bins=20,density=True)
    plt.plot(xs, y1s * lam1*np.exp(-lam1) + y2s * lam2*np.exp(-lam2) )
    plt.plot(xs, y1s * lam1*np.exp(-lam1)*np.exp(-lam2) + y2s * lam2*np.exp(-lam1)*np.exp(-lam2) )
    plt.yscale('log')
    plt.ylim(bottom=1.e-10)
    plt.savefig('2shells_log.png')
    plt.show()

    plt.hist(Js,bins=20,density=True)
    # plt.plot(xs,ys*( np.exp(-lam2)*lam1*np.exp(-lam1) + np.exp(-lam1)*lam2*np.exp(-lam2) ))
    plt.plot(xs,y1s* lam1*np.exp(-lam1) + y2s * lam2*np.exp(-lam2) )
    plt.plot(xs, y1s * lam1*np.exp(-lam1)*np.exp(-lam2) + y2s * lam2*np.exp(-lam1)*np.exp(-lam2) )
    plt.xscale('log'); plt.yscale('log')
    plt.ylim(bottom=1.e-10)
    plt.savefig('2shells_log.png')
    plt.show()


# T = ascii.read('*.txt',guess=False) #没读到空行
# ascii.write(T,output='file1.txt',overwrite=True)

k1 = pow(10., 45.48-43.97)
print('k1 and k1_ueda', k1, K_AVE07(10**45.48),K_AVE20(10**45.48))
k2 = pow(10., 46.2-44.61)
print('k2 and k2_ueda', k2, K_AVE07(10**46.2),K_AVE20(10**46.2))


print('Lbol=1e13Lsun, M1450=',M1450_Lbol(1e13*Lsun))

# plt.figure(figsize=(10,8),dpi=400)
# x = np.logspace(7,15)*Lsun
# y07 = K_AVE07(x)
# y20 = K_AVE20(x)
# plt.plot(np.log10(x/Lsun),y07,label='y07')
# plt.plot(np.log10(x/Lsun),y20,label='y20')
# plt.xlim(7,15); plt.ylim(1,1000)
# plt.yscale('log')
# plt.legend(loc='best')
# plt.grid(True)
# plt.savefig('../Kx.png')

def f_obsc_U03(logLx): # Ueda 03; N_H > 22 fraction; as a func of Lx
    eta = 1.7
    phimax = (1+eta)/(3+eta)
    phi44 = .47
    beta = .1
    phi = min( phimax, max(phi44 - beta*(logLx-44), 0))
    f_sum = phi
    return f_sum

def corr_U03(M1450): # Ueda 03 + Shankar 09
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE07(L_bol)
    Lx = L_bol/f_bol
    eta = 1.7
    phimax = (1+eta)/(3+eta)
    phi44 = .47
    beta = .1
    phi = min( phimax, max(phi44 - beta*(np.log10(Lx)-44), 0))
    f_obsc_sum = phi
    return 1./(1-f_obsc_sum)
    return (1+eta/(1+eta)*phi)/(1-f_obsc_sum)

M1 = -27.2
print("correct facotr for for M1 = ",M1," 1/fobsc=",corr_U03(M1))
M1 = -20.7
print("correct facotr for for M1 = ",M1," 1/fobsc=",corr_U03(M1))

plt.figure(figsize=(10,8),dpi=400)
x = np.linspace(-19, -28)
y = np.zeros(len(x))

for i in range(len(x)):
    y[i] = corr_U03(x[i])
plt.plot(x,y,label='Ueda03')
for i in range(len(x)):
    y[i] = corr_U14D20(x[i])
plt.plot(x,y,label='Ueda14')
plt.xlim(-19, -28)
plt.xlabel("M1450",fontsize=fstick)
plt.ylabel("corr factor",fontsize=fstick)
plt.savefig('../NH_gtr22.png')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('../corr.png')

def tick_function(x):
    f_bol = KX_AVE20(10.**x)
    M1450 = M1450_Lbol(10.**x*f_bol)
    print('f_bol',f_bol,'M1450',M1450)
    return ["%.1f" % M for M in M1450] # M1450

def tick_M1450_MBH(x):
    M1450s = M1450_Lbol(L_M(x,.1))
    return ["%.1f" % M for M in M1450s]

Lx = 1e41 # erg/s
print('KX_AVE20(Lx)',KX_AVE20(Lx))
Lx = 1e43 # erg/s
print('M1450 from Lx = 1e43',M1450_Lbol(Lx*KX_AVE20(Lx)))
Lx = 3e45 # erg/s
print('M1450 from Lx = 3e45',M1450_Lbol(Lx*KX_AVE20(Lx)))

ratio = np.exp( (t_from_z(4.)-t_from_z(15.))*.3/t_Edd )
print('ratio',ratio)
print(ratio**.42)
for i in range(10):
    if i<5:
        continue
    else:
        print(i)
    print(i)
print('1e9, 0.1 Eddington ratio: L_bol=%.1e'%L_M(1e9,.1))
print('1e9, 0.1 Eddington ratio: M1450=%.1f'%M1450_Lbol(L_M(1e9,.1)))
print('M1450=-22, edd=0.1, mass:%.1e'%M_L(Lbol_M1450(-22.),.1))
print('M1450=-30, edd=0.1, mass:%.1e'%M_L(Lbol_M1450(-30.),.1))
# def kernel_MBH(Mgrow_ratio, dt, f_duty, mu, sigma_dex):
k1 = kernel_MBH1(1.001, 1e6, .1, .1, .1)
# def kernel_MBH2(M1, M0, dt, f_duty, mu, sigma_dex, eta8, delta):
k2 = kernel_MBH2(1e8,1e8/1.001,1e6,.1,.1,.1,.1,.32)
print('k1',k1,'k2',k2)
print('6-15',(t_from_z(6)-t_from_z(15))/Myr, '4-6',(t_from_z(4)-t_from_z(6))/Myr)

fig = plt.figure(figsize=(10,8),dpi=400)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

pre = '/Users/wli/plots/MF2e10_f0.70m0.21s0.15'
for a in ['e0.100d0.008alpha1.0','e0.100d0.020alpha1.0']:
    fname = pre+a
    T = ascii.read(fname, guess=False, delimiter=' ')
    x = T['M_BH']; y = T['dn_MBH']
    ax1.plot(x,y,label=a)
new_tick_locations = np.logspace(6,12,num=20)
# ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_M1450_MBH(new_tick_locations))
print(new_tick_locations,tick_M1450_MBH(new_tick_locations))

# ax1.set_xlim(41,47); ax1.set_ylim(0,1.)
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xscale('log')
ax1.set_yscale('log')
# ax2.set_xscale('log')
ax1.set_xlabel(r"$M_{BH}$",fontsize=fstick)
ax2.set_xlabel(r"M1450",fontsize=fstick)
ax1.set_ylabel("dn_MBH",fontsize=fstick)
plt.savefig('../mf.png')
exit(0)

fig = plt.figure(figsize=(10,8),dpi=400)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
# # Ueda 03 
# x = np.linspace(42,46) # log Lx
# y = np.zeros(len(x))
# for i in range(len(x)):
#     y[i] = f_obsc_U03(x[i])
# plt.plot(x,y,label='Ueda03')
# plt.xlim(42,46); plt.ylim(.3,.8)
# =========================

# # Ueda 14 ===============
x = np.linspace(41,47) # log Lx
y = np.zeros(len(x))
z = .5
for i in range(len(x)):
    y[i] = f_obsc_U14(x[i],z)
ax1.plot(x,y,label='Ueda14, z=.5')
z = 2
for i in range(len(x)):
    y[i] = f_obsc_U14(x[i],z)
ax1.plot(x,y,label='Ueda14, z=2')
z = 1.5
for i in range(len(x)):
    y[i] = f_obsc_U14(x[i],z)

print(y)
new_tick_locations = np.array([41,42,43,44,45,46,47])
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_locations))

ax1.plot(x,y,label='Ueda14, z=1.5')
# =========================
ax1.set_xlim(41,47); ax1.set_ylim(0,1.)
ax1.legend(loc='best')
ax1.grid(True)

ax1.set_xlabel(r"log10(Lx)",fontsize=fstick)
ax2.set_xlabel(r"M1450",fontsize=fstick)
ax1.set_ylabel("obscured fraction",fontsize=fstick)
plt.savefig('../NH_gtr22.png')


def lognorm_dist():
    a = np.zeros(5)
    b = np.ones(5)
    a = np.concatenate((a,b))
    b += b
    print(a,b[1:-1])
    sigma = .1
    log10Jmean = np.log10(.3)
    N = int(1e8)
    # rs = lognorm.rvs(sigma*np.log(10), scale=10**log10Jmean, size=N)
    # rs = lognorm.rvs(0.1*np.log(10), scale=0.3, size=1000) # center=0.3, scatter=0.1dex
    print('isf',lognorm(0.1*np.log(10), scale=0.3).isf(1e-8))
    # plt.hist(rs,bins=100,density=True)
    # xs = np.logspace(-1,0,num=50)
    # ys = 1./np.sqrt(2*pi)/sigma * np.exp(- (np.log10(xs)-log10Jmean)**2/2./sigma**2) / (xs*np.log(10))
    # plt.plot(xs,ys)
    # plt.plot(xs, lognorm(sigma*np.log(10), scale=10**log10Jmean).pdf(xs), '--')
    # plt.xscale('log'); plt.yscale('log')
    # plt.show()

lognorm_dist()


# Phi_M_star = 1.14e-8
# M_star = -25.13
# alpha  = -1.5; beta = -2.81
# f_bol = 4.4
# z = 6

# def Phi_M(m):
#     return Phi_M_star/ (pow(10., 0.4*(alpha+1)*(m-M_star)) + pow(10., 0.4*(beta+1)*(m-M_star)))

# def LF(l): # dn/dlogL in Mpc^-3 dex^-1
#     Phi_L_star = Phi_M_star * 2.5
#     L_star = pow(10,-.4*(M_star-34.1)) * 3e18/1450 *1e7 * f_bol
#     L_1 = pow(10,-.4*(-27.2-34.1)) * 3e18/1450 *1e7 * f_bol 
#     L_2 = pow(10,-.4*(-20.7-34.1)) * 3e18/1450 *1e7 * f_bol 
#     print('break L',L_star/W37, 'Phi_L_star', Phi_L_star)
#     t = (np.log10(l) - np.log10(L_1)) / (np.log10(L_2) - np.log10(L_1))
#     return Phi_L_star/( pow(l/L_star,-(alpha+1)) + pow(l/L_star,-(beta+1)) ) * (2*(1-t)+3*t) 

# xdata = np.array([1e45, 1e46, 1e47, 1e48])
# ydata = np.array([5e-7, 1.3e-7, 1e-8, 2e-10])

# x = np.logspace(45,48)
# y = LF(x)
# plt.plot(x/W37,y)
# # plt.scatter(xdata/W37,ydata,s=40)
# plt.xlim(10,1e4); 
# plt.xscale('log')
# plt.ylim(1e-10,1e-6); 
# plt.yscale('log')
# plt.grid(True)
# plt.savefig('./LF.png')