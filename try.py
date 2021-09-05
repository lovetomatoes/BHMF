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


# a = np.array([True,False,True])
# a = [ [[[] for i in range(4)] for j in range(2)] for k in range(5)]
# a[1][0][2].append(1)
# print(a)
print(r'$\mathrm{\Chi^2}$')
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

# read()
# main()

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