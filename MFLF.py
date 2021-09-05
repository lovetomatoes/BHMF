from PYmodule import *

def dP_dloglambda(log_lbd,mu,sigma):
    return 1./np.sqrt(2*np.pi)/sigma * np.exp(- pow(log_lbd-np.log10(mu),2)/(2.*pow(sigma,2)))

def dndlogM_fit(M,Phi_star,M_star,alpha):
    return Phi_star*pow(M/M_star,-alpha)*np.exp(-M/M_star)

def M_L(L,Edd_ratio):
    return L/Edd_ratio/1.25e38

def dndlogL_fit(L):
    f_duty = 0.75
    mu = 0.6
    sigma = 0.3
    Phi_star = 1.23e-8
    M_star = 2.24e9
    alpha = 1.03

    Intg_0 = 1.e-9
    Intg = 0
    log_lbd = np.log10(mu)-2*sigma # initial lbd
    dlog_lbd = 0.0001
    i = 0
    while pow(Intg-Intg_0,2)/pow(Intg_0,2) > 1e-10:
        Intg_0 = Intg
        # print(dndlogM_fit(M_L(L,10**log_lbd),Phi_star,M_star,alpha))
        Intg += dndlogM_fit(M_L(L,10**log_lbd),Phi_star,M_star,alpha) * dlog_lbd * dP_dloglambda(log_lbd,mu,sigma)
        log_lbd += dlog_lbd
        i+=1
    # print("final Eddington ratio:",10**log_lbd)
    return f_duty*Intg


def LF(l): # dn/dlogL in Mpc^-3 dex^-1
    Phi_M_star = 1.14e-8
    M_star = -25.13
    alpha  = -1.5; beta = -2.81
    f_bol = 4.4
    Phi_L_star = Phi_M_star * 2.5
    L_star = pow(10,-.4*(M_star-34.1)) * 3e18/1450 *1e7 * f_bol
    L_1 = pow(10,-.4*(-27.2-34.1)) * 3e18/1450 *1e7 * f_bol 
    L_2 = pow(10,-.4*(-20.7-34.1)) * 3e18/1450 *1e7 * f_bol 
    # print('break L',L_star/W37, 'Phi_L_star', Phi_L_star)
    t = (np.log10(l) - np.log10(L_1)) / (np.log10(L_2) - np.log10(L_1))
    return Phi_L_star/( pow(l/L_star,-(alpha+1)) + pow(l/L_star,-(beta+1)) ) * (2*(1-t)+3*t) 

# x = np.logspace(45,48,num=40)
# # x = np.logspace(48.9,50)
# y = LF(x)
# plt.plot(x/W37,y)
# y_fit = np.ones(len(x))
# for i in range(len(x)):
#     y_fit[i] = dndlogL_fit(x[i])
#     # print(y[i],y_fit[i])
#     # print(M_L(x[i],8))
# plt.plot(x/W37,y_fit)
# # plt.xlim(10,1e4); plt.ylim(1e-10,1e-6)
# plt.xscale('log'); plt.yscale('log')
# plt.grid(True)
# plt.savefig('./LF.png')

f_obsc = .8

def LF_M1450_DO_ana(M): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    Phi_M_star = 1.09e-8
    M_star = -24.9
    alpha  = -1.76; beta = -2.73
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) / (1-f_obsc)

def LF_M1450_DO(M): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    Phi_M_star = 1.09e-8
    M_star = -24.9
    alpha  = -1.23; beta = -2.73
    M1450_0 = M1450_Lbol(3e47)
    print('M1450_0',M1450_0)
    ff = 5 + 20*(M-M1450_0)/7.5 #if M>M1450_0 else f_obsc
    # ff = f_obsc * pow (4./3., M-M1450_0)
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) * ff

def LF_M1450_CO(M): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    Phi_M_star = 1.09e-8
    M_star = -24.9
    alpha  = -1.23; beta = -2.73
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) / (1-f_obsc)

plt.figure(figsize=(10,8),dpi=400)
x = np.linspace(41,47,num=30)
phi_4375_0 = .43
beta = .24
phi_min = .2
eps = .17
a1 = .48
phi_4375_z = phi_4375_0*pow(1.+1., a1)
y = phi_4375_z - beta*(x-43.75)
plt.plot(x,y)
plt.ylim(0,1)
plt.xlim(x[0],x[-1])
plt.grid(True)
plt.legend(loc='best',fontsize=fslegend)
plt.savefig('./obsc.png')

x = np.linspace(-30,-22,num=40)
z = 6
y = LF_M1450(x,z)
y_DO = np.zeros(len(x))
for i in range(len(x)):
    y_DO[i] = LF_M1450_DO(x[i])
print('ratio',LF_M1450_DO_ana(-22)/LF_M1450_CO(-22))
print('3 order of magnitudes ',10**(.4*7.5*.53))
plt.figure(figsize=(10,8),dpi=400)
plt.plot(x,y*1e9,c='purple')
plt.plot(x,LF_M1450_CO(x)*1e9,c='g')
plt.plot(x,LF_M1450_DO_ana(x)*1e9,'--',c='C1',label='DO ana')
plt.plot(x,y_DO*1e9,'--',label='DO')
plt.yscale('log')
plt.ylim(1e-3,1e3)
plt.xlim(x[-1],x[0])
plt.grid(True)
plt.legend(loc='best',fontsize=fslegend)
plt.savefig('./magFunction.png')

# lbd = 0.8
# for i in range(len(x)):
#     Phi_star = 1.23e-8
#     M_star = 2.24e9
#     alpha = 1.03
#     print(dndlogM_fit(M_L(x[i],.8),Phi_star,M_star,alpha))