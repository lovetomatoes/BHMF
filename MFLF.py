from PYmodule import *

def dP_dloglambda(log_lbd,mu,sigma):
    return 1./np.sqrt(2*np.pi)/sigma * np.exp(- pow(log_lbd-np.log10(mu),2)/(2.*pow(sigma,2)))

def MF_W10(M,Phi_star,M_star,alpha): # Willot+ 2010 BH mass function
    return Phi_star*pow(M/M_star,-alpha)*np.exp(-M/M_star)

def Mass_from_Lbol(L,Edd_ratio):
    return L/Edd_ratio/1.25e38

# Torres-Alba et al. 2020; constant and differential correction (Vito+18; Ueda+14)
# not used 
f_obsc_const = .8
def LF_M1450_DO_ana(M): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018; alpha changed; normalization enhanced
    Phi_M_star = 1.09e-8
    M_star = -24.9
    alpha  = -1.76; beta = -2.73
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) / (1-f_obsc_const)

x = np.linspace(-30,-22,num=40)
z = 6
y = LF_M1450(x,z)
y_DO07 = np.zeros(len(x))
y_DO20 = np.zeros(len(x))
for i in range(len(x)):
    y_DO07[i] = y[i]*corr_U14H07(x[i])
    y_DO20[i] = y[i]*corr_U14D20(x[i])
plt.figure(figsize=(10,8),dpi=400)
plt.plot(x,y*1e9,c='purple')
plt.plot(x,LF_M1450_CO(x)*1e9)
plt.plot(x,LF_M1450_DO_ana(x)*1e9,'--',label='DO ana')
plt.plot(x,y_DO07*1e9,'--',label='DO H07')
plt.plot(x,y_DO20*1e9,'--',label='DO D20')
plt.yscale('log')
plt.ylim(1e-3,1e3)
plt.xlim(x[-1],x[0])
plt.grid(True)
plt.legend(loc='best',fontsize=fslegend)
plt.savefig('../magFunction.png')