import numpy as np
from ctypes import * # c 类型库
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, vstack
import os
from scipy.stats import *
import time
import math
import numpy.ma as ma
from scipy import special

z4figpre = '../z4/figs/'
z4datapre = '../z4/data/'
z5figpre = '../z5/figs/'
z5datapre = '../z5/data/'
z6figpre = '../z6/figs/'
z6datapre = '../z6/data/'
datapre = '../data/'
figpre = '../figs/'

typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {'4':'Akiyama_18','5':'Niida_20','6':'Matsuoka_18'}

global G, h0, H0, Omega_m0, Omega_L0, m_H, mu, Ms, pi, km, pc, Myr, alpha_T
G, c, k_B, m_H = 6.67408e-8, 2.9979245e10, 1.38064852e-16, 1.66053904e-24
pi = 3.141593
mu = 1.2
Ms = 2.e33
Lsun = 3.828e33
pc = 3.e18
Mpc = 1.e6*pc
km = 1.e5
yr = 365*24*3600
Myr = 1.e6*(365*24*3600)
Omega_m0 = 0.307
Omega_L0 = 1 - Omega_m0
h0 = .677
H0 = h0*100*km/Mpc

t_Edd = 1./(4*pi*G/.4/(0.1*c))
fbol_1450 = 4.4

log10Ms = [11,12,13]
# n_base = [1.63,1.09e-01,4.02e-03,3.87e-05,1.07e-08]
n_base = [4.02e-03,3.87e-05,1.07e-08]
Mhpres = [datapre+'1e11',datapre+'1e12',datapre+'1e13']

f_bsm = [.6,.4]
Nbsm = 2 # control bsm case to use

W37 = 1e44

alpha_T = 2.324e4

fstick = 20
fstxt = 20
fslabel = 23
fstitle = 20
fslegend = 20

my_cmap = plt.get_cmap("viridis")
rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

# -----------------  binned data in Gpc^-3 mag^-1     ---------------------
# -----------------  '6': z=6 Matsuoka 2018; Table 4     ---------------------
# -----------------  '5': z=5 Niida 2020; Table 3 excluding m_i<23.1    ---------------------
# -----------------  '4': z=4 Akiyama 2018;     ---------------------    
bin_cen = {'6':np.array([-29,-27.5,-26.75,-26.25,-25.75,-25.25,-24.75,-24.25,-23.75,-23.25,-22.75,-22]),
           '5':np.append(np.arange(-28.63,-26.,.25),np.arange(-27.07,-23.5,.5)),
           '4':np.arange(-28.875, -21.625, .25)
           }
bin_wid = {'6':np.array([2, 1, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1]),
           '5':np.append(0.25*np.ones(len(np.arange(-28.63,-26.,.25))),
                         0.5*np.ones(len(np.arange(-27.07,-23.5,.5)))),
           '4':.25*np.ones(len(bin_cen['4']))
           }
bin_edg = {'6':np.append(bin_cen['6'] - bin_wid['6']/2., bin_cen['6'][-1]+bin_wid['6'][-1]/2.),
           '5':np.append(bin_cen['5'] - bin_wid['5']/2., bin_cen['5'][-1]+bin_wid['5'][-1]/2.),
           '4':np.append(bin_cen['4'] - bin_wid['4']/2., bin_cen['4'][-1]+bin_wid['4'][-1]/2.)
           }
Phi_obs = {'6':np.array([.0079, .242, .58, .9, 1.33, 4.6, 7.,  6.6, 8.3, 10.9, 23., 16.2]),
           '5':10.*np.array([.018,.018,.0092,.055,.083,.212,.277,.499,.628,.772,.639,
                             1.2,.6,1.8,2.98,7.78,5.39,10.7,12.5]),
           '4':10**np.array([-9.710,-9.534,-9.710,-9.057,-8.781,-8.534,-8.192,-7.966,-7.853,-7.642,
                             -7.552,-7.387,-7.077,-7.189,-6.899,-6.756,-6.745,-6.577,-6.487,-6.493,
                             -6.405,-6.341,-6.369,-6.297,-6.382,-6.298,-6.219,-6.088,-6.253]) * 1e9
           }
Phi_err = {'6':np.array([.0079, .061, .17, .32, .6, 1.2, 1.7, 2., 2.6, 3.6, 8.1, 16.2]),
           '5':10.*np.array([.013,.013,.0092,.023,.028,.044,.051,.068,.076,.087,.171,
                             1.58,1.39,1.75,2.01,2.81,2.46,3.2,3.4]),
           '4':np.array([0.138,0.169,0.138,0.292,0.402,0.534,0.791,1.026,1.169,1.490,1.657,2.545,
                         29.581,17.287,23.032,27.088,27.422,33.384,37.385,37.343,42.783,47.960,47.366,
                         51.538,46.855,52.208,59.367,71.892,80.733])
           }     


# !!!!!!!!!!!! eta_max should be used when 1/eta \simeq (1-eta)/eta breaks
eta_max = .5; eta_min = 0.057 # 0.057

def M1M0(M0,dt,f_duty,mu_fit,eta8,delta):
    M1 = 1e8*pow(mu_fit*f_duty*delta*dt/(eta8*10.*t_Edd)+pow(M0/1e8,delta),1./delta)

    ## eta mimicking Ueda14 empirical formula: eta = eta8*(M/M8)**delta
    eta = eta8*pow(M0/1e8,delta)
    ## M1: mass after growth following t^(1/delta) power
    M1 = 1e8*pow(mu_fit*f_duty*delta*dt/(eta8*10.*t_Edd)+pow(M0/1e8,delta),1./delta)

    ## if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_greater(eta, eta_max)
    M1[eta.mask] = M0[eta.mask]*np.exp(mu_fit*f_duty*dt/(eta_max*10.*t_Edd))
    # i = ma.argmax(eta)
    # M1[eta.mask] = M0[eta.mask]/M0[i]*M1[i]
    eta[eta.mask] = eta_max

    ## if eta < minimum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_less(eta, eta_min)
    i = ma.argmin(eta)
    M1[eta.mask] = M0[eta.mask]/M0[i]*M1[i]
    # the following exponential formula not continuous
    # !!!!!!!! M1[eta.mask] = M0[eta.mask]*np.exp(mu_fit*f_duty*dt/(eta_min*10.*t_Edd))    
    return M1


# ---------------- kernel*: kernel of calculating P(lbd) integral ----------------

# exponential growth
def kernel_MBH1(Mgrow_ratio, dt, f_duty, mu, sigma_dex):
    lbd = np.log(Mgrow_ratio)/( f_duty*dt/(0.1*10.*t_Edd) )
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# power law growth + exp extrapolation
def kernel_MBH2(M1, M0, dt, f_duty, mu, sigma_dex, eta8, delta):
    lbd = ( pow(M1/1e8, delta) - pow(M0/1e8, delta) )/(f_duty*delta*dt)*(eta8*10.*t_Edd)
    eta = eta8*pow(M0/1e8,delta)
    # if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_greater(eta, eta_max)
    # M1[eta.mask] = M0[eta.mask]*np.exp(mu*f_duty*dt/(eta_max*10.*t_Edd))
    i = ma.argmax(eta)
    lbd[eta.mask] = lbd[i] * np.log(M1/M0[eta.mask]) / np.log(M1/M0[i])
    # lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_max*10.*t_Edd))
    eta[eta.mask] = eta_max

    # if eta < minimum -> use Eddington accretion -- M(t) \propto M0
    eta = ma.masked_less(eta, eta_min)
    i = ma.argmin(eta)
    # lbd[i] \propto log(M1/M0[i]) from Eddington accretion
    lbd[eta.mask] = lbd[i] * np.log(M1/M0[eta.mask]) / np.log(M1/M0[i])
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# piece-wise lbd; 2 exp + 1 pow
def kernel_MBH3(M1, M0, dt, f_duty, mu, sigma_dex, eta8, delta):
    lbd = ( pow(M1/1e8, delta) - pow(M0/1e8, delta) )/(f_duty*delta*dt)*(eta8*10.*t_Edd)
    eta = eta8*pow(M0/1e8,delta)
    # if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_greater(eta, eta_max)
    # M1[eta.mask] = M0[eta.mask]*np.exp(mu*f_duty*dt/(eta_max*10.*t_Edd))
    lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_max*10.*t_Edd))
    eta[eta.mask] = eta_max
    # if eta < minimum -> use Eddington accretion -- M(t) \propto M0
    eta = ma.masked_less(eta, eta_min)
    lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_min*10.*t_Edd))
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# lambda from Schechter function 
def kernelS_MBH(Mgrow_ratio, dt, f_duty, l_cut):
    lbd = np.log(Mgrow_ratio)/( f_duty*dt/(0.1*10.*t_Edd) )
    return lbd/l_cut

def kernelS_M1450(M1450, MBH, l_cut):
    lbd = Lbol_M1450(M1450)/(1.25e38*MBH)
    return lbd/l_cut

def kernel_M1450(M1450, MBH, mu, sigma_dex):
    lbd = Lbol_M1450(M1450)/(1.25e38*MBH)
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# Willot+ 2010
def MF(M,z=6):
    alpha = -1.03
    Phi_star = 1.23e-8
    M_star = 2.24e9
    if z==6:
        return Phi_star*pow(M/M_star,alpha)*np.exp(-M/M_star)
    if z==4:
        M_star *= 10
        return Phi_star*pow(M/M_star,alpha)*np.exp(-M/M_star)

def L_M(M,Edd_ratio): # L_bol from M_BH in Msun
    return 1.25e38*Edd_ratio*M
def M_L(L,Edd_ratio): # M_BH in Msun from L_bol
    return L/(1.25e38*Edd_ratio)

def Mdot2M(Mdot):
    eta = 1
    beta = 2.775e-6*(1.5)**.5
    Mdot1 = 0.04
    Mdot2 = 0.1
    if Mdot<Mdot1:
        M = eta*Mdot/beta
    elif Mdot>Mdot2:
        M = (0.83*np.log10(Mdot)+2.48)*1.e5
    else:
        M1 = eta*Mdot1/beta
        M2 = (0.83*np.log10(Mdot2)+2.48)*1.e5
        t = (np.log(Mdot)-np.log(Mdot1))/(np.log(Mdot2)-np.log(Mdot1))
        M = np.exp( t*np.log(M2) + (1-t)*np.log(M1) )
    return M

def LF(l): # dn/dlogL in Mpc^-3 dex^-1
    Phi_M_star = 1.14e-8
    M_star = -25.13
    alpha  = -1.5; beta = -2.81
    Phi_L_star = Phi_M_star * 2.5
    L_star = pow(10,-.4*(M_star-34.1)) * 3e18/1450 *1e7 * fbol_1450
    L_1 = pow(10,-.4*(-27.2-34.1)) * 3e18/1450 *1e7 * fbol_1450 
    L_2 = pow(10,-.4*(-20.7-34.1)) * 3e18/1450 *1e7 * fbol_1450 
    # print('break L',L_star/W37, 'Phi_L_star', Phi_L_star)
    t = (np.log10(l) - np.log10(L_1)) / (np.log10(L_2) - np.log10(L_1))
    return Phi_L_star/( pow(l/L_star,-(alpha+1)) + pow(l/L_star,-(beta+1)) ) * (2*(1-t)+3*t) 

def LF_M1450(M,z=6): # dn/dmag in Mpc^-3 mag^-1
    if z==6: 
        # Willot 2010 CFHQS + SDSS
        Phi_M_star = 1.14e-8
        M_star = -25.13
        alpha  = -1.5; beta = -2.81
        # Matsuoka 2018
        Phi_M_star = 1.09e-8
        M_star = -24.9
        alpha  = -1.23; beta = -2.73
    elif z==5:
        # McGreer 2018 DPL; 
        Phi_M_star = pow(10., -8.97+0.47)
        M_star = -27.47
        alpha  = -1.97; beta = -4.
        # refit by Matsuoka 2018 (beta & M_star); me: (alpha & Phi_M_star)
        Phi_M_star = 3.8e-8
        M_star = -25.6
        alpha  = -1.23; beta = -3.
        # Niida 2020 Table 4
        Phi_M_star = 1.01e-7
        M_star = -25.05
        alpha  = -1.22; beta = -2.9
    elif z==4: # Akiyama 2018
        Phi_M_star = 2.66e-7
        M_star = -25.36
        alpha  = -1.3; beta = -3.11
    else:
        print("wrong redshift")
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) #* (2*(1-t)+3*t) 

def M1450_Lbol(L):
    return 34.1-2.5*np.log10(L/(fbol_1450*3e18/1450*1e7))
def Lbol_M1450(M):
    return pow(10., -0.4*(M-34.1)) * (fbol_1450*3e18/1450*1e7)

# X-ray bolometric correction; Hopkins+07 & Duras+20
def K_AVE07(Lbol):
    return 10.83*pow(Lbol/(1e10*Lsun),0.28)+6.08*pow(Lbol/(1e10*Lsun),-0.02)
def K_AVE20(Lbol):
    a = 10.96
    b = 11.93
    c = 17.79
    return a*( 1 + pow(np.log10(Lbol/Lsun)/b,c) )

def KX_AVE20(Lx):
    a = 15.33
    b = 11.48
    c = 16.2
    return a*( 1 + pow(np.log10(Lx/Lsun)/b,c) )

# obscured fraction = Type II AGN fraction
def f_obsc_U14(logLx,z): # Ueda 14; 22< log NH < 24 fraction; as a func of Lx
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+z)**a1 if z<=2. else phi4375_0*(1+2.)**a1
    phimax = .84
    phimin = .2
    beta = .24
    phi = min( phimax, max(phi4375_z - beta*(logLx-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    return f_obsc_sum

# constant obscured fraction; motivated by Vito+ 2018
f_obsc_const = .8

# correction factor including Compton thick AGNs; different fbol_Xray
def corr_U14H07(M1450): # Ueda+14 & Shankar+09
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE07(L_bol)
    Lx = L_bol/f_bol
    eta = 1.7
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+2.)**a1
    phimax = .84
    phimin = .2
    beta = .24
    phi = min( phimax, max(phi4375_z - beta*(np.log10(Lx)-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    f_CTK = phi/2.
    return 1./(1-f_obsc_sum)
def corr_U14D20(M1450): # Ueda 14
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE20(L_bol)
    Lx = L_bol/f_bol
    eta = 1.7
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+2.)**a1
    phimax = .84 #(1+eta)/(3+eta)
    phimin = .2
    beta = .24
    if isinstance(M1450,float):
        phi = min( phimax, max(phi4375_z - beta*(np.log10(Lx)-43.75), phimin))
    else:
        phi = np.zeros(len(M1450))
        for i in range(len(M1450)):
            phi[i] = min( phimax, max(phi4375_z - beta*(np.log10(Lx[i])-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    f_CTK = phi/2.
    return 1./(1-f_obsc_sum)

def LF_M1450_CO(M,z): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    return LF_M1450(M,z)/(1-f_obsc_const)

def LF_M1450_DO(M,z): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    return LF_M1450(M,z)*corr_U14D20(M)

def t_freefall(nH):
    C = np.sqrt( 32*G*(mu*m_H)/ (3*pi) )
    return 1./C/np.sqrt(nH)

def t_from_z(z): # age of universe at redshift z: tH = 2/(3Hz)
    return 2./(3*H0*np.sqrt(Omega_m0)) * pow(1+z, -1.5)

def Tv(Mh,z):
    return alpha_T * (Mh/1.e8)**(2./3.) *  (1+z)/11.

def Mh_Tv(Tv,z):
    return 1.e8*(Tv/alpha_T/(1+z)*11.)**1.5

def Omega_mz(z):
    return Omega_m0*(1+z)**3 /(Omega_m0*(1+z)**3 + Omega_L0)

def Hz(z):
    return H0*np.sqrt( Omega_m0*(1+z)**3 + Omega_L0 ) 
 
def RHO_crit(z):
    return 3*pow(H0,2)/(8*pi*G)*(1+z)**3*Omega_m0/Omega_mz(z) 

class HALO:
    def __init__(self,M,z0):
        self.Mh = M
        self.z = z0
        self.c = 18*pow(self.Mh/(1.e11*Ms), -0.13)/(1+self.z) #concentration parameter c from Dekel & Birnboim 2006 Eq(22)
        c, z = self.c, self.z
        self.d = Omega_mz(z) - 1 
        d = self.d
        self.Delta_crit = 18.0*pi*pi + 82*d - 39*d*d  # Delta_crit ~ 200, overdensity
        Delta_crit = self.Delta_crit

        self.delta0 = self.Delta_crit/3.*pow(c,3)/(-c/(1+c) + np.log(1+c)) # characteristic overdensity parameter 
        delta0 = self.delta0
        self.rho_crit = RHO_crit(z)  # mean density of DM at z
        self.rho_c = self.rho_crit * delta0 

        self.Rvir = pow( self.Mh/(4./3*pi*Delta_crit*self.rho_crit),1./3. ) 
        self.Rs = self.Rvir/self.c 
        self.Vc = np.sqrt(G*self.Mh/self.Rvir) 

        self.t_dyn = self.Rvir/self.Vc 
        self.Tvir = G*self.Mh*(mu*m_H)/(2.*k_B*self.Rvir)
        self.gc = 2*c/(np.log(1+c) - c/(1+c)) 
        self.alpha = self.Tvir/self.Mh**(2./3)

    def Rho_r(self, r):
        rho_crit, delta0, Rvir = self.rho_crit, self.delta0, self.Rvir
        c, x = self.c, r/Rvir
        return rho_crit*delta0/( c*x * (1+c*x)**2 )

    # x = r/Rvir  c = Rvir/Rs
    def F_NFW(self,x):
        c = self.c
        return -c*x/(1+c*x) + np.log(1+c*x)
        
    def M_enc(self,r):
        rho_crit, delta0, Rs, Rvir = self.rho_crit, self.delta0, self.Rs, self.Rvir
        M_r = 4*pi*rho_crit*delta0*pow(Rs,3)*self.F_NFW(r/Rvir)
        return M_r

    def Phi(self, r):
        # lim r -> 0
        #return -4*pi*G*rho_crit*delta0*Rs*Rs 
        rho_crit, delta0, Rs = self.rho_crit,  self.delta0, self.Rs
        return -4*pi*G*rho_crit*delta0*(Rs**3)/r*np.log(1+r/Rs) 
