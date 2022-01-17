from PYmodule import *

z = int(6)
tz = t_from_z(z)

# print(M1M0(np.array([1e2,1e5]),tz-t_from_z(20),1,.1,.1,0.001))
# print(M1M0(np.array([1e2,1e5]),50*Myr,1,.1,.1,0.001))

# f_duty=.7;mu_fit=.21*30;sigma_fit=.15
# M0 = 1e4; dt = 500*Myr
# print('M1=%.1e'%(M0*np.exp(mu_fit*f_duty*dt/t_Edd)), 'N_mf=%d'%N_mf)

# MF
abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
# bin_edg = np.linspace(-29, 10)
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_edg) - 1

f_duty = 1.
d_fit = .5
l_cut = .5
a = .1

dn_MBH = np.zeros(N_mf)
Mstar0 = 100
dt_seed = 500*Myr

T_k = ascii.read(z6datapre+'fort.16')
hist, bin_edges = np.histogram(pow(10., T_k['logMBH']),bins=abin_mf,density=False)
# print('hist len:',np.sum(hist))
hist = hist/len(T_k)*1e-3
Phi_k = hist/dlog10M

# to be consistent w/ Phi_k; change t_Edd to 400Myr in PYmodule
Nt = 2 # int
dt_seed /= Nt

for i in range(Nt):
    for ibin in range(N_mf):
        if i==0:
            # #----------- Schechter lbd -----------
            x0 = kernelS_MBH_M(bin_left[ibin],  Mstar0, dt_seed, f_duty, l_cut, d_fit)
            x1 = kernelS_MBH_M(bin_right[ibin], Mstar0, dt_seed, f_duty, l_cut, d_fit)
            x0 = max(x0,0.); x1 = max(x1,0.)
            # x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
            dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)
            dn_MBH[ibin] = dP_seed*1e-3 # Mpc^-3
        else:
            dn_MBH_prev = dn_MBH.copy()
            x0 = kernelS_MBH_M(M_BH[ibin], bin_right, dt_seed, f_duty, l_cut, d_fit)
            x1 = kernelS_MBH_M(M_BH[ibin], bin_left,  dt_seed, f_duty, l_cut, d_fit)
            x0[x0<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
            dn_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dn_MBH_prev)

consv_ratio = np.nansum(dn_MBH)/1e-3
print('mass consv_ratio',consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH), Phi_k],
    names=('M_BH','Phi','W10_MF','Phi_k')
)
MFname = z6datapre+'MF100'
# MFname = z6datapre+'MF100'+ \
#                     'f%.2f'%f_duty+ \
#                     'd%.2f'%d_fit+ \
#                     'l%.2f'%l_cut+ \
#                     'a%.2f'%a
ascii.write( Table([np.log10(T['M_BH']), T['Phi'], T['W10_MF'], T['Phi_k']],
            names=['M_BH','Phi','W10_MF','Phi_k']),
            MFname,
            formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e','Phi_k':'4.2e'},
            overwrite=True)


# # --------- Luminosity Function ---------
Phi = np.zeros(N_lf)
Phi_csv = 0.
T['dn_MBH'] = T['Phi']*dlog10M
for ibin in range(N_lf):
    # #----------- Schechter lbd -----------
    x0 = kernelS_M1450(bin_edg[ibin+1], T['M_BH'], l_cut)
    x1 = kernelS_M1450(bin_edg[ibin],   T['M_BH'], l_cut)
    dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

    dPhi = np.nansum(T['dn_MBH']*dP_M1450)
    Phi_csv += dPhi
    Phi[ibin] = dPhi/bin_wid[ibin]*f_duty
# print('consv of dP_M1450:',Phi_csv/1e-3)

Phi *= 1e9
Phi_DO = Phi/corr_U14D20(bin_cen)
Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)

T = Table(
    [bin_cen,Phi_obs,Phi_DO,Phi,Chi2*np.ones(N_lf)],
    names=('bin_cen','Phi_obs','Phi_DO','Phi','Chi2')
)
LFname = z6datapre+'LF100'
ascii.write(T, LFname,
            formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            overwrite=True)