from PYmodule import *

print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
print('z:32 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(32))/Myr)
print('t_Edd         : %3.2f Myr', t_Edd/Myr)

NM = 3

T_tell = 8000
eta = 0.3


T = ascii.read('../z6/data/MF')
T  = T[np.logical_and(T['bin_cen']>1e6,T['bin_cen']<1e11)]
dlog10M = np.log10(T['bin_cen'][1]/T['bin_cen'][0])

# LF test bins
bin_obs = np.linspace(-29,-22,50)
bin_wid = bin_obs[1:]-bin_obs[:-1]
bin_cen = bin_obs[:-1]+bin_wid/2.
Nbin = len(bin_cen)

for M1450 in [-29., -22.]:
    print('M1450',M1450,'Mbh = %3.2e'%M_L(Lbol_M1450(M1450),1.))

flambda = .19 # .18 .20
z = 6

f_duty = .5
sigma_fit = .12
mu_fit = flambda/f_duty
Phi = np.zeros(Nbin)

for ibin in range(Nbin): # Nbin
    for i in range(len(T)): # loop over each M_BH at z 
        M_BH = T['bin_cen'][i]
        dn_MBH = T['Phi'][i]*dlog10M
        kernel = kernel_M1450(bin_obs[ibin], M_BH, mu_fit, sigma_fit)

        if -10.<kernel<10.:
            x0 = kernel_M1450(bin_obs[ibin+1], M_BH, mu_fit, sigma_fit)
            x1 = kernel_M1450(bin_obs[ibin], M_BH, mu_fit, sigma_fit)
            dP_M1450 = .5*(math.erfc(x0) - math.erfc(x1))
            dP = dn_MBH*dP_M1450
            Phi[ibin] += dP/bin_wid[ibin]*f_duty

T = Table(
    [bin_cen, Phi, LF_M1450_DO(bin_cen,z)],
    names=('bin_cen','Phi','Phi_DO')
)
# ascii.write(T, z6datapre+'allmassMF_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm'+str(Nbsm-1)+'alpha1'+'NM'+str(NM),formats={'bin_left':'6.2e','Phi':'4.2e'},overwrite=True)
ascii.write(T, z6datapre+'LF11',formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_DO':'4.2e'},overwrite=True)