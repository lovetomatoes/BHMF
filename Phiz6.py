from PYmodule import *

print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
print('z:32 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(32))/Myr)
print('t_Edd         : %3.2f Myr', t_Edd/Myr)

N_Mh = 1 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 1000
eta = 0.3
z = 6


Ts = [] # bsm=0,1 two files
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        T['Mstar_z'] = np.zeros(len(T))
        T['Lbol_z'] = np.zeros(len(T))
        T['M1450_z'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['type'] = np.ones(len(T))
        T['Edd_ratio'] = np.ones(len(T))
        TM.append(T)
    Ts.append(TM)

# MF
abin_mf =  np.logspace(6,12,num=100) # default endpoint=True
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# LF test bins
bin_obs = np.linspace(-29,-22,50)
bin_wid = bin_obs[1:]-bin_obs[:-1]
bin_cen = bin_obs[:-1]+bin_wid/2.
N_lf = len(bin_cen)

for M1450 in [-29., -22.]:
    print('M1450',M1450,'Mbh = %3.2e'%M_L(Lbol_M1450(M1450),1.))



for flambda in [.19]: # .18 .20
    for f_duty in [0.5]: # .6 .4 
        for sigma_fit in [.12]: # .10  .14
            mu_fit = flambda/f_duty

            dn_MBH = np.zeros(N_mf)
            for ibin in range(N_mf): # N_mf
                for iM in range(N_Mh):
                    for i_bsm in range(Nbsm):
                        T = Ts[iM][i_bsm]
                        dP = 0
                        for i in range(Ntr):  
                            x0 = kernel_MBH(abin_mf[ibin]/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit)
                            x1 = kernel_MBH(abin_mf[ibin+1]/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit)
                            dP_MBH = .5*(math.erfc(x0) - math.erfc(x1))
                            dP += dP_MBH/Ntr
                        dn_MBH[ibin] += dP*n_base[iM]*f_bsm[i_bsm]
            T = Table(
                [abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0]), dn_MBH],
                names=('M_BH','dn_MBH')
            )

            T  = T[np.logical_and(T['M_BH']>1e6,T['M_BH']<1e11)] # select M_BH range
            
            Phi = np.zeros(N_lf)
            for ibin in range(N_lf): # N_lf
                for i in range(len(T)): # loop over each M_BH at z 
                    M_BH = T['M_BH'][i]
                    kernel = kernel_M1450(bin_obs[ibin], M_BH, mu_fit, sigma_fit)
                    if -10.<kernel<10.:
                        x0 = kernel_M1450(bin_obs[ibin+1], M_BH, mu_fit, sigma_fit)
                        x1 = kernel_M1450(bin_obs[ibin], M_BH, mu_fit, sigma_fit)
                        dP_M1450 = .5*(math.erfc(x0) - math.erfc(x1))
                        dP = T['dn_MBH'][i]*dP_M1450
                        Phi[ibin] += dP/bin_wid[ibin]*f_duty

            T = Table(
                [bin_cen, Phi, LF_M1450_CO(bin_cen,z), LF_M1450_DO(bin_cen,z)],
                names=('bin_cen','Phi','Phi_CO','Phi_DO')
            )
            # ascii.write(T, z6datapre+'allmassMF_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm'+str(Nbsm-1)+'alpha1'+'N_Mh'+str(N_Mh),formats={'bin_left':'6.2e','Phi':'4.2e'},overwrite=True)
            ascii.write(T, z6datapre+'LF11_',formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_CO':'4.2e','Phi_DO':'4.2e'},overwrite=True)