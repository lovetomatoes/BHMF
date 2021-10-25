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

Ts = [] # bsm=0,1 two files
for iM in range(NM):
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

Nbin = len(abin_mf)-1

flambda = .19 # .18 .20
z = 6
N = 10000

for f_duty in [0.5]: # .6 .4 
    for sigma_fit in [.12]: # .10  .14
        mu_fit = flambda/f_duty
        Phi = np.zeros(Nbin)

        for ibin in range(Nbin): # Nbin
            for iM in range(NM):
                for i_bsm in range(Nbsm):
                    T = Ts[iM][i_bsm]
                    dP = 0
                    for i in range(N):  
                        x0 = kernel_MBH(abin_mf[ibin]/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit)
                        x1 = kernel_MBH(abin_mf[ibin+1]/T['Mstar0'][i],t_from_z(z)-t_from_z(T['z_col'][i]),f_duty, mu_fit, sigma_fit)
                        dP_MBH = .5*(math.erfc(x0) - math.erfc(x1))
                        dP += dP_MBH/N
                    Phi[ibin] += dP*n_base[iM]*f_bsm[i_bsm]/dlog10M
        T = Table(
            [abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0]), Phi],
            names=('bin_cen','Phi')
        )
        # ascii.write(T, z6datapre+'allmassMF_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm'+str(Nbsm-1)+'alpha1'+'NM'+str(NM),formats={'bin_left':'6.2e','Phi':'4.2e'},overwrite=True)
        ascii.write(T, z6datapre+'MF',formats={'bin_cen':'6.2e','Phi':'4.2e'},overwrite=True)