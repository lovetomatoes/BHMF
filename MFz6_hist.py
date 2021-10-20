from PYmodule import *

print('z:45, t_Hubble: ', t_from_z(45)/Myr)
print('z:30, t_Hubble: ', t_from_z(30)/Myr)
print('z:10, t_Hubble: ', t_from_z(10)/Myr)
print('z:5,  t_Hubble: ', t_from_z(5)/Myr)

iM = 0

T_tell = 8000
eta = 0.3

Ts = [] # bsm=0,1 two files
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
    Ts.append(T)

# MF
abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
flambda = .19 # .18 .20

N_concatenate = int(1e3)
z = 6

for f_duty in [0.5]: # .6 .4 
    for sigma_fit in [.12]: # .10  .14
        mu_fit = flambda/f_duty
        Phi = np.zeros(len(abin_mf)-1)
        for i_concatenate in range(N_concatenate):
            for i_bsm in range(Nbsm):
                T = Ts[i_bsm]
                kernel = norm.rvs(size=len(T))
                T['Edd_ratio'] = mu_fit * np.exp(kernel*sigma_fit*np.log(10))
                T['Mstar_z'] = T['Mstar0']* np.exp( (t_from_z(z)-t_from_z(T['z_col'])) / t_Edd * f_duty* T['Edd_ratio'])
                # MF
                hist_mf, bin_edges = np.histogram(T['Mstar_z'],bins=abin_mf,density=False)
                Phi += hist_mf*n_base[iM]*f_bsm[i_bsm]/(1e4*N_concatenate)/dlog10M
        T = Table(
            [abin_mf[:-1], Phi],
            names=('bin_left','Phi')
        )
        # ascii.write(T, z6datapre+'normMF_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm'+str(Nbsm-1)+'alpha1'+'N'+str(int(np.log10(N_concatenate))),formats={'bin_left':'%6.2f','Phi':'4.2e'},overwrite=True)
        ascii.write(T, z6datapre+'histMF',formats={'bin_left':'%6.2f','Phi':'4.2e'},overwrite=True)