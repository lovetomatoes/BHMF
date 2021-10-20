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
    T['Edd_ratio_now'] = np.ones(len(T))
    Ts.append(T)

# luminosity bin same as Matsuoka18 setting
bin_cen = [-29,-27.5,-26.75,-26.25,-25.75,-25.25,-24.75,-24.25,-23.75,-23.25,-22.75,-22]
bin_wid = [2, 1, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1]
bin_obs = np.zeros(len(bin_cen)+1)
for i in range(len(bin_cen)):
    bin_obs[i] = bin_cen[i] - bin_wid[i]/2.
bin_obs[-1] =  bin_cen[i] + bin_wid[i]/2.

# test bins
bin_obs = np.linspace(-29,-19,50)
bin_wid = bin_obs[1:]-bin_obs[:-1]
bin_cen = bin_obs[:-1]+bin_wid/2.


flambda = .19 # .18 .20

N_concatenate = int(1e3)
z = 6

for f_duty in [0.5]: # .6 .4 
    for sigma_fit in [.12]: # .10  .14
        mu_fit = flambda/f_duty
        Phi = np.zeros(len(bin_cen))
        for i_concatenate in range(N_concatenate):
            for i_bsm in range(Nbsm):
                T = Ts[i_bsm]
                kernel = norm.rvs(size=len(T))
                T['Edd_ratio'] = mu_fit * np.exp(kernel*sigma_fit*np.log(10))
                T['Mstar_z'] = T['Mstar0']* np.exp( (t_from_z(z)-t_from_z(T['z_col'])) / t_Edd * f_duty* T['Edd_ratio'])
                T['Edd_ratio_now'] = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=len(T)) # 使用rvs随机生成lbd 未用本身的Edd_ratio
                T['Lbol_z'] = L_M(T['Mstar_z'], T['Edd_ratio_now'])
                T['M1450_z'] = M1450_Lbol(T['Lbol_z'])
            # LF
                hist_lf, bin_edges = np.histogram(T['M1450_z'],bins=bin_obs,density=False)
                Phi += hist_lf*n_base[iM]*f_bsm[i_bsm]/(1e4*N_concatenate)/bin_wid*f_duty
        T = Table(
            [bin_cen, Phi, LF_M1450(bin_cen,z)],
            names=('bin_cen','Phi','Phi_matsu')
        )
        ascii.write(T, z6datapre+'histLF',formats={'bin_cen':'%6.2f','Phi':'4.2e','Phi_matsu':'4.2e'},overwrite=True)