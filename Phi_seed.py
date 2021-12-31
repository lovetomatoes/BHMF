from PYmodule import *

print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
print('z:32 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(32))/Myr)
print('t_Edd         : %3.2f Myr', t_Edd/Myr)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3


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
iM = 0; i_bsm = 0
T = ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
T['Mstar0'] = np.zeros(len(T))
for i in range(len(T)):
    T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
T['dt'] = (t_from_z(6.) - t_from_z(T['z_col']))/Myr

ascii.write(Table([t_from_z(T['z_col'])/Myr, T['Mstar0']], names =['t_Myr','M_ini']), '../M_ini_Mh1e11bsm0alpha1',
formats={'t_Myr':'4.2e','M_ini':'4.2e'}, overwrite=True)
print(T.info)
print('mean growth time to z=6',np.mean(T['dt']),'mean Mstar0',np.mean(T['Mstar0']))

abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

hist_mf, bin_edges = np.histogram(T['Mstar0'],bins=abin_mf,density=False)
dn_MBH =  hist_mf/len(T) * n_base[0]

ascii.write( Table([M_BH, dn_MBH/dlog10M],names=['M_BH','dn_dlog10M']),
            z6datapre+'IMF',
            formats={'M_BH':'4.2e','dn_dlog10M':'4.2e'},
            overwrite=True)