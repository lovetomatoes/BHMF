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

ascii.write(Table([t_from_z(T['z_col'])/Myr, T['Mh_col']], names =['t_Myr','M_ini']), '../M_ini_Mh1e11bsm0alpha1',
formats={'t_Myr':'4.2e','M_ini':'4.2e'}, overwrite=True)