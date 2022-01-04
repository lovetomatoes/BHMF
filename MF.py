from PYmodule import *
# several MF by best fitting parameters
# models: 'Schechter best fit','Log-normal best fit','Schechter + mass dependent ER'

Ts={}
alpha = 1.
f_duty = {}

z = int(6)
f_duty[str(z)] = .7; mu_fit = .21; sigma_fit = .15; eta8 = 0.1; delta_fit = 0.001
fname=z6datapre+'MF2e10_'+'f%3.2f'%f_duty[str(z)]+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

plt.figure(figsize=(10,10),dpi=400)
fname = z6datapre + 'MF_SC_t8.0e+02f0.8l1.0e-01a0.500alpha1.0'
T = ascii.read(fname, guess=False,delimiter=' ')
M_BH = T['M_BH']
PhiM = T['dn_dlog10M']
plt.plot(np.log10(M_BH),MF(M_BH),c='black',label='Willot10')
plt.plot(np.log10(M_BH),PhiM,label='Schechter best fit')

fname = z6datapre + 'MF_LN_t9.0e+02f0.8m5.0e-02s0.300alpha1.0'
T = ascii.read(fname, guess=False,delimiter=' ')
M_BH = T['M_BH']
PhiM = T['dn_dlog10M']
plt.plot(np.log10(M_BH),PhiM,label='Log-normal best fit')

fname = z6datapre + 'MF_SC_t4.0e+02f0.6d1.0l8.0e-01a0.100alpha1.0'
T = ascii.read(fname, guess=False,delimiter=' ')
M_BH = T['M_BH']
PhiM = T['dn_dlog10M']
plt.plot(np.log10(M_BH),PhiM,label='Schechter + mass dependent ER')

plt.xlim(2,15); plt.ylim(1e-14,1e0)
plt.xlim(7,10.5); plt.ylim(2e-9,1e-4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right',fontsize=fslabel)
plt.xlabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel);
plt.ylabel(r'$\mathrm{Mpc^{-3} dex^{-1}}$',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig('../Phi_M.png')