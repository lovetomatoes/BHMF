from PYmodule import *

Ts={}
alpha = 1.
f_duty = {}

z = int(6)
f_duty[str(z)] = .7; mu_fit = .21; sigma_fit = .15; eta8 = 0.1; delta_fit = 0.001
fname=z6datapre+'MF2e10_'+'f%3.2f'%f_duty[str(z)]+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

z = int(4)
f_duty[str(z)]= 1; mu_fit = .1; sigma_fit = .1; eta8 = 0.05; delta_fit = 0.1 # best
f_duty[str(z)]= .5; mu_fit = .2; sigma_fit = .05; eta8 = 0.05; delta_fit = 0.1 # best

fname=z4datapre+'MF_'+'z%d'%z+'f%3.2f'%f_duty[str(z)]+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

plt.figure(figsize=(10,10),dpi=400)

for z in [6, 4]:
    T = Ts[str(int(z))]
    M_BH = T['M_BH']
    PhiM = T['dn_dlog10M']*f_duty[str(z)]
    plt.plot(np.log10(M_BH),PhiM,label='z'+str(z))


plt.xlim(7,10.5)
plt.yscale('log'); plt.ylim(2e-9,1e-4)
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.savefig('../Phi_M.png')