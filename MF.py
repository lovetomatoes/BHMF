from PYmodule import *

Ts={}
alpha = 1.

z = int(6)
f_duty = .7; mu_fit = .21; sigma_fit = .15; eta8 = 0.1; delta_fit = 0.001
fname=z6datapre+'MF2e10_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
print(fname)
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

z = int(4)
f_duty = .6; mu_fit = .1; sigma_fit = .29; eta8 = 0.07; delta_fit = 0.25 # best
# f_duty = .4; mu_fit = .2; sigma_fit = .24; eta8 = 0.08; delta_fit = 0.2 # for mu closer to .3
fname=z4datapre+'MF_'+'z%d'%z+'z6_2e10f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
print(fname)
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

plt.figure(figsize=(10,10),dpi=400)

for z in [6, 4]:
    T = Ts[str(int(z))]
    M_BH = T['M_BH']
    if z==4:
        PhiM = T['dn_dlog10M']/2.
    else:
        PhiM = T['dn_dlog10M']
    plt.plot(np.log10(M_BH),PhiM,label='z'+str(z))
    # if z==4:
    #     plt.plot(np.log10(M_BH/3.),PhiM,label='z'+str(z)+' M/3')

z = int(4)
f_duty = .6; mu_fit = .1; sigma_fit = .29; eta8 = 0.07; delta_fit = 0.25 # best
# f_duty = .4; mu_fit = .2; sigma_fit = .24; eta8 = 0.08; delta_fit = 0.2 # for mu closer to .3
fname=z4datapre+'MF_'+'z%d'%z+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
print(fname)
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

for z in [4]:
    T = Ts[str(int(z))]
    M_BH = T['M_BH']
    PhiM = T['dn_dlog10M']/2.
    plt.plot(np.log10(M_BH),PhiM,label='z'+str(z))
    # if z==4:
    #     plt.plot(np.log10(M_BH/3.),PhiM,label='z'+str(z)+' M/3')

plt.xlim(7,10.5)
plt.yscale('log'); plt.ylim(2e-9,1e-4)
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.savefig('../Phi_M.png')