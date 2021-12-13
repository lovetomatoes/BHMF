from PYmodule import *

Ts={}
alpha = 1.
f_duty = {}

z = int(6)
f_duty[str(z)] = .7; mu_fit = .21; sigma_fit = .15; eta8 = 0.1; delta_fit = 0.001
fname=z6datapre+'MF2e10_'+'f%3.2f'%f_duty[str(z)]+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')

plt.figure(figsize=(10,10),dpi=400)
T = Ts[str(int(z))]
M_BH = T['M_BH']
PhiM = T['dn_dlog10M']*f_duty[str(z)]
plt.plot(np.log10(M_BH),PhiM,label='z'+str(z))

z = int(4)
# f_duty[str(z)]= 1; mu_fit = .1; sigma_fit = .1; eta8 = 0.05; delta_fit = 0.1 # best
# f_duty[str(z)]= .5; mu_fit = .2; sigma_fit = .05; eta8 = 0.05; delta_fit = 0.1 # best
# f_duty[str(z)]= 10.; mu_fit = .1; sigma_fit = .1; eta8 = 0.05; delta_fit = 0.1 # best

d_range = [.1]
e_range = [.05]
f_range = [1, .5, 10., 100., 1e3]
m_range =[.1, .2, .01, .001, 1e-4]
s_range = [.1]

f_range = [1, 1e3]
m_range =[.1, 1e-4]

i = 0
for f_duty[str(z)] in f_range:
    for mu_fit in m_range:
        for sigma_fit in s_range:
            for eta8 in e_range:
                for delta_fit in d_range:
                    fname=z4datapre+'MF_'+'z%d'%z+'f%3.2f'%f_duty[str(z)]+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    if os.path.isfile(fname):
                        # i += 1; continue
                        Ts[str(z)] = ascii.read(fname, guess=False,delimiter=' ')
                    else:
                        # print('nofile',fname)
                        continue
                    i += 1 
                    T = Ts[str(int(z))]
                    M_BH = T['M_BH']
                    PhiM = T['dn_dlog10M']#*f_duty[str(z)]
                    plt.plot(np.log10(M_BH),PhiM,label='f=%.1e'%f_duty[str(z)]+'; m=%.1e'%mu_fit)

print('z=4 i=%d files plotted'%i)

plt.xlim(2,15); plt.ylim(1e-14,1e0)
plt.xlim(7,10.5); plt.ylim(2e-9,1e-4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.savefig('../Phi_M.png')