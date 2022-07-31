from PYmodule import *
from PYmodule.models_logd import *
# tau=10, 18(f2/ best fit)

f_seed = 0.01
ts = [10, 18, 50]

f_seed = 0.1
ts = [10, 20, 50]

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))

readprex = z6datapre+f_seedlabel+'Phi_easy'
prex = '../' + f_seedlabel

labels = [r'$\mathrm{t_{life}}$', r'$\log \delta$', r'$\lambda_0$', r'$\alpha$']


colors = ['c','purple','y']
figm = plt.figure(num='MF',figsize=(10, 10))
axm = figm.add_subplot(1, 1, 1)
figl = plt.figure(num='LF',figsize=(10, 10))
axl = figl.add_subplot(1, 1, 1)

for i in range(len(ts)):
    t = ts[i]
    MFname = readprex +'MFt{:d}'.format(int(t))
    LFname = readprex +'LFt{:d}'.format(int(t))
    TM = ascii.read(MFname, guess=False, delimiter=' ')
    TL = ascii.read(LFname, guess=False, delimiter=' ')
    axm.plot(M_BH, TM['Phi'],c=colors[i], label='t=%d Myr'%t)
    axm.plot(M_BH, TM['W10_MF'],c='C0')

    axl.plot(M1450, TL['Phi_DO'],c=colors[i], label='t=%d Myr'%t)
    axl.errorbar(bin_cen, Phi_obs, yerr=Phi_err,color='C0',fmt='o',capsize=10)



axm.set_xlim(1e7,1e10); axm.set_xscale('log')
axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
axm.tick_params(labelsize=fstick)
figm.savefig('../%sMFtau.png'%f_seedlabel)

axl.set_xlim(np.max(M1450),np.min(M1450))
axl.set_ylim(5e-3,1e2)
axl.set_yscale('log')
# axl.set_xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
# axl.set_ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
axl.tick_params(labelsize=fstick)
figl.savefig('../%sLFtau.png'%f_seedlabel)