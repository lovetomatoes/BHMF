from PYmodule import *

t1 = time.time()

f_seed = 0.01
fMname = '../rhoM_evol.txt'# %int(abs(np.log10(f_seed)))
fLname = '../rhoL_evol.txt'
prex = '../ndraw100'
fMname = prex+'rhoM_evol.txt'# %int(abs(np.log10(f_seed)))
fLname = prex+'rhoL_evol.txt'

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
prex = '../4p/M0r8_' + f_seedlabel
ndraw = 3
fMname = prex+'ndraw{:d}rhoM_evol.txt'.format(ndraw)
fLname = prex+'ndraw{:d}nL_evol.txt'.format(ndraw)

TM = ascii.read(fMname, guess=False, delimiter=' ')
TL = ascii.read(fLname, guess=False, delimiter=' ')

zs = TM['zs']
rhoM = TM['rhoM']; rhoM_med = TM['rhoM_med']; rhoM_std = TM['rhoM_std']
nL = TL['nL']; nL_med = TL['nL_med']; nL_std = TL['nL_std']

figl, axl = plt.subplots(num='LF',figsize=(10, 10))
z_ = np.arange(5.6,7.5,.1)
dash_k = 0.39*pow(10., -.78*(z_-6.7))
axl.plot(z_,dash_k, '--',c='black')
axl.errorbar(6.,1.33,yerr=0.33,fmt='D',elinewidth=2,capsize=4,label='Jiang+ 2016',
            c='red',ms=10,mfc='none')
axl.errorbar(6.7,0.39,yerr=0.11,fmt='s',elinewidth=2,capsize=4,label='Wang+ 2019',
            c='red',ms=10,mfc='none')

for i in range(len(zs)):
    z = zs[i]
    axl.errorbar(z,nL[i],yerr=nL_std[i],fmt='o',ms=10,elinewidth=2,capsize=4,
            c='C%d'%z,label='z=%d'%z,mfc='none')
axl.plot(zs,nL,c='grey')
axl.set_xlim(5.5,10.5)
# axl.set_ylim(1e-10,1e-4); 
axl.set_yscale('log')
axl.legend(fontsize=fslegend)
# axl.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axl.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axl.tick_params(labelsize=fslabel)
figl.savefig(prex+'nL.png',dpi=300,bbox_inches='tight')


figm, axm = plt.subplots(num='MF',figsize=(10, 10))
for i in range(len(zs)):
    z = zs[i]
    axm.errorbar(z, rhoM[i],yerr=rhoM_std[i],fmt='o',elinewidth=2,capsize=4,
            color='C%d'%z,label='z=%d'%z)
axl.plot(zs,rhoM,c='grey')
axm.set_xlim(5.5,10.5)
# axm.set_ylim(1e-10,1e-4); 
axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fslabel)
figm.savefig(prex+'rhoM.png',dpi=300,bbox_inches='tight')
