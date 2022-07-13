from PYmodule import *
from PYmodule.models_logd import *

t1 = time.time()

corr = 'U'
f_seed = 0.01
z = 6
zmax = 11
M_up  = 1e7
L_up  = -26
# L_ups = [-23,-24,-25,-26,-27] # 可用集合
ndraw = 3

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
readprex = '../4p/M0r8_' + f_seedlabel
readprex = '../M0r8_' + f_seedlabel
fname = readprex + '.h5'
prex = readprex + 'ln'
# fname = '../4p/3n2f_logd_4pr8_f2.h5'

labels = [r'$\mathrm{t_{life}}$', r'$\log \delta$', r'$\lambda_0$', r'$\alpha$']

reader = emcee.backends.HDFBackend(fname)
ndim = len(labels)

tau = reader.get_autocorr_time(tol=1)
tau = np.max(tau); print('max tau:',tau)
Nburnin = int(3*tau)
Nthin = int(tau/2)

# Nburnin = 0
Nthin = 1

samples = reader.get_chain(discard=Nburnin)
probs = reader.get_log_prob(discard=Nburnin)
print('len of samples:', len(samples))

samples = reader.get_chain(discard=Nburnin, thin=Nthin, flat=True)
probs = reader.get_log_prob(discard=Nburnin, thin=Nthin, flat=True)
print('len of extracted samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('best paras:',theta_max,np.max(probs))

figm, axm = plt.subplots(num='MF',figsize=(10, 10))
figl, axl = plt.subplots(num='LF',figsize=(10, 10))
curve_name = 'MF'
zs = np.arange(z,zmax,1)
rhoM = []; rhoM_min = []; rhoM_max = []
nL = []; nL_min = []; nL_max = []
# xs of models
Mx = M_BH[::int(N_mf/100)] # ~100 points
dlog10Mx = np.log10(Mx[1]/Mx[0])
for z in zs:
    #-------   MF   -------#
    draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
    thetas = samples[draw]
    model_thetas = [model(theta_i,z,f_seed,corr) for theta_i in thetas]
    mod_list = []
    for i in range(ndraw):
        mod = model_thetas[i]['MF'][::int(N_mf/100)]
        mod_list.append(mod)
        # axm.plot(Mx, mod, c='grey',label='_',alpha=.2)
    lnspread = np.std(np.log(mod_list),axis=0)
    lnmed = np.median(np.log(mod_list),axis=0)
    print('where lnspread nan?',Mx[np.isnan(lnspread)])
    print('where lnmed nan?',Mx[np.isnan(lnmed)])
    print('where lnspread inf?',Mx[np.isinf(lnspread)])
    print('where lnmed inf?',Mx[np.isinf(lnmed)])

    axm.fill_between(Mx,np.exp(lnmed-lnspread),np.exp(lnmed+lnspread),color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    #--- bestfit curve
    best_model = model(theta_max,z,f_seed,corr)
    y_best = best_model['MF'][::int(N_mf/100)]
    axm.plot(Mx, y_best, color='C%d'%(z%10), label='z=%d'%z)
    #--- mass integration above M_up
    dn_best = y_best*dlog10Mx
    index = np.where(Mx>M_up)
    rhoM.append(np.nansum(dn_best[index]*Mx[index]))
    dn_min = np.exp(lnmed-lnspread) *dlog10Mx
    rhoM_min.append(np.nansum(dn_min[index]*Mx[index]))
    dn_max = np.exp(lnmed+lnspread) *dlog10Mx
    rhoM_max.append(np.nansum(dn_max[index]*Mx[index]))
    # print(dlog10Mx*np.nansum((lnmed-lnspread)[Mx>M_up]),dlog10Mx*np.nansum(np.exp(lnmed+lnspread)[xs>M_up]))

    #-------   LF   -------#
    mod_list = []
    for i in range(ndraw):
        mod = model_thetas[i]['LF']
        mod_list.append(mod)
        # axl.plot(M1450, mod, c='grey',label='_',alpha=.2)
    lnspread = np.std(np.log(mod_list),axis=0)
    lnmed = np.median(np.log(mod_list),axis=0)
    axl.fill_between(M1450,np.exp(lnmed-lnspread),np.exp(lnmed+lnspread),color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    #--- bestfit curve
    y_best = best_model['LF']
    axl.plot(M1450, y_best, color='C%d'%(z%10), label='z=%d'%z)
    #--- integration brighter than L_up
    dn_best = y_best*dmag
    index = np.where(M1450<L_up)
    nL.append(np.nansum(dn_best[index]))
    dn_min = np.exp(lnmed-lnspread)*dmag
    nL_min.append(np.nansum(dn_min[index]))
    dn_max = np.exp(lnmed+lnspread)*dmag
    nL_max.append(np.nansum(dn_max[index]))
    print('after z{:d}, running time: {:.1f} min'.format(z,(time.time()-t1)/60))
# print(nL_min,nL_max)

axl.set_xlim(np.max(M1450),np.min(M1450))
axl.set_ylim(5e-3,1e2)
axl.set_yscale('log')
# axl.set_xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
# axl.set_ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
axl.tick_params(labelsize=fstick)
# plt.tight_layout()
figl.savefig(prex+'ndraw{:d}LF_spread.png'.format(ndraw),dpi=300,bbox_inches='tight')

axm.set_xlim(1e7,1e10); axm.set_xscale('log')
axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fstick)
figm.savefig(prex+'ndraw{:d}MF_spread.png'.format(ndraw),dpi=300,bbox_inches='tight')

# rhoM_min=rhoM; rhoM_max=rhoM
# fname = prex+'rhoM_evol.txt'
fname = prex+'ndraw{:d}rhoM_evol.txt'.format(ndraw)
ascii.write( Table([zs, rhoM, rhoM_min, rhoM_max],
            names=['zs','rhoM','rhoM_min','rhoM_max']),
            fname,
            formats={'zs':'5d','rhoM':'4.2e','rhoM_min':'4.2e','rhoM_max':'4.2e'},
            overwrite=True)

# nL_min=nL; nL_max=nL
# fname = prex+'nL_evol.txt'
fname = prex+'ndraw{:d}nL_evol.txt'.format(ndraw)
ascii.write( Table([zs, nL, nL_min, nL_max],
            names=['zs','nL','nL_min','nL_max']),
            fname,
            formats={'zs':'5d','nL':'4.2e','nL_min':'4.2e','nL_max':'4.2e'},
            overwrite=True)

# print('time=',time.time()-t1)