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
ndraw = 60

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
prex = '../4p/M0r8_' + f_seedlabel
fname = prex + '.h5'
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
rhoM = []; rhoM_med = []; rhoM_std = []
nL = []; nL_med = []; nL_std = []
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
    spread = np.std(mod_list,axis=0)
    med_model = np.median(mod_list,axis=0)
    axm.fill_between(Mx,med_model-spread,med_model+spread,color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    #--- bestfit curve
    best_model = model(theta_max,z,f_seed,corr)
    y_best = best_model['MF'][::int(N_mf/100)]
    axm.plot(Mx, y_best, color='C%d'%(z%10), label='z=%d'%z)
    #--- mass integration above M_up
    dn_best = y_best*dlog10Mx
    index = np.where(Mx>M_up)
    rhoM.append(np.sum(dn_best[index]*Mx[index]))
    dn_med = med_model*dlog10Mx
    rhoM_med.append(np.sum(dn_med[index]*Mx[index]))
    dn_std = spread*dlog10Mx
    rhoM_std.append(np.sum(dn_std[index]*Mx[index]))
    # print(dlog10Mx*np.sum((med_model-spread)[Mx>M_up]),dlog10Mx*np.sum((med_model+spread)[xs>M_up]))

    #-------   LF   -------#
    mod_list = []
    for i in range(ndraw):
        mod = model_thetas[i]['LF']
        mod_list.append(mod)
        # axl.plot(M1450, mod, c='grey',label='_',alpha=.2)
    spread = np.std(mod_list,axis=0)
    med_model = np.median(mod_list,axis=0)
    axl.fill_between(M1450,med_model-spread,med_model+spread,color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    #--- bestfit curve
    y_best = best_model['LF']
    axl.plot(M1450, y_best, color='C%d'%(z%10), label='z=%d'%z)
    #--- integration brighter than L_up
    dn_best = y_best*dmag
    nL.append(np.sum(dn_best[M1450<L_up]))
    dn_med = med_model*dmag
    nL_med.append(np.sum(dn_med[M1450<L_up]))
    dn_std = spread*dmag
    nL_std.append(np.sum(dn_std[M1450<L_up]))
    print('after z{:d}, running time: {:.1f} min'.format(z,(time.time()-t1)/60))

axl.set_xlim(np.max(M1450),np.min(M1450))
axl.set_ylim(5e-3,1e2)
axl.set_yscale('log')
# axl.set_xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
# axl.set_ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
axl.tick_params(labelsize=fstick)
# plt.tight_layout()
figl.savefig(prex+'ndraw{:d}LF_spreadz{:d}.png'.format(ndraw,z),dpi=300,bbox_inches='tight')

axm.set_xlim(1e7,1e10); axm.set_xscale('log')
axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fstick)
figm.savefig(prex+'ndraw{:d}MF_spreadz{:d}.png'.format(ndraw,z),dpi=300,bbox_inches='tight')

# rhoM_med=rhoM; rhoM_std=rhoM
# fname = prex+'rhoM_evol.txt'
fname = prex+'ndraw{:d}rhoM_evol.txt'.format(ndraw)
ascii.write( Table([zs, rhoM, rhoM_med, rhoM_std],
            names=['zs','rhoM','rhoM_med','rhoM_std']),
            fname,
            formats={'zs':'5d','rhoM':'4.2e','rhoM_med':'4.2e','rhoM_std':'4.2e'},
            overwrite=True)

# nL_med=nL; nL_std=nL
# fname = prex+'nL_evol.txt'
fname = prex+'ndraw{:d}nL_evol.txt'.format(ndraw)
ascii.write( Table([zs, nL, nL_med, nL_std],
            names=['zs','nL','nL_med','nL_std']),
            fname,
            formats={'zs':'5d','nL':'4.2e','nL_med':'4.2e','nL_std':'4.2e'},
            overwrite=True)

# print('time=',time.time()-t1)