from PYmodule import *
from PYmodule.models_logd import *

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
fname = '../M0r8_'+f_seedlabel+'.h5'
prex = '../' + f_seedlabel

labels = [r'$\mathrm{t_{life}}$', r'$\log \delta$', r'$\lambda_0$', r'$\alpha$']

reader = emcee.backends.HDFBackend(fname)
ndim = len(labels)

tau = reader.get_autocorr_time(tol=1)
tau = np.max(tau); print('max tau:',tau)
Nburnin = int(3*tau)
Nthin = int(tau/2)

Nburnin = 0
Nthin = 1

samples = reader.get_chain(discard=Nburnin)
probs = reader.get_log_prob(discard=Nburnin)
print('len of samples:', len(samples))

samples = reader.get_chain(discard=Nburnin, thin=Nthin, flat=True)
probs = reader.get_log_prob(discard=Nburnin, thin=Nthin, flat=True)
print('len of extracted samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('best paras:',theta_max,np.max(probs))

# fig = corner.corner(samples,show_titles=True,title_kwargs={'fontsize':15},
# label_kwargs={'fontsize':20},max_n_ticks=4,labels=labels,plot_datapoints=True,
# quantiles=[0.16, 0.5, 0.84])
# axes = np.array(fig.axes).reshape((ndim, ndim))
# for i in range(ndim):
#     for j in range(ndim):
#         axes[i][j].tick_params(labelsize=12)
# plt.savefig(prex+'_corner.png',dpi=400,rasterized=True)
# exit(0)


ndraw = 60
fig, ax = plt.subplots(figsize=(10, 10))
curve_name = 'MF'
best_model = model(theta_max)
xs = best_model['M_BH'][::int(N_mf/100)] # ~100 points
y_data = best_model[curve_name+'_data'][::int(N_mf/100)]
y_logdata = np.log10(y_data)
y_best = best_model[curve_name][::int(N_mf/100)]
y_err = best_model[curve_name+'_data_err'][::int(N_mf/100)]
ax.plot(xs, y_data, c='C0',label='W10')
# error band of W10
ax.fill_between(xs,pow(10.,y_logdata-y_err/2.),pow(10.,y_logdata+y_err/2.),color='C0',alpha=0.3,label='data error')

draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
model_thetas = [model(theta_i) for theta_i in thetas]
mod_list = []
for i in range(ndraw):
    mod = model_thetas[i][curve_name][::int(N_mf/100)]
    mod_list.append(mod)
    # ax.plot(xs, mod, c='grey',label='_',alpha=.2)
spread = np.std(mod_list,axis=0)
med_model = np.median(mod_list,axis=0)
ax.fill_between(xs,med_model-spread,med_model+spread,color='purple',alpha=0.3,label=r'$1\sigma$ Posterior Spread')
ax.plot(xs, y_best, c='purple', label='Highest Likelihood Model')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
ax.legend(fontsize=fslegend)
plt.xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
plt.tick_params(labelsize=fstick)
plt.savefig(prex+'ndraw%dMF_spread.png'%ndraw,dpi=300,bbox_inches='tight')

fig, ax = plt.subplots(figsize=(10, 10))
curve_name = 'LF'
xs = best_model['M1450']
x_data = best_model['M1450_data']
y_data = best_model[curve_name+'_data']
y_data_err = best_model[curve_name+'_data_err']
y_best = best_model[curve_name]
ax.scatter(x_data, y_data, label='_')
plt.errorbar(x_data, y_data, yerr=y_data_err,fmt='o',capsize=10)
# print('y_data_err',y_data_err)
mod_list = []
for i in range(ndraw):
    mod = model_thetas[i][curve_name]
    mod_list.append(mod)
    # ax.plot(xs, mod, c='grey',label='_',alpha=.2)
spread = np.std(mod_list,axis=0)
med_model = np.median(mod_list,axis=0)
plt.fill_between(xs,med_model-spread,med_model+spread,color='purple',alpha=0.3,label='_')
ax.plot(xs, y_best, c='purple', label='_')
ax.text(-26,5, f_seedlabel+'\n'+ \
labels[0]+' = %.2f Myr\n'%(theta_max[0])+labels[1]+' = %.2f\n'%(theta_max[1])\
+labels[2]+' = %.2f\n'%(theta_max[2])+labels[3]+' = %.2f\n'%(theta_max[3])
, fontsize=20)
ax.set_xlim(np.max(xs),np.min(xs))
ax.set_ylim(5e-3,1e2)
ax.set_yscale('log')
plt.xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
plt.tick_params(labelsize=fstick)
# plt.tight_layout()
plt.savefig(prex+'ndraw%dLF_spread.png'%ndraw,dpi=300,bbox_inches='tight')

'''
ndraw = 2
fig, ax = plt.subplots(figsize=(10, 10))
curve_name = 'MF'
best_model = model(theta_max)
xs = best_model['M_BH']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.plot(xs, y_data, label='data')
draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
model_thetas = [model(theta_i) for theta_i in thetas]
for i in range(ndraw):
    mod = model_thetas[i][curve_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='Highest Likelihood Model')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
ax.legend(fontsize=fslegend)
plt.xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
plt.tick_params(labelsize=fstick)
plt.savefig(prex+'MF_spread.png',dpi=300,bbox_inches='tight')

# def sample_walkers(nsamples,flattened_chain):
#     models = []
#     draw = np.floor(np.random.uniform(0,len(flattened_chain),size=nsamples)).astype(int)
#     thetas = flattened_chain[draw]
#     for i in thetas:
#         mod = model(i)
#         models.append(mod)
#     spread = np.std(models,axis=0)
#     med_model = np.median(models,axis=0)
#     return med_model,spread
# med_model, spread = sample_walkers(100,)
# plt.fill_between(xs,med_model-spread,med_model+spread,color='grey',alpha=0.3,label=r'$1\sigma$ Posterior Spread')

fig, ax = plt.subplots(figsize=(10, 10))
curve_name = 'LF'
xs = best_model['M1450']
x_data = best_model['M1450_data']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.scatter(x_data, y_data, label='_')
# models = []
for i in range(ndraw):
    mod = model_thetas[i][curve_name]
    # models.append(mod)
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='_')
ax.text(-26,10,
labels[0]+' = %.2f Myr\n'%(theta_max[0])+labels[1]+' = %.2f\n'%(theta_max[1])\
+labels[2]+' = %.2f\n'%(theta_max[2])+labels[3]+' = %.2f\n'%(theta_max[3])
, fontsize=20)
ax.set_xlim(-22,-29)
ax.set_ylim(1e-2,1e2)
ax.set_yscale('log')
plt.xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
plt.tick_params(labelsize=fstick)
# plt.tight_layout()
plt.savefig(prex+'LF_spread.png',dpi=300,bbox_inches='tight')
'''