from PYmodule import *
from PYmodule.models import *

fname = '../1p/M30LF_1p_r_4even_ns5e+03.h5'
reader = emcee.backends.HDFBackend(fname)
prex = fname[:-3]
labels = ['t_life', 'prob']
ndim = len(labels) - 1


tau = reader.get_autocorr_time()
print(tau)
Nburnin = int(3*tau)
Nthin = int(tau/2)


fig, axes = plt.subplots(ndim+1, figsize=(10, 7), sharex=True)
samples = reader.get_chain(discard=Nburnin)
probs = reader.get_log_prob(discard=Nburnin)

for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
i += 1
ax = axes[i]
ax.plot(probs[:, :], "k", alpha=0.3)
ax.set_xlim(0, len(samples))
ax.set_ylabel(labels[i])
ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")
plt.savefig(prex+'_anachain.png')

# best theta_max --> prob_max
samples = reader.get_chain(flat=True)
probs = reader.get_log_prob(flat=True)
print('len of samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('best paras:',labels,theta_max,np.max(probs))

# 
samples = reader.get_chain(discard=Nburnin, thin=Nthin, flat=True)
probs = reader.get_log_prob(discard=Nburnin, thin=Nthin, flat=True)
print('len of cut samples:', len(samples))

theta_max = samples[np.argmax(probs)]
print('PLOT: best paras:',labels,theta_max,np.max(probs))


all_samples = np.concatenate(
    (samples, probs[:, None]), axis=1
)
fig = corner.corner(all_samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
plt.savefig(prex+'_anacorner.png')



def sample_walkers(nsamples,flattened_chain,mod_name):
    mod_collection = []
    draw = np.floor(np.random.uniform(0,len(flattened_chain),size=nsamples)).astype(int)
    thetas = flattened_chain[draw]
    for i in thetas:
        mod = model(i)[mod_name]
        mod_collection.append(mod)
    spread = np.std(mod_collection,axis=0)
    med_model = np.median(mod_collection,axis=0)
    return med_model,spread


fig, axes = plt.subplots(1,2, figsize=(12, 6),dpi=400)
ax = axes[0]
mod_name = 'MF'
best_model = model(theta_max)[mod_name]
xs = model(theta_max)['M_BH']
y_data = model(theta_max)[mod_name+'_data']
med_model, spread = sample_walkers(10,samples,mod_name)
ax.plot(xs, y_data, label='data')
ax.plot(xs, best_model, c='C1', label='Highest Likelihood Model')
ax.fill_between(xs,med_model-spread,med_model+spread,color='grey',alpha=0.5,label=r'$1\sigma$ Posterior Spread')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
# ax.yaxis.set_label_coords(-0.1, 0.5)

ax = axes[1]
mod_name = 'LF'
best_model = model(theta_max)[mod_name]
xs = model(theta_max)['M1450']
y_data = model(theta_max)[mod_name+'_data']
med_model, spread = sample_walkers(10,samples,mod_name)
ax.scatter(xs, y_data, label='data')
ax.plot(xs, best_model, c='C1', label='Highest Likelihood Model')
ax.fill_between(xs,med_model-spread,med_model+spread,color='grey',alpha=0.5,label=r'$1\sigma$ Posterior Spread')
ax.set_xlim(-22,-29)
ax.set_yscale('log')
plt.savefig('../1p/model_spread.png')