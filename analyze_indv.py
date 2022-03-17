from PYmodule import *
from PYmodule.models import *

# 1p
fname = '../1p/M30LF_1p_r_4even_ns5e+03.h5'
labels = ['t_life', 'prob']

# 2p
fname = '../2p/MLFnorm2prange1_l0_9.0e-01_a_1.0e-01.h5'
# fname = '../2p/MLF2prange1_r_4even_ns5.0e+03.h5'
# fname = '../2p/MLF2prange1_l0_2.0e+00_a_5.0e-02.h5'
# fname = '../2p/MLF2prange1_l0_5.0e-01_a_3.0e-01.h5'
print(fname, 'l_cut: ',l_cut, 'a: ',a)
t1 = time.time()

labels = ['t_life', 'd_fit', 'prob']

reader = emcee.backends.HDFBackend(fname)
prex = fname[:-3]
ndim = len(labels) - 1


tau = reader.get_autocorr_time()
tau = np.max(tau); print(tau)
Nburnin = int(3*tau)
Nthin = int(tau/2)

# samples, probs & theta_max --> prob_max
samples = reader.get_chain(discard=Nburnin, thin=Nthin, flat=True)
probs = reader.get_log_prob(discard=Nburnin, thin=Nthin, flat=True)
print('len of cut samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('PLOT: best paras:',labels,theta_max,np.max(probs))

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
ax = axes[0]; mod_name = 'MF'
best_model = model(theta_max)[mod_name]
xs = model(theta_max)['M_BH']
y_data = model(theta_max)[mod_name+'_data']
ax.plot(xs, y_data, label='data')
ndraw = 10
draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
for i in thetas:
    mod = model(i)[mod_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, best_model, c='C1', label='Highest Likelihood Model')

# ax.fill_between(xs,med_model-spread,med_model+spread,color='grey',alpha=0.5,label=r'$1\sigma$ Posterior Spread')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
ax.legend()
# ax.yaxis.set_label_coords(-0.1, 0.5)
ax = axes[1]; mod_name = 'LF'
best_model = model(theta_max)[mod_name]
xs = model(theta_max)['M1450']
y_data = model(theta_max)[mod_name+'_data']
ax.scatter(xs, y_data, label='_')
for i in thetas:
    mod = model(i)[mod_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
# ax.fill_between(xs,med_model-spread,med_model+spread,color='grey',alpha=0.5,label=r'$1\sigma$ Posterior Spread')
ax.plot(xs, best_model, c='C1', label='_')
ax.text(-26,1e2,r'$t_{life}=$'+'{0:.1e}Myr\n'.format(theta_max[0])+r'$\delta=$'+'{0:.2f}'.format(theta_max[1]))
ax.set_xlim(-22,-29)
ax.set_yscale('log')
# plt.savefig(fname[:6]+'model_spread_indv.png')
ax.legend()
plt.savefig(fname[:6]+'norml0_{0:.1e}_a_{1:.1e}'.format(l_cut,a)+'spread.png')
print(time.time()-t1)