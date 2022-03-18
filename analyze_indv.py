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
t1 = time.time()

# prex='../MFcur7to10_err_quad_LFcur22.5to29_err1'
# prex='../MFcur7to10_err1_LFcur22.5to29_err1'
# prex='../MFcur7to10_err1_LFcur22.5to29_err0.5'
prex='../2p/2prange1_logM0=8_lcut=.9_a=.1/MLF2prange1_M8_l0_9.0e-01_a_1.0e-01'
fname = prex + '.h5'

print(fname,'logM0: ',logM0, 'l_cut: ',l_cut, 'a: ',a)

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

ndraw = 100
fig, axes = plt.subplots(1,2, figsize=(12, 6),dpi=400)
ax = axes[0]; curve_name = 'MF'
best_model = model(theta_max)
xs = best_model['M_BH']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.plot(xs, y_data, label='data')
draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
for i in thetas:
    mod = model(i)[curve_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='Highest Likelihood Model')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
ax.legend()
ax = axes[1]; curve_name = 'LF'
xs = best_model['M1450']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.scatter(xs, y_data, label='_')
for i in thetas:
    mod = model(i)[curve_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='_')
ax.text(-26,3e1,r'$t_{life}=$'+'{0:.1e}Myr\n'.format(theta_max[0])+r'$\delta=$'+'{0:.2f}'.format(theta_max[1]))
ax.set_xlim(-22,-29)
ax.set_ylim(1e-2,1e2)
ax.set_yscale('log')
ax.legend()
plt.savefig(prex+'_spread.png')

print('time:',time.time()-t1)