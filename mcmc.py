from PYmodule import *
from PYmodule.mcz6 import *
from emcee import EnsembleSampler as EnsembleSampler
import corner

t_life = 120.
f_0  =  1.
d_fit = .25
l_cut = .9
a = .1

# x = (t_life, f_0, d_fit, l_cut, a)
# print(lnlike(x))
# exit(0)

t1 = time.time()
initial = np.array([t_life, f_0, d_fit, l_cut, a])
ndim = len(initial)
nwalkers = 10
p0 = [initial + 1e-2*initial*np.random.randn(ndim) for i in range(nwalkers)]
# print(p0)

t1 = time.time()
sampler = EnsembleSampler(nwalkers, ndim, lnprobab)
sampler.run_mcmc(p0, 10, progress=True)

print('nwalkers=%d'%nwalkers,'Delta time=',time.time()-t1)

labels = ['t_life', 'f_0', 'd_fit', 'l_cut', 'a']

fig, axes = plt.subplots(5, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
print('samples',samples)
probs = sampler.get_log_prob()
print('probs',probs)


labels = ['t_life', 'f_0', 'd_fit', 'l_cut', 'a']
# for i in range(ndim):
#     ax = axes[i]
#     ax.plot(samples[:, :, i], "k", alpha=0.3)
#     ax.set_xlim(0, len(samples))
#     ax.set_ylabel(labels[i])
#     ax.yaxis.set_label_coords(-0.1, 0.5)

# axes[-1].set_xlabel("step number")
# plt.savefig('paras.png')

samples = sampler.flatchain
print('len of samples:', len(samples))
theta_max = samples[np.argmax(sampler.flatlnprobability)]
print(np.max(sampler.flatlnprobability),labels,theta_max)
print('time',time.time()-t1)

tau = sampler.get_autocorr_time(tol=1)
print(tau)