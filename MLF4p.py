# Gauss prior for parameter l_cut and a;

from PYmodule import *
from PYmodule.MLF4p import *
from emcee import EnsembleSampler as EnsembleSampler
import corner
import os
os.environ["OMP_NUM_THREADS"] = "1"
from schwimmbad import MPIPool

t1 = time.time()

# initial guess
t_life = 120.
d_fit = 0.25
l_cut = l_mean # defined in PYmodule/__init__.py, l_mean and a_mean
a = a_mean

initial = np.array([t_life,d_fit,l_cut,a])

ndim = len(initial)
nwalkers = 100
nsteps = 5000
rball = 1e-4

prex='MLF4prange1_Gs.2_l0_{0:.1e}_a_{1:.1e}'.format(l_cut,a)

fname =prex+'.h5'

# nsteps = 10000
# prex += '_xu'

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    np.random.seed(42)

    backend = emcee.backends.HDFBackend(fname)

    # --initialize-- clear  output and reset
    backend.reset(nwalkers, ndim)
    p0 = [initial + rball*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = EnsembleSampler(nwalkers, ndim, lnprobab, pool=pool, backend=backend)
    sampler.run_mcmc(p0, nsteps, progress=True)

    # # --resume-- 
    # sampler = EnsembleSampler(nwalkers, ndim, lnprobab, pool=pool, backend=backend)
    # sampler.run_mcmc(None, nsteps, progress=True)


fig, axes = plt.subplots(ndim+1, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
# print('samples',samples)
probs = sampler.get_log_prob()
# print('probs',probs)

labels = ['t_life', 'd_fit', 'l_cut', 'a', 'prob']
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

plt.savefig(prex+'_chain.png')


samples = sampler.flatchain
probs = sampler.flatlnprobability
print('len of samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('best paras:',labels,theta_max,np.max(probs))

all_samples = np.concatenate(
    (samples, probs[:, None]), axis=1
)
fig = corner.corner(all_samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
plt.savefig(prex+'_corner.png')

# pyname = sys.argv[0][:-3] # current .py file name
print('nwalkers={0:d}, nsteps={1:.0e}, rball={2:.0e}'.format(int(nwalkers),int(nsteps),rball))

print(
    "Mean acceptance fraction: {0:.3f}".format(
        np.mean(sampler.acceptance_fraction))
)

tau = sampler.get_autocorr_time(tol=1)
print(tau)

print('l_mean={0:.1e}, sigma_l={1:.1e}, a_mean={2:.1e}, sigma_a={3:.1e}'.format(l_mean,sigma_l,a_mean,sigma_a))
print('running time: {0:.1f} hrs'.format((time.time()-t1)/3600))