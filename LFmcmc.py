from PYmodule import *
from PYmodule.LFmcz6 import *
from emcee import EnsembleSampler as EnsembleSampler

t_life = 120.
f_0  =  1.
d_fit = .25
l_cut = .9
a = .1
# -3.72

# t_life = 2.82186941e+02 
# f_0 = 9.82415200e-01 
# d_fit = 4.74430660e-01 
# l_cut = 9.81936736e-01
# a = 2.40364785e-02
# # -1.48

# t_life = 133.97
# f_0  = 1.07
# d_fit = .19
# l_cut = .88
# a = .09
# # -6.06

# x = (t_life, f_0, d_fit, l_cut, a)
# print(lnlike(x))
# exit(0)

t1 = time.time()
initial = np.array([t_life, f_0, d_fit, l_cut, a])
ndim = len(initial)
nwalkers = 10
nsteps  = 10
p0 = [initial + 1e-2*initial*np.random.randn(ndim) for i in range(nwalkers)]

t1 = time.time()
sampler = EnsembleSampler(nwalkers, ndim, lnprobab)
sampler.run_mcmc(p0, nsteps, progress=True)

print('nwalkers=%d'%nwalkers,'Delta time=',time.time()-t1)

print('samples',sampler.flatchain)
probs = sampler.flatlnprobability
print('probs',probs)
print('time',time.time()-t1)

tau = sampler.get_autocorr_time(tol=1)
print(tau)