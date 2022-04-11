from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

z0 = 35.
t0 = t_from_z(z0)
z1 = 6.
t_end = t_from_z(z1)

Dt = t_end - t0
M0 = 1e3


t_life, d_fit, l_cut, a = 1.68814272e+01, 7.31406767e-03, 1.02157385e+00, 1.46099557e-01 # f_seed = .01
t_life *= Myr

# table stores the cdf of lambda 
I_toinf = integral_toinf(a)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x)/I_toinf

N_BH = int(1e7)
N_BH = int(1e1)

Nt = int(Dt/t_life+1)

M1s = [M0*np.ones(N_BH)]
zs = [z0]
L1s = [L_M(M0,.01)*np.ones(N_BH)]

L_limit = 1e46 # 1e47

t = t0
for i in range(Nt):
    if t + t_life > t_end:
        dt = t_end - t
        t = t_end
    else:
        t = t + t_life
        dt = t_life
    uniform_a = np.random.uniform(size=N_BH)
    ls = np.zeros(N_BH)
    for i in range(N_BH):
        ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
    ls = ls*l_cut
    # xx,yy = np.meshgrid(uniform_a,Pa_x)
    # # argmax: the first yy>xx
    # ls = x[(yy>xx).argmax(axis=0)]

    M1 = M1M0_d(M0,ls,dt,d_fit)
    M0 = M1
    L1 = L_M(M1,ls)

print('exp:%.2e'%M1M0_e(1e3,Dt,1.))
print('d=%.2e'%d_fit,'into loop: %.2e'%np.mean(M1))

# for d_fit in np.logspace(-6,-1,num=6):
#     t = t0
#     M0 = 1.e3
#     for i in range(Nt):
#         if t + t_life > t_end:
#             dt = t_end - t
#             t = t_end
#         else:
#             t = t + t_life
#             dt = t_life
#         uniform_a = np.random.uniform(size=N_BH)
#         ls = np.zeros(N_BH)
#         for i in range(N_BH):
#             ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
#         ls = ls*l_cut
#         # xx,yy = np.meshgrid(uniform_a,Pa_x)
#         # # argmax: the first yy>xx
#         # ls = x[(yy>xx).argmax(axis=0)]
#         ## ---------------------- ##
#         ls = np.ones(len(ls))
#         ## ---------------------- ##
#         M1 = M1M0_d(M0,ls,dt,d_fit)
#         M0 = M1
#         L1 = L_M(M1,ls)
#     print('d=%f'%d_fit,'into loop: %.2e'%np.mean(M1))

# # z=6 BH mass, λ, L_bol
# ascii.write(Table([M1, ls, L1]),'../BHatz6.dat',names=['M1','ls','L1'],formats={'M1':'10.2e','ls':'10.2e','L1':'10.2e'},overwrite=True)
# print('time after evol:',time.time()-t1)

# # all BHs, all time: lambda distribution
# abin = np.log10(x)
# hist, bin_edges = np.histogram(np.log10(ls),bins=abin,density=False)
# plt.figure(figsize=(10,8),dpi=400)
# plt.scatter( bin_edges[:-1],hist/len(ls))
# # print(np.sum(hist)/len(ls))
# # print(np.sum((Pa_x[1:]-Pa_x[:-1])))
# plt.plot(np.log10(x[:-1])/2.,(Pa_x[1:]-Pa_x[:-1]),c='C1')
# plt.yscale('log')
# plt.savefig('../Plambda_all_MBHevol.png')