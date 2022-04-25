from PYmodule import *
from PYmodule.l_intg import *

N_BH = int(1e5)
prex = z6datapre+'/f{0:.0e}N{1:d}'.format(f_seed,int(np.log10(N_BH)))
prex = '../'

# z=6 BH mass, λ, L_bol
T_z6 = ascii.read(prex+'BHatz6.dat', guess=False, delimiter=' ')
M1, ls, L1 = T_z6['M1'], T_z6['ls'], T_z6['L1']

T_evol = ascii.read(prex+'tN_evol.dat', guess=False, delimiter=' ')
zs, ts, N_bri = T_evol['zs'],T_evol['ts'],T_evol['N_bri']

T_M = np.loadtxt(prex+'Mevol.dat') # each line a BH growth track; all end with M>1e8
plt.figure(figsize=(10,10),dpi=400)
for i in range(len(T_M)):
    plt.plot(zs,T_M[i])
M1 = T_M[i]
i = np.argmax(M1); print('max M1[i]%.1e ?'%M1[i])
plt.plot(zs,T_M[i])
plt.yscale('log')
plt.xlim(np.max(zs),np.min(zs))
plt.ylim(1e2,1e11)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'_Mevol.png')

L_limit = 1e46
strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
T_hist = ascii.read(prex+strhist_+'.dat', guess=False, delimiter=' ')
log_l, hist_, hist_tot, ana = T_hist['log_l'],T_hist[strhist_],T_hist['hist_tot'],T_hist['ana']


exit(0)

## ----------------------------------------------------------------------------------------------
## wli: read z=6 file
# # plot λ dist.
# abin = np.log10(x)
# hist, bin_edges = np.histogram(np.log10(ls),bins=abin,density=False)
# plt.figure(figsize=(10,8),dpi=400)
# plt.scatter( bin_edges[:-1],hist/len(ls))
# # print(np.sum(hist)/len(ls))
# # print(np.sum((Pa_x[1:]-Pa_x[:-1])))
# plt.plot(np.log10(x[:-1]),(Pa_x[1:]-Pa_x[:-1]),c='C1')
# plt.yscale('log')
# plt.savefig('../Plambda_z6MBHevol.png')


lbin = np.linspace(-2,1.2,num=20)
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
x = np.logspace(np.log10(lambda_0),1.2,num=len(lbin))
# wli: temporary
l_cut, a = 1., .1
I_toinf = integral_toinf(a,lambda_0/l_cut)
Pana = integral(a,x)/I_toinf
ascii.write(Table([bin_edges[:-1],hist/len(ls),Pana[1:]-Pana[:-1]]),prex+'hist_tot.dat',
names=['log_l','hist','ana'],formats={'log_l':'10.2f','hist':'10.2e','ana':'10.2e'},overwrite=True)


# M1450=-26, Lbol=1e47; M1450=-22, Lbol=2e45
L_limit = 1e46
print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))

index = np.where(L_limit<L1)
M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
print('len(ls) after selection:',len(ls_),' / ',N_BH)
hist, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)

# histogram file
ascii.write(Table([bin_edges[:-1],hist/len(ls)]), prex+'hist_L%d.dat'%(int(np.log10(L_limit))),
names=['log_l','hist'],formats={'log_l':'10.2f','hist':'10.2e'},overwrite=True)
## ----------------------------------------------------------------------------------------------

## wli: change to read all time file
# all time samples
M1s = np.array(M1s)
L1s = np.array(L1s)
# print('N_act at zs', N_act)
ascii.write(Table([zs,ts,N_act]),'../Nact_evol.dat',names=['zs','ts','N_act'],formats={'zs':'10.2f','ts':'10.2e','N_act':'10.2e'},overwrite=True)

# BH mass evol.
M_low = 1e8
index = np.where(M1>M_low)
M1s = M1s.transpose()[index]
L1s = L1s.transpose()[index]
M1 = M1[index]

plt.figure(figsize=(10,10),dpi=400)
# print('M1s.transpose()[index]',M1s.transpose()[index].shape)
# print('M1s.transpose()[0]',M1s.transpose()[0].shape)
for i in range(len(M1)):
    plt.plot(zs,M1s[i])
i = np.argmax(M1); print('max M1[i]%.1e ?'%M1[i])
plt.plot(zs,M1s[i])
plt.yscale('log')
plt.xlim(z0,z1)
plt.ylim(1e2,1e11)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'_Mevol.png')

# # L_bol evol.
# index = np.logical_and(1e9<M0,M0<1e10)
# index = np.nonzero(index)[0]
# # print(len(index))

# plt.figure(figsize=(10,10),dpi=400)
# for i in index:
#     plt.plot(zs,L1s.transpose()[i])
# index = np.argmax(M0)
# plt.plot(zs,L1s.transpose()[index])
# plt.yscale('log')
# plt.xlim(10,z1)
# plt.ylim(1e41,1e48)
# plt.grid(True)
# plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
# plt.ylabel(r'$\mathrm{L_{bol}}$',fontsize=fslabel)
# # plt.legend(loc='best',fontsize=fslabel)
# plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
# plt.savefig(prex+'_L.png')