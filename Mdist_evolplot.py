from PYmodule import *
from PYmodule.l_intg import *
# analyse BHMF, ERDF from sampling in Mdist_evol
# later plot hist file from Mdist_evol directly 

T_seed = Ts[0][0]
tz = t_from_z(6.)/Myr
t0 = np.min(T_seed['t_col']/Myr)
# print(tz,t_life,t0);exit(0)

# wli: paras must set the same as Mdist_evol
prex = '../4p/distf1N5_07241632_'
t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
prex = '../4p/distf2N4_06221827_'
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01

# z=6 BH mass, λ, L_bol
T_z6 = ascii.read(prex+'BHatz6.txt', guess=False, delimiter=' ')
f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))

M1, ls, L1 = T_z6['M1'], T_z6['ls'], T_z6['L1']

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x,x0)/I_toinf

lbin = np.linspace(-2,1.2,num=20)
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf

SMl1 = ls[M1>1e9] # super-massive BHs >1e9

# T_evol = ascii.read(prex+'tN_evol.dat', guess=False, delimiter=' ')
# zs, ts, N_bri = T_evol['zs'],T_evol['ts'],T_evol['N_bri']

T_M = np.loadtxt(prex+'Mevol.txt') # each line a BH growth track; all end with M>1e8
T_l = np.loadtxt(prex+'levol.txt') # each line a lambda track; all end with M>1e8

# growth tracks for all BHs(>1e8)
plt.figure(figsize=(10,10),dpi=400)
for i in range(len(T_M)):
    # T_M[i] finally stay at final mass with dt=0; extract growing tracks
    len_series = np.argmax(T_M[i]) + 1
    # checkpoint series, back to t_seed (ideal)
    ts = np.arange(tz,t0-t_life,-t_life)[:len_series]
    zs = z_tH(ts)[::-1]
    plt.plot(zs,T_M[i][:len_series], c='grey', alpha=0.5)

M1 = T_M.transpose()[-1]
# print('M1',M1,T_M.shape,'len zs',len(zs))
# i = np.argmax(M1); print('max M1 = %.1e'%M1[i])
# plt.plot(zs,T_M[i], '--',c='black',lw=2)

SMM1 = M1[M1>1e9]
# print('number of >1e9 MBH:',sum(M1>1e9))
# same as prev. specially for >1e9 SMBHs
for i in range(len(SMM1)):
    len_series = np.argmax(T_M[M1>1e9][i]) + 1
    ts = np.arange(tz,t0-t_life,-t_life)[:len_series]
    zs = z_tH(ts)[::-1]
    ts = ts[::-1]
    plt.plot(zs,T_M[M1>1e9][i][:len_series],c='C%d'%i)
    plt.plot(zs,M1M0_e(SMM1[i],-(ts[-1]-ts)*Myr,SMl1[i]),'--',c='C%d'%i)
    plt.plot(zs,M1M0_e(SMM1[i],-(ts[-1]-ts)*Myr,1),'.',c='C%d'%i)
    plt.scatter(6,SMM1[i],marker='D',s=30,c='C%d'%i)

plt.yscale('log')
plt.xlim(np.max(zs),np.min(zs))
plt.xlim(30,5.8)
plt.ylim(1e2,1e10)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'Mevol.png', dpi=300, bbox_inches='tight')

# f_duty (λ>1) in growth history, v.s. final M_BH
plt.figure(figsize=(10,10),dpi=400)
f_duty = np.zeros(len(T_M))
f_2 = 0
t_start = t_from_z(10.)/Myr
for i in range(len(T_M)):
    len_series = np.argmax(T_M[i]) + 1
    ts = np.arange(tz,t0-t_life,-t_life)[:len_series]
    l1s = T_l[i][:len_series]
    # f_duty: Δt of λ>1 cycles, devided by total Δt
    l_crit = 1.
    # l_crit = l_cut
    f_duty[i] = t_life*np.sum(l1s>l_crit)/(ts[0]-ts[-1])
    f_2 += t_life * np.sum(l1s[ts>t_start]>l_crit)/(tz-t_start)
print(f_seedlabel+' duty cycle l>{:.2f} from z=10: {:.3f}'.format(l_crit,f_2/len(T_M)))

plt.scatter(M1,f_duty)
plt.xscale('log')
# plt.ylim(1e2,1e10)
plt.grid(True)
plt.xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{f_{duty}~(\lambda>1)}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'f_duty.png', dpi=300, bbox_inches='tight')


# dP/dlnλ, in growth history; selected by final M_BH>1e8, 1e9
# biased to high λ than Schechter model Pλ 
plt.figure(figsize=(10,10),dpi=400)
lbin = np.linspace(-2,1.2,num=20)
hist_M8 = np.zeros(len(lbin)-1)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf

for i in range(len(T_M)):
    len_series = np.argmax(T_M[i]) + 1
    ts = np.arange(tz,t0-t_life,-t_life)[:len_series]
    zs = z_tH(ts)[::-1]
    l1s = T_l[i][:len_series]
    hist, bin_edges = np.histogram(np.log10(l1s),bins=lbin,density=False)
    hist_M8 += hist/len_series/len(T_M)
hist_M9 = np.zeros(len(lbin)-1)
print('SMBHs>1e9, #=%d with l at z=6:'%len(SMM1))
print(T_l[M1>1e9].transpose()[-1])
for i in range(len(SMM1)):
    len_series = np.argmax(T_M[M1>1e9][i]) + 1
    ts = np.arange(tz,t0,-t_life)[:len_series]
    zs = z_tH(ts)[::-1]
    l1s = T_l[M1>1e9][i][:len_series]
    hist, bin_edges = np.histogram(np.log10(l1s),bins=lbin,density=False)
    hist_M9 += hist/len_series/len(SMM1)

plt.plot(lbin[:-1],Pana[1:]-Pana[:-1],c='black')
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
plt.bar(lbin[:-1], hist/len(ls),width=lbin[1]-lbin[0],fill=0,alpha=0.5,label='M9')
plt.bar(lbin[:-1], hist_M9,width=lbin[1]-lbin[0],color='C'+str(0),alpha=0.5,label='M9')
plt.bar(lbin[:-1], hist_M8,width=lbin[1]-lbin[0],color='C'+str(1),alpha=0.5,label='M8')

plt.yscale('log')
plt.xlabel(r'$\mathrm{\log\lambda}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{dP/d\log\lambda}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'lambda_hist.png', dpi=300, bbox_inches='tight')

exit(0)

# plot z=6 λ dist.
L_limit = 1e45
strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
lbin = np.linspace(-2,1.2,num=20)
T_hist = ascii.read(prex+strhist_+'.dat', guess=False, delimiter=' ')
left_of_lbin, hist_, hist_tot, ana = T_hist['log_l'],T_hist[strhist_],T_hist['hist_tot'],T_hist['ana']

print('sum(hist_tot)={0:.2e},sum(hist_):{1:.2e}'.format(np.sum(hist_tot),np.sum(hist_)))

plt.figure(figsize=(10,9))
x = (lbin[:-1]+lbin[1:])/2.
wid = lbin[1:]-lbin[:-1]
plt.bar(x, hist_tot,width=wid,color='C'+str(0),alpha=0.5,label='total')
plt.bar(x, 10*hist_,   width=wid,color='C'+str(1),alpha=0.5,label=r'$L>10^{45}$ erg/s')
# plt.bar(x, hist_,   width=wid,color='C'+str(1),alpha=0.5,label=r'$L>10^{46}$ erg/s')
# plt.yscale('log'); plt.ylim(bottom=1e-5)

# plot z=6 λ dist.
L_limit = 1e46
strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
lbin = np.linspace(-2,1.2,num=20)
T_hist = ascii.read(prex+strhist_+'.dat', guess=False, delimiter=' ')
left_of_lbin, hist_, hist_tot, ana = T_hist['log_l'],T_hist[strhist_],T_hist['hist_tot'],T_hist['ana']
plt.bar(x, 100*hist_,   width=wid,color='C'+str(2),alpha=0.5,label=r'$L>10^{46}$ erg/s')

# index = np.logical_and(1e7<x, x<1e10)
plt.plot(x, ana, '--',c='black',lw=2, label=r'$\mathrm{P(\lambda)}$')
plt.tick_params(labelsize=fstick)
plt.xlabel(r'$\log\lambda$',fontsize=fslabel)
# plt.ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{\left(Mpc^{-3}dex^{-1}\right)}$',fontsize=fslabel)
# plt.xlim(1e2,1e10); plt.ylim(1e-10,1e-2)
plt.legend(fontsize=fslegend,loc='best')
plt.savefig(prex+'l_hist.png', dpi=300, bbox_inches='tight')