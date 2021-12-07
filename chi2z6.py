from PYmodule import *
# 读取z=6 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)

Chi2_min = {'no_corr':1e10, 'CO':1e10, 'DO':1e10}
Chi2_min = 1e10
z = int(6)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
M1450_min = -30.; M1450_max = -23.5
Phi_obs = ma.masked_where(np.logical_or(bin_cen>M1450_max,bin_cen<M1450_min),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]
# print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)

# # -----------------     obs data w/ error      --------------------
# Phi_obs_CO = Phi_obs/(1-f_obsc_const)
# Phi_obs_DO = Phi_obs*corr_U14D20(bin_cen)
# Phi_err_CO = Phi_err/(1-f_obsc_const)
# Phi_err_DO = Phi_err*corr_U14D20(bin_cen)

alpha = 1.
find_min = False

d_range = [0., 0.001, .01, .1] # .001 in [0., .001, .01, .1]
e_range = [0., .05, .1, .15, .2] # .1 in [0., .05, .1, .15, .2]
f_range = np.arange(.5, 1., .1) # .7 in np.arange(.5, 1., .1)
m_range = np.arange(.2, .5, .01) # .21 in np.arange(.2, .5, .01)
s_range = np.arange(.01, 0.2, .01) # .15 in np.arange(.01, 0.2, .01)

d_range = [0.001,.01,.1]
e_range = [.1]
f_range = np.arange(.3, 1., .1)
m_range = np.logspace(-2, 0, num=5)
s_range = np.arange(.1, 0.5, .1)

i = 0
for delta_fit in d_range:
    for eta8 in e_range:
        if (eta8 == 0. and delta_fit == 0.):
            pass
        elif eta8*delta_fit == 0.:
            continue
        for f_duty in f_range:
            for mu_fit in m_range:
                for sigma_fit in s_range:
                    # i += 1; continue
                    fname = z6datapre+'LF_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    if os.path.isfile(fname):
                        # i += 1; continue
                        T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
                    else:
                        # print('nofile',fname)
                        continue
                    Chi2 = np.nansum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
                    # print('chi2:',Chi2)
                    if np.nanmin([Chi2, Chi2_min]) == Chi2:
                        find_min = True
                        Chi2_min = Chi2
                        f_min = f_duty
                        m_min = mu_fit
                        s_min = sigma_fit
                        e_min = eta8
                        d_min = delta_fit
                    
# print('i=',i); exit(0)

if find_min:
    print('f_min%3.2f'%f_min+'m_min%3.2f'%m_min,'s_min%3.2f'%s_min+'e_min%.3f'%e_min+'d_min%.3f'%d_min)

f_duty = f_min; mu_fit = m_min; sigma_fit = s_min; eta8 = e_min; delta_fit = d_min
# f_duty = .7; mu_fit = .21
# for eta8 in [.05,.1,.2,]: 
#     for delta_fit in [.01,.02,.03,.04,.05,.1,.2,.3,.4,.5]: # f*mu .18, .19, .20
#         for sigma_fit in [.15]: # .10  .14

fname = z6datapre+'LF_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
Chi2 = np.sum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
# print(T['Phi_DO'])
plt.figure(figsize=(10,10),dpi=200)
plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
plt.plot(bin_cen, LF_M1450(bin_cen,z)*1e9, label=lfnames[str(z)], c='C0')
plt.scatter(bin_cen, T['Phi_DO'], c='C1', label='histogram')
plt.text(-26,1e2,'f=%3.2f'%f_duty+'\n'+r'$\mu=$'+'%3.2f'%mu_fit+'\n'+r'$\sigma=$'+'%3.2f'%sigma_fit+'\ne=%.3f'%eta8+'\n'+r'$\delta=$'+'%.3f'%delta_fit,fontsize = fstxt)
plt.xlim(bin_edg[-1],bin_edg[0]); plt.ylim(5e-3,1e4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%.2e'%Chi2,fontsize=fstitle)
plt.savefig(z6figpre+'chi2_Phi_DO_'+'M1450min%.0f'%M1450_min+'max%.1fz'%M1450_max+str(z)+
                'f%3.2f'%f_duty+
                'm%3.2f'%mu_fit+
                's%3.2f'%sigma_fit+
                'e%.3f'%eta8+
                'd%.3f'%delta_fit+
                'alpha%.1f'%alpha+'.png')
print('chi2_Phi_DO_'+'M1450min%.0f'%M1450_min+'max%.1fz'%M1450_max+str(z)+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha+'.png')
print('searched %d'%i)