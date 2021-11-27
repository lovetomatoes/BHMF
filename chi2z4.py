from PYmodule import *
# 读取z=4 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)

Chi2_min = 1e10
z = int(4)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
M1450_min = -30.; M1450_max = -22.
Phi_obs = ma.masked_where(np.logical_or(bin_cen>M1450_max,bin_cen<M1450_min),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]
# print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)

alpha = 1.
find_min = False
i = 0
for delta_fit in [0., 28, .29, .3, .31, .32]: # 0.3 among [0., .05, .1, .15, .2, .25, .3, .35, .4]
    for eta8 in [0., .08, .1, .12]:# .1
        if (eta8 == 0. and delta_fit == 0.):
            pass
        elif eta8*delta_fit == 0.:
            continue

        for f_duty in np.arange(.2, 1., .1): # .4
            for mu_fit in np.arange(.2, .6, .1): # .4
                for sigma_fit in np.arange(.01, 0.25, .01): # .1
                    i += 1
                    fname = z4datapre+'LF_'+'z%d'%z+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    if os.path.isfile(fname):
                        T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
                    else:
                        continue
                    Chi2 = np.sum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
                    if np.nanmin([Chi2, Chi2_min]) == Chi2:
                        find_min = True
                        Chi2_min = Chi2
                        f_min = f_duty
                        m_min = mu_fit
                        s_min = sigma_fit
                        e_min = eta8
                        d_min = delta_fit
                    # if i==2:
                    #     break

if find_min:
    print('f_min:',f_min,'m_min',m_min, 's_min',s_min, 'chi2_min:',Chi2_min)

f_duty = f_min; mu_fit = m_min; sigma_fit = s_min; eta8 = e_min; delta_fit = d_min
fname = z4datapre+'LF_'+'z%d'%z+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
Chi2 = np.sum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
print(T['Phi_DO'])
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
plt.savefig(z4figpre+'chi2_Phi_DO_'+'M1450min%.0f'%M1450_min+'max%.1fz'%M1450_max+str(z)+'MBHz6_2e10'+
                'f%3.2f'%f_duty+
                'm%3.2f'%mu_fit+
                's%3.2f'%sigma_fit+
                'e%.3f'%eta8+
                'd%.3f'%delta_fit+
                'alpha%.1f'%alpha+'.png')