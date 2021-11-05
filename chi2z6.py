from PYmodule import *
# 读取z=6 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)

Chi2_min = {'no_corr':1e10, 'CO':1e10, 'DO':1e10}
Chi2_min = 1e10
z = int(6)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
Phi_obs = ma.masked_where(np.logical_or(bin_cen>-23,bin_cen<-30),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]
# print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)

# # -----------------     obs data w/ error      --------------------
# Phi_obs_CO = Phi_obs/(1-f_obsc_const)
# Phi_obs_DO = Phi_obs*corr_U14D20(bin_cen)
# Phi_err_CO = Phi_err/(1-f_obsc_const)
# Phi_err_DO = Phi_err*corr_U14D20(bin_cen)

find_min = False
for f_duty in np.arange(.2, 1., .1):
    for mu_fit in np.arange(.01, .5, .01):
        for sigma_fit in np.arange(.01, 0.2, .01):
            fname = z6datapre+'LF2e10_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'alpha1'
            if os.path.isfile(fname):
                T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
            else:
                continue
            Chi2 = np.sum(pow( (np.log(T['Phi']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
            if np.nanmin([Chi2, Chi2_min]) == Chi2:
                find_min = True
                Chi2_min = Chi2
                f_min = f_duty
                m_min = mu_fit
                s_min = sigma_fit

if find_min:
    print('f_min:',f_min,'m_min',m_min, 's_min',s_min, 'chi2_min:',Chi2_min)

f_duty = f_min; mu_fit = m_min; sigma_fit = s_min
fname = z6datapre+'LF11_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'alpha1'
T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
Chi2 = np.sum(pow( (np.log(T['Phi']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
print(T['Phi'])
plt.figure(figsize=(10,10),dpi=200)
plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
plt.plot(bin_cen, LF_M1450(bin_cen,z)*1e9, label=lfnames[str(z)], c='C0')
plt.scatter(bin_cen, T['Phi'], c='C1', label='histogram')
plt.text(-26,1e3,'f=%3.2f'%f_duty+'\n'+r'$\mu=$'+'%3.2f'%mu_fit+'\n'+r'$\sigma=$'+'%3.2f'%sigma_fit,fontsize = fstxt)
plt.xlim(bin_edg[-1],bin_edg[0]); plt.ylim(5e-3,1e4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%.2e'%Chi2,fontsize=fstitle)
plt.savefig(z6figpre+'chi2_no_z'+str(z)+'MBH2e10'+
                'f%3.2f'%f_duty+
                'm%3.2f'%mu_fit+
                's%3.2f'%sigma_fit+
                'alpha1.png')