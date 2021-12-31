from PYmodule import *
# 读取z=6 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)

Chi2_min = 1e10
z = int(6)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
M1450_min = -30.; M1450_max = -23.5
Phi_obs = Phi_obs[str(z)]
# Phi_obs = ma.masked_where(np.logical_or(bin_cen>M1450_max,bin_cen<M1450_min),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]
# print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)

# # -----------------     obs data w/ error      --------------------
# Phi_obs_CO = Phi_obs/(1-f_obsc_const)
# Phi_obs_DO = Phi_obs*corr_U14D20(bin_cen)
# Phi_err_CO = Phi_err/(1-f_obsc_const)
# Phi_err_DO = Phi_err*corr_U14D20(bin_cen)

fname = z6datapre+'LF_SC_t8.0e+02f0.8l1.0e-01a0.500alpha1.0'
# fname = z6datapre+'LF_LN_t2.0e+02f0.7m1.0e-02s1.000alpha1.0'
fname = z6datapre+'LF_LN_t9.0e+02f0.8m5.0e-02s0.300alpha1.0'
T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
Phi_DO = T['Phi_DO']
print(Phi_obs,Phi_DO,Phi_err)
Chi2 = np.sum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
print(Chi2, T['Chi2'][0])
Chi2 = T['Chi2'][0]

plt.figure(figsize=(10,10),dpi=200)
plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
plt.plot(bin_cen, LF_M1450(bin_cen,z)*1e9, label=lfnames[str(z)], c='C0')
plt.scatter(bin_cen, T['Phi_DO'], c='C1', label='histogram')
# plt.text(-26,1e2,'f=%3.2f'%f_duty+'\n'+r'$\mu=$'+'%3.2f'%mu_fit+'\n'+r'$\sigma=$'+'%3.2f'%sigma_fit+'\ne=%.3f'%eta8+'\n'+r'$\delta=$'+'%.3f'%delta_fit,fontsize = fstxt)
plt.xlim(bin_edg[-1],bin_edg[0]); plt.ylim(5e-3,1e4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='lower left',fontsize=fslabel)
plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%.2e'%Chi2,fontsize=fstitle)
plt.savefig(z6figpre+'ln.png')