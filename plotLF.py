from PYmodule import *
# 画几个best fitting para拟合的LF v.s. Matsu2018 data+fitting curve
# models: 'Schechter best fit','Log-normal best fit','Schechter + mass dependent ER'

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

fnames = [z6datapre+'LF_SC_t8.0e+02f0.8l1.0e-01a0.500alpha1.0',
          z6datapre+'LF_LN_t9.0e+02f0.8m5.0e-02s0.300alpha1.0',
          z6datapre+'LF_SC_t4.0e+02f0.6d1.0l8.0e-01a0.100alpha1.0']
labels = ['Schechter best fit','Log-normal best fit','Schechter + mass dependent ER']
# T = ascii.read(fname, guess=False, delimiter=' ')
# Phi_DO = T['Phi_DO']
# Chi2 = np.sum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
# assert Chi2 == T['Chi2'][0]

plt.figure(figsize=(10,10),dpi=200)
plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='black')
plt.plot(bin_cen, LF_M1450(bin_cen,z)*1e9, label=lfnames[str(z)], c='black')

for i in range(len(fnames)):
    fname = fnames[i]; lab = labels[i]
    T = ascii.read(fname, guess=False, delimiter=' ')
    plt.scatter(bin_cen, T['Phi_DO'], label=lab)

# plt.text(-26,1e2,'f=%3.2f'%f_duty+'\n'+r'$\mu=$'+'%3.2f'%mu_fit+'\n'+r'$\sigma=$'+'%3.2f'%sigma_fit+'\ne=%.3f'%eta8+'\n'+r'$\delta=$'+'%.3f'%delta_fit,fontsize = fstxt)
plt.xlim(bin_edg[-1],bin_edg[0]); plt.ylim(5e-3,1e4)
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right',fontsize=fslabel)
plt.xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel);
plt.ylabel(r'$\mathrm{Gpc^{-3} mag^{-1}}$',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
# plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%.2e'%Chi2,fontsize=fstitle)
plt.savefig(z6figpre+'SC_M.png')