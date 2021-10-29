from PYmodule import *
# 读取z=6 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)

Chi2_min = {'no_corr':1e10, 'CO':1e10, 'DO':1e10}
Chi2_min = 1e10
z = int(6)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
Phi_obs = ma.masked_where(bin_cen>-23,Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]
print('bin_cen',bin_cen,'Phi_obs',Phi_obs,'Phi_err',Phi_err)

# # -----------------     obs data w/ error      --------------------
# Phi_obs_CO = Phi_obs/(1-f_obsc_const)
# Phi_obs_DO = Phi_obs*corr_U14D20(bin_cen)
# Phi_err_CO = Phi_err/(1-f_obsc_const)
# Phi_err_DO = Phi_err*corr_U14D20(bin_cen)

for f_duty in np.arange(.2, 1., .1): # .6 .4 
    for mu_fit in np.arange(.01, 1., .1): # f*mu .18, .19, .20
        for sigma_fit in np.arange(.08, 0.2, .01): # .10  .14
            fname = z6datapre+'LF11_'+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'alpha1'
            if os.path.isfile(fname):
                T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
                print(fname)
            else:
                continue

            Chi2 = np.sum(pow( (np.log(T['Phi']) - np.log(Phi_obs))/np.log(Phi_err), 2))
            print('f_duty%3.2f'%f_duty,'mu_fit%3.2f'%mu_fit, 'sigma_fit%3.2f'%sigma_fit, 'Chi2%.1e'%Chi2)

            if Chi2 < Chi2_min:
                Chi2_min = Chi2
                f_min = f_duty
                m_min = mu_fit
                s_min = sigma_fit
            Chi2_min = Chi2_min/(len(Phi_obs)-1)
            exit(0)
            # plt.figure(figsize=(10,10),dpi=200)
            # plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
            # plt.plot(bin_cen, LF_M1450_DO(bin_cen,z)*1e9, label=lfnames[z], c='C0')
            # plt.scatter(bin_cen, T['Phi'], c='C1', label='histogram')
            # plt.text(-26,1e3,r'$f='+str(f_duty)+'\n'+r'$\mu=$'+str(mu_fit)+'\n'+r'$\sigma=$'+str(sigma_fit),fontsize = fstxt)
            # plt.xlim(bin_cen[-1]+bin_wid[-1]/2.,bin_cen[0]-bin_wid[0]/2.); plt.ylim(5e-3,1e4)
            # plt.yscale('log')
            # plt.grid(True)
            # plt.legend(loc='lower left',fontsize=fslabel)
            # plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%6.1f'%Chi2,fontsize=fstitle)
            # plt.savefig(z6figpre+'chi2z'+str(z)+'f'+str(int(f_duty*10))+'s'+str(int(sigma*100))+'.png')

            # print('f_min:',f_min, 's_min',s_min, 'chi2_min:',Chi2_min)