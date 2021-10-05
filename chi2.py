from PYmodule import *

log10Ms = [9,10,11,12]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {4:'Akiyama_18',5:'McGreer_18',6:'Matsuoka_18'}
pres = ['./data/1e9','./data/1e10','./data/1e11','./data/1e12','./data/1e13']

i_bsm = 0
iM = 2

T_tell = 8000
eta = 0.3
Chi2_min = 1e10

# 读取z=4 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2)
z = int(6)

flambda = .18
f_duty = .5
for f_duty in [.4,.5,.6]:
    for sigma in [.1, .12, .14]: # .12, .15, .18
        fname = z6datapre+'Phi_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma*100))+'bsm01alpha1N4'
        if os.path.isfile(fname):
            T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
        else:
            continue
        
    # --------------------    CUT    -----------------------
    #    need to do in dist_grow ? wli
        # T = T[T['Mstar_z']<1e10]
    # --------------------    CUT    -----------------------

        # LF
        bin_cen = T['bin_cen']
        Phi_lf = T['Phi']
        bin_wid = np.array([2, 1, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1])

    # -----------------     obs data w/ error      --------------------
        Phi_obs = np.array([16.2, 23., 10.9, 8.3, 6.6, 7., 4.6, 1.33, .9, .58, .242, .0079])
        # print(Phi_obs)
        Phi_obs = Phi_obs[::-1]
        Phi_obs_CO = Phi_obs/(1-f_obsc_const)
        Phi_obs_DO = Phi_obs*corr_U14D20(bin_cen)
        Phi_err = np.array([16.2, 8.1, 3.6, 2.6, 2., 1.7, 1.2, .6, .32, .17, .061, .0079])
        Phi_err = Phi_err[::-1]
        Phi_err_CO = Phi_err/(1-f_obsc_const)
        Phi_err_DO = Phi_err*corr_U14D20(bin_cen)

        Phi_obs = Phi_obs_DO
        Phi_err = Phi_err_DO

        index = np.where(bin_cen<-24)
        index = np.logical_and(bin_cen<-24,bin_cen>-28)
        Phi_obs = Phi_obs[index]
        Phi_err = Phi_err[index]
        bin_cen = bin_cen[index]
        Phi_lf = Phi_lf[index]
        bin_wid = bin_wid[index]

        Chi2 = 0
        err = 1
        for i in range(len(Phi_obs)-1):
            if Phi_lf[i] != 0:
                # Chi2 += pow( (np.log(Phi_lf[i]) - np.log(Phi_obs[i]))*Phi_obs[i]/Phi_err[i], 2)
                Chi2 += pow( (np.log(Phi_lf[i]) - np.log(Phi_obs[i]))/np.log(Phi_err[i]), 2)
        print('flambda',flambda, 'sigma',sigma, Chi2)

        if Chi2 < Chi2_min:
            Chi2_min = Chi2
            f_min = flambda
            s_min = sigma 
        Chi2_min = Chi2_min/(len(Phi_obs)-1)
        # print(Phi_lf)
        plt.figure(figsize=(10,10),dpi=200)

        plt.errorbar(bin_cen, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
        plt.plot(bin_cen, LF_M1450_DO(bin_cen,z)*1e9, label=lfnames[z], c='C0')
        plt.scatter(bin_cen, Phi_lf, c='C1', label='histogram')
        plt.text(-26,1e3,r'$f\mu=$'+str(flambda)+';\nf='+str(f_duty)+'\n'+r'$\sigma=$'+str(sigma),fontsize = fstxt)
        plt.xlim(bin_cen[-1]+bin_wid[-1]/2.,bin_cen[0]-bin_wid[0]/2.); plt.ylim(5e-3,1e4)
        plt.yscale('log')
        plt.grid(True)
        plt.legend(loc='lower left',fontsize=fslabel)
        plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{\log{\sigma_i}^2}=$'+'%6.1f'%Chi2,fontsize=fstitle)
        plt.savefig(z6figpre+'DOmatsu_z'+str(z)+'fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma*100))+'.png')

        print('f_min:',f_min, 's_min',s_min, 'chi2_min:',Chi2_min)