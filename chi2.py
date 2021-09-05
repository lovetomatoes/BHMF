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
z = int(4)
Nsite = int(1e7)

for flambda in [0.11, .12, .13]:
    for sigma in [.15]: # .12, .15, .18
        fname = z4datapre+'LF_z'+str(z)+'f'+str(int(flambda*100))+'s'+str(int(sigma*100))+'N3'

        T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
        Nsite = len(T)

    # --------------------    CUT    -----------------------
        T = T[T['Mstar_z']<1e10] # wli: added
        # LF
        abin_lf =np.linspace(-29, -21.75,num=30)
        wid_lf = abs(abin_lf[1:]-abin_lf[:-1])
        dmag = abin_lf[1]-abin_lf[0] # print(dmag)

        Phi_lf, bin_edges = np.histogram(T['M1450_z'],bins=abin_lf,density=False)
        f_duty = .5
        Phi_lf = Phi_lf * n_base[iM]/float(Nsite)/dmag*f_duty*1e9
    # -----------------     obscured fraction correction      --------------------
        # #                   canceled 
        # x = (abin_lf[:-1]+abin_lf[1:])/2.
        # M_1 = -27.2
        # M_2 = -20.7
        # t = (x - M_1) / (M_2 - M_1)
        # Phi_lf = Phi_lf/(2*(1-t)+3*t)

        Phi_obs = 10**np.array([-6.253,-6.088,-6.219,-6.298,-6.382,-6.297,-6.369,-6.341,-6.405,-6.493,-6.487,-6.577,-6.745,-6.756,-6.899,-7.189,-7.077,
        -7.387,-7.552,-7.642,-7.853,-7.966,-8.192,-8.534,-8.781,-9.057,-9.710,-9.534,-9.710]) * 1e9
        # print(Phi_obs)
        Phi_obs = Phi_obs[::-1]
        Phi_err = np.array([80.733,71.892,59.367,52.208,46.855,51.538,47.366,47.960,42.783,37.343,37.385,33.384,27.422,27.088,23.032,17.287,29.581,
        2.545,1.657,1.490,1.169,1.026,0.791,0.534,0.402,0.292,0.138,0.169,0.138])
        Phi_err = Phi_err[::-1]

        Chi2 = 0
        err = 1
        for i in range(len(abin_lf)-1):
            if Phi_lf[i] != 0:
                Chi2 += pow( (np.log(Phi_lf[i]) - np.log(Phi_obs[i]))*Phi_obs[i]/Phi_err[i], 2)
            # Chi2 += pow( (Phi_lf[i] - Phi_obs[i])/Phi_err[i], 2)
        print('flambda',flambda, 'sigma',sigma, Chi2)

        if Chi2 < Chi2_min:
            Chi2_min = Chi2
            f_min = flambda
            s_min = sigma 
        # print(Phi_lf)
        plt.figure(figsize=(10,10),dpi=200)
        x = (abin_lf[1:]+abin_lf[:-1])/2

        plt.errorbar(x, Phi_obs, yerr=Phi_err, fmt="o", c='C0')
        plt.plot(x, LF_M1450(x,z)*1e9, label=lfnames[z], c='C0')
        plt.scatter(x, Phi_lf, c='C1', label='histogram')
        plt.text(-26,1e3,r'$f\mu=$'+str(flambda)+'; '+r'$\sigma=$'+str(sigma),fontsize = fstxt)
        plt.xlim(-22,-29); plt.ylim(1e-2,1e4)
        plt.yscale('log')
        plt.legend(loc='lower left',fontsize=fslabel)
        plt.title(r'$\mathrm{\chi^2}=\sum_i \frac{(\log{E_i}-\log{O_i})^2}{(\sigma_i/O_i)^2}=$'+str(int(Chi2)),fontsize=fstitle)
        plt.savefig(z4figpre+'akiz'+str(z)+'f'+str(int(flambda*100))+'s'+str(int(sigma*100))+'.png')

print('f_min:',f_min, 's_min',s_min, 'chi2_min:',Chi2_min)