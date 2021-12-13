from PYmodule import *
# 读取z=4 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 
# 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2) 
# 存 chi2 用1维数组  输出各个para和Chi2文件

Chi2_min = 1e10
z = int(4)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
M1450_min = -30.; M1450_max = -22. # z=4
Phi_obs = ma.masked_where(np.logical_or(bin_cen>M1450_max,bin_cen<M1450_min),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]

alpha = 1.
find_min = False

d_range = [0.001, .01, .1] # .001 in [0., .001, .01, .1]
e_range = [.05, .1, .15, .2] # .1 in [0., .05, .1, .15, .2]
f_range = np.arange(.5, 1., .1) # .7 in np.arange(.5, 1., .1)
m_range = np.arange(.2, .5, .01) # .21 in np.arange(.2, .5, .01)
s_range = np.arange(.01, 0.2, .01) # .15 in np.arange(.01, 0.2, .01)

# d_range = [0.001,.01,.1]
# e_range = [.1]
# f_range = np.arange(.3, 1., .1)
# m_range = np.logspace(-2, 0, num=5)
# s_range = np.arange(.1, 0.5, .1)

d_range = [.1]
e_range = [.05]
f_range = [10., 1.]
m_range =[.1, .01]
s_range = [.1]

Nd = len(d_range); Ne = len(e_range); Nf = len(f_range); Nm = len(m_range); Ns = len(s_range)
N_tot = Nd*Ne*Nf*Nm*Ns

Chi2s = np.ones(N_tot)
ds = np.ones(N_tot); es = np.ones(N_tot); fs = np.ones(N_tot); ms = np.ones(N_tot); ss = np.ones(N_tot)

i = -1
for delta_fit in d_range:
    for eta8 in e_range:
        for f_duty in f_range:
            for mu_fit in m_range:
                for sigma_fit in s_range:
                    i += 1
                    # continue
                    fname = z4datapre+'LF_'+'z'+str(z)+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    if os.path.isfile(fname):
                        T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
                    else:
                        # print('nofile',fname)
                        Chi2s[i] = np.nan
                        continue
                    Chi2 = np.nansum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
                    Chi2s[i] = Chi2
                    ds[i] = delta_fit; es[i] = eta8; fs[i] = f_duty; ms[i] = mu_fit; ss[i] = sigma_fit

                    if np.nanmin([Chi2, Chi2_min]) == Chi2:
                        find_min = True
                        Chi2_min = Chi2
                        f_min = f_duty
                        m_min = mu_fit
                        s_min = sigma_fit
                        e_min = eta8
                        d_min = delta_fit
                        i_min = i
print('i=',i)

if find_min:
    print('f_min%3.2f'%f_min+'m_min%3.2f'%m_min,'s_min%3.2f'%s_min+'e_min%.3f'%e_min+'d_min%.3f'%d_min)

T = Table([ds,es,fs,ms,ss,Chi2s],names=['ds','es','fs','ms','ss','Chi2s'])
ascii.write(T,z4datapre+'Chi2_para',formats={'ds':'5.1e','es':'5.1e','fs':'5.1e','ms':'5.1e','ss':'5.1e','Chi2s':'5.1e'},overwrite=True)

# range_min =  np.where(Chi2s<2e3*Chi2_min)
range_min =  np.where(Chi2s<1.1*Chi2_min)
print(range_min)