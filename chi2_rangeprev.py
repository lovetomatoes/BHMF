from PYmodule import *
# 读取z=4 的M1450 分布 （文件名带着参数信息） 对于任意参数组合 
# 返回chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2) 
# 存 chi2 用5维数组  输出？分析？

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

Nd = len(d_range)
Ne = len(e_range)
Nf = len(f_range)
Nm = len(m_range)
Ns = len(s_range)
N_tot = Nd*Ne*Nf*Nm*Ns

i = -1
Chi2s = np.ones(N_tot)
for id in range(Nd):
    delta_fit = d_range[id]
    for ie in range(Ne):
        eta8 = e_range[ie]
        for i_f in range(Nf):
            f_duty = f_range[i_f]
            for im in range(Nm):
                mu_fit = m_range[im]
                for i_s in range(Ns):
                    sigma_fit = s_range[i_s]
                    i += 1
                    fname = z4datapre+'LF_'+'z'+str(z)+'f%3.2f'%f_duty+'m%3.2f'%mu_fit+'s%3.2f'%sigma_fit+'e%.3f'%eta8+'d%.3f'%delta_fit+'alpha%.1f'%alpha
                    if os.path.isfile(fname):
                        T = ascii.read(fname, guess=False, delimiter=' ') #  None has np.where(T['z_col']==-1)
                    else:
                        # print('nofile',fname)
                        # Chi2s[id][ie][i_f][im][i_s] = np.nan
                        Chi2s[i] = np.nan
                        continue
                    Chi2 = np.nansum(pow( (np.log(T['Phi_DO']) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
                    Chi2s[i] = Chi2
                    print('Chi2=',Chi2)
                    print('i=%d;'%(i-1),id,ie,i_f,im,i_s,'argmin=',id+ie+i_f+im+i_s)
                    if np.nanmin([Chi2, Chi2_min]) == Chi2:
                        find_min = True
                        Chi2_min = Chi2
                        f_min = f_duty
                        m_min = mu_fit
                        s_min = sigma_fit
                        e_min = eta8
                        d_min = delta_fit
                        # print(id,ie,i_f,im,i_s,'argmin=',id+ie+i_f+im+i_s)
                        # i_min = np.argmin(Chi2s)
                        # print('argmin:',i_min) # !!! degeneracy
                        # print(i_min%Ns, (i_min-i_min%Ns)%Nm, i_min-(i_min-i_min%Ns))
                    

for i in range(N_tot):
    pass
# print('i=',i, Nd*Ne*Nf*Nm*Ns)
print(Chi2s.shape)

# if find_min:
#     print('f_min%3.2f'%f_min+'m_min%3.2f'%m_min,'s_min%3.2f'%s_min+'e_min%.3f'%e_min+'d_min%.3f'%d_min)

# print('argmin:',np.argmin(Chi2s)) # !!! degeneracy
range_min =  np.where(Chi2s<2e3*Chi2_min)
print(range_min)
print([range_min])
# print( np.mean(range_min[0]),np.std(range_min[0]) )
# print( np.mean(range_min[1]),np.std(range_min[1]) )
# print( np.mean(range_min[2]),np.std(range_min[2]) )
# print( np.mean(range_min[3]),np.std(range_min[3]) )
# print( np.mean(range_min[4]),np.std(range_min[4]) )

# print(np.where(Chi2s<2e3*Chi2_min), len(np.where(Chi2s<2e3*Chi2_min)[0]))
# print(np.where(Chi2s<2e4*Chi2_min), len(np.where(Chi2s<2e4*Chi2_min)[0]))