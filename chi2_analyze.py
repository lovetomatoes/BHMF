from PYmodule import *
# 对于任意参数组合 有z=4的chi2 (我的定义 sum( pow( (Phi-Phi_obs)/sigma, 2) 
# 读取+分析 chi2 output文件

z = int(4)

bin_cen = bin_cen[str(z)]
bin_wid = bin_wid[str(z)]
bin_edg = bin_edg[str(z)]
# ------------------- select magnitude range
M1450_min = -30.; M1450_max = -22. # z=4
Phi_obs = ma.masked_where(np.logical_or(bin_cen>M1450_max,bin_cen<M1450_min),Phi_obs[str(z)])
Phi_err = Phi_err[str(z)]

alpha = 1.


# T = Table([ds,es,fs,ms,ss,Chi2s],names=['ds','es','fs','ms','ss','Chi2s'])
fname = z4datapre+'Chi2_para'
T = ascii.read(fname, guess=False, delimiter=' ')
ds=T['ds']; es=T['es']; fs=T['fs']; ms=T['ms']; ss=T['ss']
Chi2s = T['Chi2s']
i_min = np.argmin(Chi2s)
Chi2_min = Chi2s[i_min]

range_min =  np.where(Chi2s<3*Chi2_min); print('len range_min',len(ds[range_min]))
# print(ds[range_min],es[range_min],fs[range_min],ms[range_min],ss[range_min])


print(ds[i_min],'%.3f'%np.mean(ds[range_min]),'%.3f'%np.median(ds[range_min]),'%.3f'%(np.std(ds[range_min])) )
print(es[i_min],'%.3f'%np.mean(es[range_min]),'%.3f'%np.median(es[range_min]),'%.3f'%(np.std(es[range_min])) )
print(fs[i_min],'%.2f'%np.mean(fs[range_min]),'%.2f'%np.median(fs[range_min]),'%.2f'%(np.std(fs[range_min])) )
print(ms[i_min],'%.2f'%np.mean(ms[range_min]),'%.2f'%np.median(ms[range_min]),'%.2f'%(np.std(ms[range_min])) )
print(ss[i_min],'%.2f'%np.mean(ss[range_min]),'%.2f'%np.median(ss[range_min]),'%.2f'%(np.std(ss[range_min])) )