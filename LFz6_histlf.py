from PYmodule import *
# calibrating LFz6.py

print('z:45, t_Hubble: %3.2f Myr', t_from_z(45)/Myr)
print('z:30, t_Hubble: %3.2f Myr', t_from_z(30)/Myr)
print('z:10, t_Hubble: %3.2f Myr', t_from_z(10)/Myr)
print('z:5,  t_Hubble: %3.2f Myr', t_from_z(5)/Myr)
print('z:32 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(32))/Myr)
print('t_Edd         : %3.2f Myr', t_Edd/Myr)

NM = 3

T_tell = 8000
eta = 0.3

# for num in np.linspace(1,100):
#     print(num,math.erfc(num))
# exit(0)

T = ascii.read('../z6/data/MF')
# T  = T[np.logical_and(T['bin_cen']>1e6,T['bin_cen']<1e12)]
dlog10M = np.log10(T['bin_cen'][1]/T['bin_cen'][0])

# LF test bins
bin_obs = np.linspace(-29,-22,50)
bin_wid = bin_obs[1:]-bin_obs[:-1]
bin_cen = bin_obs[:-1]+bin_wid/2.
Nbin = len(bin_cen)

for M1450 in [-29., -22.]:
    print('M1450',M1450,'Mbh = %3.2e'%M_L(Lbol_M1450(M1450),1.))

flambda = .19 # .18 .20
z = 6

# for i in np.linspace(1,50):
#     print('i',i,' ',math.erfc(i))
# exit(0)

f_duty = .5
sigma_fit = .12
mu_fit = flambda/f_duty
Phi = np.zeros(Nbin)
Nsamp = 1000
for i in range(len(T)):
    M_BH = T['bin_cen'][i]
    dn_MBH = T['Phi'][i]*dlog10M
    Edd_ratio_now = lognorm.rvs(sigma_fit*np.log(10), scale=mu_fit, size=Nsamp) # 使用rvs随机生成lbd 未用本身的Edd_ratio
    Lbol_z = L_M(M_BH, Edd_ratio_now)
    M1450_z = M1450_Lbol(Lbol_z)
    hist_lf, bin_edges = np.histogram(M1450_z,bins=bin_obs,density=False)
    Phi += hist_lf/Nsamp/bin_wid*dn_MBH*f_duty

T = Table(
    [bin_cen, Phi],
    names=('bin_cen','Phi')
)
# ascii.write(T, z6datapre+'allmassMF_fl'+str(int(flambda*100))+'f'+str(int(f_duty*10))+'s'+str(int(sigma_fit*100))+'bsm'+str(Nbsm-1)+'alpha1'+'NM'+str(NM),formats={'bin_left':'6.2e','Phi':'4.2e'},overwrite=True)
ascii.write(T, z6datapre+'LF_histlf',formats={'bin_cen':'6.2e','Phi':'4.2e'},overwrite=True)