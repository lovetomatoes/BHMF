from PYmodule import *
# produce LF(z=6) with MF & parameterized P(λ)
# MF, LF Willot 2010 checked right
# 但reproduce 的LF 和Willott 2010 有一点偏差 不管什么mass range for BHMF都不行
# Willott breaking point: M1450_star = -25; MBH_star = 2e9
#                         Phi1450_star = 1.14e-8; PhiM_star = 1.23e-8
# 我认为这个Phi1450_star 和PhiM_star 不match; 应该 2.5*Phi1450_star = PhiM_star

z = int(6)
tz = t_from_z(z)
z6datapre = '../z6/data2/'
# MF setting; Willot+ 2010
for M in [2.24e9]:
    print('MBH_star=%.1e'%M)
    print('M1450_star=%.1f'%M1450_Lbol(L_M(M,.6)))
print(.5*(special.erfc(-np.inf) - special.erfc(np.inf)))
print('Phi_ Lfrom MBH_star:', LF(L_M(M,.6)))
M_1 = -27.2
M_2 = -20.7
t = (M1450_Lbol(L_M(M,.6))-M_1) / (M_2-M_1)
t = 2*(1-t) + 3*t
print(t)
print('PhiM_star L from MBH_star 一圈后:', 2.5*(LF(L_M(M,.6)))/.75*t)

# l = np.logspace(45,48)
# plt.figure(figsize=(10,10),dpi=400)
# plt.plot(l/1e44,LF(l))
# plt.xscale('log'); plt.yscale('log')
# plt.xlim(10,1e4); plt.ylim(1e-10,1e-6)
# plt.savefig('../LF_L_W10.png')
# print(LF(L_M(M,.6)))
# exit(0)
abin_mf =  np.logspace(7,10,num=100)
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
M_left = abin_mf[:-1]; M_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]); #print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1
dn_MBH = MF(M_BH)*dlog10M # number density per M_BH bin

# plt.figure(figsize=(10,10),dpi=400)
# plt.plot(M_BH,dn_MBH/dlog10M)
# plt.xscale('log'); plt.yscale('log')
# plt.xlim(1e7,1e10); plt.ylim(1e-10,1e-4)
# plt.grid(True)
# plt.savefig('../MF_W10.png')

# LF data bins same w/ Matsu18
abin_lf = np.linspace(-28,-24,num=30)
abin_lf = np.linspace(-30,-18,num=30) # P_1450 conserve checked
dmag = abin_lf[1]-abin_lf[0] #; print(dmag)
L_left = abin_lf[:-1]; L_right = abin_lf[1:]
M1450  = (L_left+L_right)/2.
N_lf = len(M1450)

f_range = [.75] # 只影响normalization
m_range = [.6]
s_range = [.3]

# f_range = np.arange(.1, 1., .1)
# m_range = np.append( [.01,.05], np.arange(.1,2.,.1))
# s_range = np.arange(.1,1.,.1) # a>0 total P convergent
# print(len(f_range)*len(m_range)*len(s_range) )
# exit(0)

# obscured correction W10
M_1 = -27.2
M_2 = -20.7

i = 0
Chi2_min = 1e10; find_min = False
for f_duty in f_range:
    for mu_fit in m_range:
        for sigma_fit in s_range:
            i = i+1
        # # --------- Luminosity Function ---------
            Phi = np.zeros(N_lf)
            Phi_csv = 0.
            for ibin in range(N_lf):
                # #----------- Schechter lbd -----------
                x0 = kernel_M1450(L_right[ibin], M_BH, mu_fit, sigma_fit)
                x1 = kernel_M1450(L_left[ibin],  M_BH, mu_fit, sigma_fit)
                dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))

                dPhi = np.nansum(dn_MBH*dP_M1450)
                Phi_csv += dPhi
                t = (M1450[ibin]-M_1) / (M_2-M_1)
                t = 2*(1-t) + 3*t
                Phi[ibin] = dPhi/dmag*f_duty/t

            print('consv of dP_M1450:',Phi_csv/np.sum(dn_MBH))

            Phi_Matsu = LF_M1450(M1450) # both Phi & Phi_Matsu in Mpc^-3 mag^-1
            Chi2 = np.nansum(pow(np.log(Phi) - np.log(Phi_Matsu),2))/N_lf

            T = Table(
                [M1450,Phi_Matsu,Phi,Phi_csv/np.sum(dn_MBH)*np.ones(N_lf),Chi2*np.ones(N_lf)],
                names=('M1450','Phi_Matsu','Phi','P1450_csv','Chi2')
            )
            # MF mass range; f_duty, mu, sigma
            LFname = z6datapre+'L_ME_LN_'+ \
                    'M%.2e'%np.log10(np.min(M_BH))+ \
                    'f%.2f'%f_duty+ \
                    'm%.1f'%mu_fit+ \
                    's%.3f'%sigma_fit
            ascii.write(T, LFname,
                        formats={'M1450':'6.2f','Phi_Matsu':'4.2e','Phi':'4.2e','P1450_csv':'4.2e','Chi2':'4.2e'},
                        overwrite=True)

            if np.nanmin([Chi2, Chi2_min]) == Chi2:
                find_min = True
                Chi2_min = Chi2
                LFname_min = LFname

T = ascii.read(LFname_min,guess=False, delimiter=' ')
plt.figure(figsize=(10,10),dpi=400)
plt.plot(abin_lf,1e9*LF_M1450(abin_lf,z=6,W10=True),label='W10 LF')
# plt.plot(T['M1450'], 1e9*T['Phi_Matsu'], '--', label='Matsu')
plt.plot(T['M1450'], 1e9*T['Phi'], '--',label='MF+P'+r'${(\lambda)}$')
plt.yscale('log')
plt.xlim(-21,-29)
plt.ylim(1e-2,1e2)
plt.grid(True)
plt.xlabel('M1450',fontsize=fslabel);plt.ylabel(r'$\mathrm{\Phi \~Gpc^{-3}mag^{-1}}$',fontsize=fslabel)
plt.legend(loc='best',fontsize=fslabel)
plt.savefig('../LF_W10.png')

print(i)
if find_min:
    print(LFname_min,Chi2_min)