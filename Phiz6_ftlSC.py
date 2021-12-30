from PYmodule import *
# Phiz6 paras: f_duty, t_life, lbd ~Schechter(l_cut,a)

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.


Ts = [] # bsm=0,1 two files
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        T['Lbol_z'] = np.zeros(len(T))
        T['M1450_z'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['t_col'] = t_from_z(T['z_col'])
        TM.append(T)
        # print(np.min(T['z_col']),np.max(T['z_col']),np.min(T['Mstar0']),np.max(T['Mstar0']))
    Ts.append(TM)

# print(M1M0(np.array([1e2,1e5]),tz-t_from_z(20),1,.1,.1,0.001))
# print(M1M0(np.array([1e2,1e5]),50*Myr,1,.1,.1,0.001))

# f_duty=.7;mu_fit=.21*30;sigma_fit=.15
# M0 = 1e4; dt = 500*Myr
# print('M1=%.1e'%(M0*np.exp(mu_fit*f_duty*dt/t_Edd)), 'N_mf=%d'%N_mf)

# MF
abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# LF bins same w/ Matsu18
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
N_lf = len(bin_cen)

t_range = np.array([100, 200, 500, 1000])*Myr
f_range = np.arange(.1, 1., .3)
l_range = np.logspace(-2, 1., num=4) # l_cut
a_range = np.array([.1, 1., 2., 3.]) # a>0 total P convergent
# print(len(t_range)*len(f_range)*len(l_range)*len(a_range) )
# print(a_range); exit(0)

t_range = [1000.*Myr]
f_range = [1.]
l_range = [.1]
a_range = [.4]

i = 0
for t_life in t_range:
    for f_duty in f_range:
        for l_cut in l_range:
            for a in a_range:
                # # LF P_tot
                bin_edg = np.linspace(-50,10); N_lf = len(bin_edg)-1
                # bin_left = bin_edg[:-1]; bin_right = bin_edg[1:]
                # N_lf = len(bin_left)
                # M_BH = np.logspace(5,12,num=400)
                # Phi = np.zeros(N_lf)
                # # for ibin in range(N_lf):
                # #     x0 = kernelS_M1450(bin_edg[ibin+1], M_BH, l_cut)
                # #     x1 = kernelS_M1450(bin_edg[ibin],   M_BH, l_cut)
                # #     Phi[ibin] = np.nansum(special.gammainc(a,x1) - special.gammainc(a,x0))
                # # print('sum of dP_M1450:',np.nansum(Phi)/len(M_BH))
                # M_BH = 1e6
                # x0 = kernelS_M1450(bin_right, M_BH, l_cut)
                # x1 = kernelS_M1450(bin_left,  M_BH, l_cut)
                # Phi = np.nansum(special.gammainc(a,x1) - special.gammainc(a,x0))
                # print('sum of dP_M1450:',np.nansum(Phi))

                i = i+1
                # continue
                dn_MBH = np.zeros(N_mf)
                for iM in range(N_Mh):
                    for i_bsm in range(Nbsm):
                        T = Ts[iM][i_bsm]
                        Nt = np.max((tz-T['t_col'])//t_life)
                        Nmax = Nt
                        dP_MBH = np.zeros(N_mf)
                        while Nt>=0:
                            t_point = tz - Nt*t_life
                            T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
                            dt_seed = t_point - T_seed['t_col']
                            dP_MBH_prev = dP_MBH.copy()
                            # M_BHt = M1M0(M_BH,dt,f_duty,mu_fit,eta8,delta_fit) # if not exp growth, may needed
                            for ibin in range(N_mf):
                                # new seeds 
                                if len(T_seed):
                                    # #----------- log-norm lbd -----------
                                    # x0 = kernel_MBH1(bin_left[ibin] /T_seed['Mstar0'],dt_seed,f_duty, mu_fit, sigma_fit)
                                    # x1 = kernel_MBH1(bin_right[ibin]/T_seed['Mstar0'],dt_seed,f_duty, mu_fit, sigma_fit)
                                    # # x0 = ma.masked_invalid(x0); x1 = ma.masked_invalid(x1) # 必删! or not conserved
                                    # dP_seed = -.5*(special.erfc(x1) - special.erfc(x0)) # P(<x) = 1-.5*erfc(x)
                                    # #----------- Schechter lbd -----------
                                    x0 = kernelS_MBH(bin_left[ibin] /T_seed['Mstar0'], dt_seed, f_duty, l_cut)
                                    x1 = kernelS_MBH(bin_right[ibin]/T_seed['Mstar0'], dt_seed, f_duty, l_cut)
                                    x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                    dP_seed = special.gammainc(a,x1) - special.gammainc(a,x0)

                                    dP_seed = np.nansum(dP_seed)/len(T)
                                else:
                                    dP_seed = 0.
                                # prev BHMF
                                # #----------- log-norm lbd -----------
                                # x0 = kernel_MBH1(M_BH[ibin]/bin_right,t_life,f_duty, mu_fit, sigma_fit)
                                # x1 = kernel_MBH1(M_BH[ibin]/bin_left, t_life,f_duty, mu_fit, sigma_fit)
                                # dP_MBH[ibin] = np.nansum(.5*(special.erfc(x0) - special.erfc(x1)) * dP_MBH_prev) + dP_seed
                                # #----------- Schechter lbd -----------
                                x0 = kernelS_MBH(M_BH[ibin]/bin_right, t_life, f_duty, l_cut)
                                x1 = kernelS_MBH(M_BH[ibin]/bin_left,  t_life, f_duty, l_cut)
                                x0[x0<0] = 0.; x1[x1<0] = 0. # let P(growth_ratio<1)=0, must! or not conserved!
                                dP_MBH[ibin] = np.nansum((special.gammainc(a,x1) - special.gammainc(a,x0)) * dP_MBH_prev) + dP_seed
                            Nt -= 1
                        dn_MBH += dP_MBH*n_base[iM]*f_bsm[i_bsm]

                T = Table(
                    [M_BH, dn_MBH],
                    names=('M_BH','dn_MBH')
                )
                MFname = z6datapre+'MF_SC_'+'t%.1e'%(t_life/Myr)+ \
                        'f%.1f'%f_duty+ \
                        'l%.1e'%l_cut+ \
                        'a%.3f'%a+ \
                        'alpha%.1f'%alpha
                ascii.write( Table([T['M_BH'], T['dn_MBH']/dlog10M],
                            names=['M_BH','dn_dlog10M']),
                            MFname,
                            formats={'M_BH':'4.2e','dn_dlog10M':'4.2e'},
                            overwrite=True)
                consv_ratio = np.nansum(dn_MBH)/np.sum(n_base)
                print('conserved fraction=%.10f'%consv_ratio)
                if consv_ratio<.9:
                    print('conserved fraction=%.10f'%consv_ratio)
                # exit(0)
                T  = T[np.logical_and(True,T['M_BH']<2e10)] # select M_BH range

                Phi = np.zeros(N_lf)
                Phi_csv = 0.
                M_BH = T['M_BH']
                for ibin in range(N_lf):
                    #----------- log-norm lbd -----------
                    mu_fit = .21; sigma_fit = .15
                    x0 = kernel_M1450(bin_edg[ibin+1], M_BH, mu_fit, sigma_fit)
                    x1 = kernel_M1450(bin_edg[ibin], M_BH, mu_fit, sigma_fit)
                    dP_M1450 = .5*(special.erfc(x0) - special.erfc(x1))
                    #----------- Schechter lbd -----------
                    x0 = kernelS_M1450(bin_edg[ibin+1], M_BH, l_cut)
                    x1 = kernelS_M1450(bin_edg[ibin],   M_BH, l_cut)
                    dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

                    dPhi = np.nansum(T['dn_MBH']*dP_M1450)
                    Phi_csv += dPhi
                    # Phi[ibin] = dPhi/bin_wid[ibin]*f_duty
                print('consv of dP_M1450:',Phi_csv/np.sum(n_base))
                exit(0)

                T = Table(
                    [bin_cen, Phi*1e9, Phi*1e9*(1.-f_obsc_const), Phi*1e9/corr_U14D20(bin_cen), Phi_obs],
                    names=('bin_cen','Phi','Phi_CO','Phi_DO','Phi_obs')
                )
                LFname = z6datapre+'LF_SC_'+'t%.1e'%(t_life/Myr)+ \
                        'f%.1f'%f_duty+ \
                        'l%.1e'%l_cut+ \
                        'a%.3f'%a+ \
                        'alpha%.1f'%alpha
                ascii.write(T, LFname,
                            formats={'bin_cen':'6.2f','Phi':'4.2e','Phi_CO':'4.2e','Phi_DO':'4.2e','Phi_obs':'4.2e'},
                            overwrite=True)
                # exit(0)
print(i)