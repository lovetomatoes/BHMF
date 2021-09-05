from PYmodule import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

log10Ms = [9,10,11,12,13]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']


print('z:35, t_Hubble: ', t_from_z(35)/Myr)
print('z:20, t_Hubble: ', t_from_z(25)/Myr)
print('t_edd: ', t_Edd/Myr)

flambda = .2
print('efolding growth ratio: ', (t_from_z(25)-t_from_z(35)) / t_Edd * flambda )
print('growth ratio: ', np.exp((t_from_z(25)-t_from_z(35)) / t_Edd * flambda) )

Tzs = [ [[],[]], [[],[]], [[],[]], [[],[]], [[],[]]]
Ts_iso = [] # n_tell = 1.e4, T_tell = 4000K
Ts_H2 = []
Ts_isofail = []
Ts_isoOK = []
pres = ['../data/1e9','../data/1e10','../data/1e11','../data/1e12','../data/1e13']

zs = np.arange(5,51,5)

Nz = len(zs)-1
NM = 5
Nbsm = 2

Tzs = [ [[0 for i in range(NM)] for j in range(Nbsm)] for k in range(Nz)]

T_tell = 8000
print(t_Edd/Myr)

eta = 0.3


for iM in range(NM):
    for i_bsm in range(Nbsm):
        if iM ==2 and i_bsm ==1:
            T=ascii.read(pres[iM]+'Jcol_'+str(i_bsm)+'alpha1.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        else:
            T=ascii.read(pres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar'] = np.zeros(len(T))
        print('log10Ms[iM],i_bsm,mean z_col',log10Ms[iM],i_bsm,np.mean(T['z_col']))
        for iz in range(Nz):
            Tzs[iz][i_bsm][iM] = T[np.logical_and(T['z_col']>=zs[iz], T['z_col']<zs[iz+1])]
        #  #------------- print check selected table correct ----------------
        #     if Tzs[iz][i_bsm][iM]:
        #         print('iz',zs[iz],np.max(Tzs[iz][i_bsm][iM]['z_col']),np.min(Tzs[iz][i_bsm][iM]['z_col']),'ibsm',i_bsm,np.mean(Tzs[iz][i_bsm][iM]['i_bsm']),
        #         'iM',iM,log10Ms[iM] )
        for i in range(len(T)):
            T['Mstar'][i] = Mdot2M(T['Mdot'][i]*eta)
        print('M* range and mean',np.min(T['Mstar']),np.max(T['Mstar']), np.mean(T['Mstar']))

abin = np.logspace(np.log10(1e2),np.log10(1e6),num=30)  # endpoint=True
wid = abin[1:]-abin[:-1]
dlog10M = np.log10(abin[1]/abin[0])
pre = ''

# # ----------- seperate by z_col ---------------
# for iz in range(Nz):
#     h_0 = np.zeros(len(abin)-1)
#     h_1 = np.zeros(len(abin)-1)
#     h_2 = np.zeros(len(abin)-1)
#     for i_bsm in range(Nbsm):
#         h00 = np.zeros(len(abin)-1)
#         h11 = np.zeros(len(abin)-1)
#         h22 = np.zeros(len(abin)-1)
#         for iM in range(NM):
#             T = Tzs[iz][i_bsm][iM]
#             for i in range(len(T)):
#                 T['Mstar'][i] = Mdot2M(T['Mdot'][i]*eta)
#             T_H2 = T[T['Tg_max']<=T_tell] 
#             T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
#             T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
#             # if T_isoOK:
#             #     print(np.min(T_isoOK['Mstar']), np.max(T_isoOK['Mstar']))
#             hist0, bin_edges = np.histogram(T_H2['Mstar'],bins=abin,density=False)
#             hist1, bin_edges = np.histogram(T_isofail['Mstar'],bins=abin,density=False)
#             hist2, bin_edges = np.histogram(T_isoOK['Mstar'],bins=abin,density=False)
#             h00 = h00 + hist0*n_base[iM]/1e4
#             h11 = h11 + hist1*n_base[iM]/1e4
#             h22 = h22 + hist2*n_base[iM]/1e4
#         h_0 = h_0 + h00*f_bsm[i_bsm]
#         h_1 = h_1 + h11*f_bsm[i_bsm]
#         h_2 = h_2 + h22*f_bsm[i_bsm]
#     # whole+= np.sum(h_0+h_1+h_2)

#     plt.figure(figsize=(10,8),dpi=400)
#     plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_0/dlog10M,width=wid,color='C'+str(0),alpha=1,label=pre+typenames[0])
#     plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_1/dlog10M,width=wid,bottom=h_0/dlog10M,color='C'+str(1),alpha=1,label=pre+typenames[1])
#     plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_2/dlog10M,width=wid,bottom=(h_0+h_1)/dlog10M,color='C'+str(2),alpha=1,label=pre+typenames[2])
#     plt.tick_params(labelsize=fstick)
#     # s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
#     plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
#     plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
#     plt.xscale('log'); plt.yscale('log')
#     plt.title(r'$\eta=$'+'0.'+str(int(10*eta))+'  z='+str(int(zs[iz]))+r'$\to$'+str(int(zs[iz]+5)),fontsize=fslabel)
#     plt.xlim(1e2,1e6); plt.ylim(1.e-4,10)
#     if i_bsm==0:
#         plt.legend(fontsize=fslegend,loc='best')
#     plt.savefig(figpre+'eta'+str(int(10*eta))+'z'+str(int(zs[iz]))+'.png')


# h_0 = np.zeros(len(abin)-1)
# h_1 = np.zeros(len(abin)-1)
# h_2 = np.zeros(len(abin)-1)
# for iz in range(Nz):
#     for i_bsm in range(Nbsm):
#         h00 = np.zeros(len(abin)-1)
#         h11 = np.zeros(len(abin)-1)
#         h22 = np.zeros(len(abin)-1)
#         for iM in range(NM):
#             T = Tzs[iz][i_bsm][iM]
#             for i in range(len(T)):
#                 T['Mstar'][i] = Mdot2M(T['Mdot'][i]*eta)
#             T_H2 = T[T['Tg_max']<=T_tell] 
#             T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
#             T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
#             # if T_isoOK:
#             #     print(np.min(T_isoOK['Mstar']), np.max(T_isoOK['Mstar']))
#             hist0, bin_edges = np.histogram(T_H2['Mstar'],bins=abin,density=False)
#             hist1, bin_edges = np.histogram(T_isofail['Mstar'],bins=abin,density=False)
#             hist2, bin_edges = np.histogram(T_isoOK['Mstar'],bins=abin,density=False)
#             h00 = h00 + hist0*n_base[iM]/1e4
#             h11 = h11 + hist1*n_base[iM]/1e4
#             h22 = h22 + hist2*n_base[iM]/1e4
#         h_0 = h_0 + h00*f_bsm[i_bsm]
#         h_1 = h_1 + h11*f_bsm[i_bsm]
#         h_2 = h_2 + h22*f_bsm[i_bsm]
# plt.figure(figsize=(10,8),dpi=400)
# plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_0,width=wid,color='C'+str(0),alpha=1,label=pre+typenames[0])
# plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_1,width=wid,bottom=h_0,color='C'+str(1),alpha=1,label=pre+typenames[1])
# plt.bar( np.sqrt(abin[:-1]*abin[1:]),h_2,width=wid,bottom=(h_0+h_1),color='C'+str(2),alpha=1,label=pre+typenames[2])
# plt.tick_params(labelsize=fstick)
# s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
# plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
# plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
# plt.xscale('log'); plt.yscale('log')
# plt.title(r'$\eta=$'+'0.'+str(int(10*eta))+'  z='+str(int(zs[iz]))+r'$\to$'+str(int(zs[iz]+5)),fontsize=fslabel)
# plt.xlim(1e2,1e6); plt.ylim(1.e-6,10)
# plt.legend(fontsize=fslegend,loc='best')
# plt.savefig(figpre+'eta'+str(int(10*eta))+'allz.png')


plt.figure(figsize=(10,8),dpi=400)

# ----------- seperate by Mbase ---------------
for iM in [2]:
    # continue
    print(iM)
    h_0 = np.zeros(len(abin)-1)
    h_1 = np.zeros(len(abin)-1)
    h_2 = np.zeros(len(abin)-1)


    for i_bsm in range(Nbsm):
        h00 = np.zeros(len(abin)-1)
        h11 = np.zeros(len(abin)-1)
        h22 = np.zeros(len(abin)-1)
        h33 = np.zeros(len(abin)-1)
        for iz in range(Nz):
            T = Tzs[iz][i_bsm][iM]
            for i in range(len(T)):
                T['Mstar'][i] = Mdot2M(T['Mdot'][i]*eta)
            T_H2 = T[T['Tg_max']<=T_tell] 
            T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
            T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
            # if T_isoOK:
            #     print(np.min(T_isoOK['Mstar']), np.max(T_isoOK['Mstar']))
            hist0, bin_edges = np.histogram(T_H2['Mstar'],bins=abin,density=False)
            hist1, bin_edges = np.histogram(T_isofail['Mstar'],bins=abin,density=False)
            hist2, bin_edges = np.histogram(T_isoOK['Mstar'],bins=abin,density=False)
            h00 = h00 + hist0
            h11 = h11 + hist1
            h22 = h22 + hist2
        # h_0 = h_0 + h00 #*f_bsm[i_bsm]
        # h_1 = h_1 + h11 #*f_bsm[i_bsm]
        # h_2 = h_2 + h22 #*f_bsm[i_bsm]
        
        sign = 1 if i_bsm==1 else -1
        plt.bar( np.sqrt(abin[:-1]*abin[1:]),sign*h00/1.e4,color=None,edgecolor='black',width=wid,alpha=0.3,label=pre+typenames[0]+str(i_bsm))
        plt.bar( np.sqrt(abin[:-1]*abin[1:]),sign*h11/1.e4,color=None,edgecolor='black',width=wid,bottom=sign*h_0/1e4,alpha=0.3,label=pre+typenames[1]+str(i_bsm))
        plt.bar( np.sqrt(abin[:-1]*abin[1:]),sign*h22/1.e4,color=None,edgecolor='black',width=wid,bottom=sign*(h_0+h_1)/1e4,alpha=0.3,label=pre+typenames[2]+str(i_bsm))
        print('sign ibsm and sum of histograms: ',sign, np.sum(h00+h11+h22))

    plt.tick_params(labelsize=fstick)
    s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel('fraction',fontsize=fslabel)
    plt.xscale('log'); plt.yscale('symlog')
    plt.title(r'$\eta=0.3$',fontsize=fslabel)
    plt.xlim(1e2,1e6); plt.ylim(-.2,.2)
    if i_bsm==0:
        plt.legend(fontsize=fslegend,loc='best')

    plt.savefig(figpre+'eta3M'+str(int(log10Ms[iM]))+'bsm1_2side.png')