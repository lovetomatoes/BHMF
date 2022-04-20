from PYmodule import *
# seed BHMF for all Mh, v_bsm; migrated from src_prev/Mdist_seed.py
# no big diff: all seed MF compared w/ Mh=1e11,vbsm=0 seed MF (generated from src_prev/Phi_seed.py)
# current plot directly, better write hist file later
log10Ms = [11,12,13]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']


print('z:35, t_Hubble: ', t_from_z(35)/Myr)
print('z:20, t_Hubble: ', t_from_z(25)/Myr)
print('t_edd: ', t_Edd/Myr)

# print(len(Ts[0])); exit(0)
# Ts shape: [3][2]
Ts_iso = [] # n_tell = 1.e4, T_tell = 4000K
Ts_H2 = []
Ts_isofail = []
Ts_isoOK = []
pres = ['../data/1e11','../data/1e12','../data/1e13']
figprefix = '../'

dz = 5
zs = np.arange(5,51,dz)
# [ 5 10 15 20 25 30 35 40 45 50]
dz = 10
zs = np.arange(10,51,dz)
# [10 20 30 40 50]

Nz = len(zs)-1
NM = 3
Nbsm = 2

# check eta=0.6, 0.3 difference
for iM in range(NM):
    for i_bsm in range(Nbsm):
        for eta in [.6,.3]:
            T = Ts[iM][i_bsm]
            for j in range(len(T)):
                T['Mstar0'][j] = Mdot2M(T['Mdot'][j]*eta)
            print('eta={0:.1f} iM={1:d} i_bsm={2:d} ;'.format(eta,iM,i_bsm), 'min, max of Mstar0: {0:.1e}, {1:.1e} Msun'.format(np.min(T['Mstar0']),np.max(T['Mstar0'])))


Tzs = [ [[0 for i in range(NM)] for j in range(Nbsm)] for k in range(Nz)]

whole = 0
for iM in range(NM):
    for i_bsm in range(Nbsm):
        T = Ts[iM][i_bsm]
        for iz in range(Nz):
            Tzs[iz][i_bsm][iM] = T[np.logical_and(zs[iz]<=T['z_col'], T['z_col']<zs[iz+1])]
            # whole += len(Tzs[iz][i_bsm][iM])
        #  #------------- print check selected table correct ----------------
        #     if Tzs[iz][i_bsm][iM]:
        #         print(zs[iz],' to ',zs[iz+1],
        #         np.max(Tzs[iz][i_bsm][iM]['z_col']),np.min(Tzs[iz][i_bsm][iM]['z_col']),
        #         'ibsm',i_bsm,np.mean(Tzs[iz][i_bsm][iM]['i_bsm']),
        #         'iM',iM,log10Ms[iM] )
        # if iM==0 and i_bsm==0:
        #     print(np.max(T['Mstar0'])); exit(0)
# print(whole); exit(0)

T_tell = 8000
x = (bin_left+bin_right)/2.


# ----------- seperate by fbsm ---------------
for i_bsm in range(Nbsm):
    h_0 = np.zeros(N_mf)
    h_1 = np.zeros(N_mf)
    h_2 = np.zeros(N_mf)
    for iM in range(NM):
        T = Ts[iM][i_bsm]

        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
        # if T_isoOK:
        #     print(np.min(T_isoOK['Mstar0']), np.max(T_isoOK['Mstar0']))
        hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
        hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
        hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)

        h_0 += hist0*n_base[iM]; h_1 += hist1*n_base[iM]; h_2 += hist2*n_base[iM]
    h_0 *= f_bsm[i_bsm]; h_2 *= f_bsm[i_bsm]; h_2 *= f_bsm[i_bsm]
    
    plt.figure(figsize=(10,8),dpi=400)
    plt.bar(x,h_0/1.e4,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    plt.bar(x,h_1/1.e4,width=wid_mf,bottom=h_0/1e4,color='C'+str(1),alpha=0.5,label=typenames[1])
    plt.bar(x,h_2/1.e4,width=wid_mf,bottom=(h_0+h_1)/1e4,color='C'+str(2),alpha=0.5,label=typenames[2])
    plt.tick_params(labelsize=fstick)
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    plt.xscale('log'); plt.yscale('log')
    plt.title(r'$\mathrm{v_{bsm}}=$'+str(i_bsm)+r' $\sigma$',fontsize=fslabel)
    plt.xlim(1e2,1e6); plt.ylim(1.e-12,1e-2)
    plt.tight_layout()
    plt.savefig(figprefix+'eta3bsm'+str(i_bsm)+'.png')

# ----------- seperate by z_col ---------------
for iz in range(Nz):
    h = np.zeros(N_mf)
    h_0 = np.zeros(N_mf); h_1 = np.zeros(N_mf); h_2 = np.zeros(N_mf)
    for i_bsm in range(Nbsm):
        hh = np.zeros(N_mf)
        hh_0 = np.zeros(N_mf); hh_1 = np.zeros(N_mf); hh_2 = np.zeros(N_mf)
        for iM in range(NM):
            T = Tzs[iz][i_bsm][iM]
            hist, bin_edges = np.histogram(T['Mstar0'],bins=abin_mf,density=False)
            hh = hh + hist*n_base[iM]/1e4

            T_H2 = T[T['Tg_max']<=T_tell] 
            T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
            T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
            # if T_isoOK:
            #     print('im,ibsm,itr:',iM,i_bsm)
            #     print(np.min(T_isoOK['Mstar0']), np.max(T_isoOK['Mstar0']))
            hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
            hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
            hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)
            hh_0 = hh_0 + hist0*n_base[iM]/1e4
            hh_1 = hh_1 + hist1*n_base[iM]/1e4
            hh_2 = hh_2 + hist2*n_base[iM]/1e4
        h_0 = h_0 + hh_0*f_bsm[i_bsm]
        h_1 = h_1 + hh_1*f_bsm[i_bsm]
        h_2 = h_2 + hh_2*f_bsm[i_bsm]
        h += hh*f_bsm[i_bsm]

    plt.figure(figsize=(10,5),dpi=400)
    plt.bar(x,h/dlog10M,width=wid_mf,color='black',alpha=1,fill=0,label=typenames[0])
    plt.bar(x,h_0/dlog10M,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    plt.bar(x,h_1/dlog10M,width=wid_mf,bottom=h_0/dlog10M,color='C'+str(1),alpha=0.5,label=typenames[1])
    plt.bar(x,h_2/dlog10M,width=wid_mf,bottom=(h_0+h_1)/dlog10M,color='C'+str(2),alpha=0.5,label=typenames[2])
    plt.tick_params(labelsize=fstick)
    # s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    plt.xscale('log'); plt.yscale('log')
    plt.text(1e5,1e-3,r'$\eta=$'+'0.'+str(int(10*eta))+'\nz='+str(int(zs[iz]))+r'$\to$'+str(int(zs[iz]+dz)),fontsize=fslabel)
    plt.xlim(1e2,1e6); plt.ylim(1.e-8,1e-2)
    if i_bsm==0:
        plt.legend(fontsize=fslegend,loc='best')
    plt.tight_layout()
    plt.savefig(figprefix+'eta'+str(int(10*eta))+'z'+str(int(zs[iz]))+'.png')

# ----------- all z_col sample ---------------
h_0 = np.zeros(N_mf)
h_1 = np.zeros(N_mf)
h_2 = np.zeros(N_mf)
for i_bsm in range(Nbsm):
    hh_0 = np.zeros(N_mf)
    hh_1 = np.zeros(N_mf)
    hh_2 = np.zeros(N_mf)
    for iM in range(NM):
        T = Ts[iM][i_bsm]

        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]

        hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
        hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
        hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)
        hh_0 = hh_0 + hist0*n_base[iM]/1e4
        hh_1 = hh_1 + hist1*n_base[iM]/1e4
        hh_2 = hh_2 + hist2*n_base[iM]/1e4
    h_0 = h_0 + hh_0*f_bsm[i_bsm]
    h_1 = h_1 + hh_1*f_bsm[i_bsm]
    h_2 = h_2 + hh_2*f_bsm[i_bsm]
plt.figure(figsize=(10,8),dpi=400)
plt.bar(x,h_0,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
plt.bar(x,h_1,width=wid_mf,bottom=h_0,color='C'+str(1),alpha=0.5,label=typenames[1])
plt.bar(x,h_2,width=wid_mf,bottom=(h_0+h_1),color='C'+str(2),alpha=0.5,label=typenames[2])
plt.tick_params(labelsize=fstick)
s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
plt.xscale('log'); plt.yscale('log')
plt.title(r'$\eta=$'+'0.'+str(int(10*eta)),fontsize=fslabel)
plt.xlim(1e2,1e6); plt.ylim(1.e-8,1e-2)
plt.legend(fontsize=fslegend,loc='best')
plt.tight_layout()
plt.savefig(figprefix+'eta'+str(int(10*eta))+'allz.png')
# all z,Mh,vbsm seed BHMF
ascii.write(Table([x,(h_0+h_1+h_2)/dlog10M]),'../IMFhist.dat',names=['M','Phi'],
formats={'M':'10.2e','Phi':'10.2e'},overwrite=True)

# ----------- seperate by Mbase ---------------
for iM in range(NM):
    h_0 = np.zeros(N_mf)
    h_1 = np.zeros(N_mf)
    h_2 = np.zeros(N_mf)
    for i_bsm in range(Nbsm):
        # hh_0 = np.zeros(N_mf)
        # hh_1 = np.zeros(N_mf)
        # hh_2 = np.zeros(N_mf)
        
        T = Ts[iM][i_bsm]

        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
        # if T_isoOK:
        #     print(np.min(T_isoOK['Mstar0']), np.max(T_isoOK['Mstar0']))
        hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
        hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
        hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)
        hh_0 = hh_0 + hist0*f_bsm[i_bsm]
        hh_1 = hh_1 + hist1*f_bsm[i_bsm]
        hh_2 = hh_2 + hist2*f_bsm[i_bsm]
        h_0 += hist0*f_bsm[i_bsm]
        h_1 += hist1*f_bsm[i_bsm]
        h_2 += hist2*f_bsm[i_bsm]

    h_0 *= n_base[iM]
    h_1 *= n_base[iM]
    h_2 *= n_base[iM]
    # h_2 = h_2 + hh_2*n_base[iM]
    plt.figure(figsize=(10,8),dpi=400)
    plt.bar(x,h_0/1.e4,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    plt.bar(x,h_1/1.e4,width=wid_mf,bottom=h_0/1e4,color='C'+str(1),alpha=0.5,label=typenames[1])
    plt.bar(x,h_2/1.e4,width=wid_mf,bottom=(h_0+h_1)/1e4,color='C'+str(2),alpha=0.5,label=typenames[2])
    plt.tick_params(labelsize=fstick)
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    plt.xscale('log'); plt.yscale('log')
    plt.title(r'$\mathrm{M_h=1e}$'+str(int(log10Ms[iM])),fontsize=fslabel)
    plt.xlim(1e2,1e6); plt.ylim(1.e-12,1e-2)
    plt.tight_layout()
    plt.savefig(figprefix+'eta3M'+str(int(log10Ms[iM]))+'.png')