from PYmodule import *
from scipy.optimize import curve_fit

t1 = time.time()

f_seed = 0.01
f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
prex = '../4p/M0r8_' + f_seedlabel + 'lin' # on sk1
# prex = '../M0r8_' + f_seedlabel + 'lin' # local

Mmin = 1e9
outMprex = prex +'Mmin%d'%(int(np.log10(Mmin)))
Lmax = -26
outLprex = prex +'Lmax%d'%abs(int(Lmax))


figm, axm = plt.subplots(figsize=(10, 10))

rhoM = []; rhoM_min = []; rhoM_max = []
nL = []; nL_min = []; nL_max = []
# labels = {'NIRCam_deep','NIRCam_med','Roman_deep','Roman_wide','Euclid_deep','Euclid_wide'}
styles = ['-','--','-.','.']
colors = ['k','red','orange']
marks = ['o','v','P','D','p','+','s']
markcs = ['b','g','c','m','y']
markcs = ['grey']*10
labels = ['LSST_deep','Roman_wide','Euclid_wide','Roman_deep','Euclid_deep']
figlabs = ['LSST deep','Roman wide','Euclid wide']
N_detect = {'LSST_deep':[], 'Roman_wide':[], 'Euclid_wide':[]}

zs = np.arange(7,11)
# zs = [7]
for z in zs:
    figl, axl = plt.subplots(figsize=(10, 10))
    fMname = prex+'BHMFz{:d}'.format(z)
    fLname = prex+'QLFz{:d}'.format(z)
    TM = ascii.read(fMname, guess=False, delimiter=' ')
    TL = ascii.read(fLname, guess=False, delimiter=' ')

#   plot best and spread for BHMF/QLF at z
    y_best, med, spread = TM['y_best'], TM['med'], TM['spread']
    axm.fill_between(TM['Mx'],med-spread,med+spread,color='C%d'%9,alpha=0.3,label=r'_$1\sigma$ Posterior Spread',zorder=3)
    axm.plot(TM['Mx'], y_best,linewidth=lw, color='C%d'%9, label='z=%d'%z,zorder=3)
    y_best, med, spread = TL['y_best'], TL['med'], TL['spread']
    y_noobs = y_best * corr_U14D20(M1450)
    axl.fill_between(TL['M1450'],med-spread,med+spread,color='C%d'%9,alpha=0.3,label=r'_$1\sigma$ Posterior Spread',zorder=3)
    axl.plot(TL['M1450'], y_best, linewidth=lw,color='C%d'%9, label='z=%d'%z,zorder=3)
    axl.plot(TL['M1450'], y_noobs, '--',linewidth=lw,color='C%d'%9, label='z=%d'%z,zorder=3)
    # z=9, add Galaxy LF curve (data dots)
    if z==9:
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='Schechter'), '--',c='purple')
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='DPL'), c='purple')
        obslabels = ['Donnan','Harikane','Bowler','Stefanon','Morishita']
        obslabel = 'Donnan'; iobs = obslabels.index(obslabel)
        T_data = ascii.read('../z_9_LF_data/z9_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e3,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Harikane'; iobs = obslabels.index(obslabel)
        T_data = ascii.read('../z_9_LF_data/z9_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e4, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['dPhi_-']*1e4,T_data['dPhi_+']*1e4),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Bowler'; iobs = obslabels.index(obslabel)
        T_data = ascii.read('../z_9_LF_data/z9_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e3,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Stefanon'; iobs = obslabels.index(obslabel)
        T_data = ascii.read('../z_9_LF_data/z9_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag'],
            alpha=.5,yerr = (T_data['Phi_-']*1e3,T_data['Phi_+']*1e3),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Morishita'; iobs = obslabels.index(obslabel)
        T_data = ascii.read('../z_9_LF_data/z9_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], pow(10,T_data['logPhi'])*1e9, xerr = T_data['dmag']/2,
            alpha=.5,
            yerr = ( pow(10,T_data['logPhi'])*1e9 - pow(10,T_data['logPhi']-T_data['logPhi_-'])*1e9,
            pow(10,T_data['logPhi']+T_data['logPhi_+'])*1e9 - pow(10,T_data['logPhi'])*1e9 ),
            fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
    
    if z==8:
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='Schechter'), '--',c='purple')
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='DPL'), c='purple')
        obslabels = ['Donnan','Bouwens','McLure','Bowler','Stefanon'] 
        obslabel = 'Donnan'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_8_LF_data/z8_'+obslabel+'.dat',guess=False,delimiter=' ',)
        # T_data = ascii.read('../z_8_LF_data/z8_Donnan.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e3,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Bouwens'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_8_LF_data/z8_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e9, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e9,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'McLure'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_8_LF_data/z8_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e9, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e9,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Bowler'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_8_LF_data/z8_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e3,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Stefanon'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_8_LF_data/z8_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['Phi_-']*1e3,T_data['Phi_+']*1e3),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
       
    if z==7:
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='Schechter'), '--',c='purple')
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='DPL'), c='purple')
        obslabels = ['Bouwens','Harikane'] 
        obslabel = 'Bouwens'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z-7_LF_data/z7_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e9, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e9,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'Harikane'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z-7_LF_data/z7_'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['Phi_-']*1e3,T_data['Phi_+']*1e3),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
   
    if z==10:
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='Schechter'), '--',c='purple')
        axl.plot(TL['M1450'],LF_Gal(TL['M1450'],z, mod='DPL'), c='purple')
        obslabels = ['z10-12_Donnan','z10-12_Harikane','z13_Harikane','z13_Naidu','z10_Morishita','z11_GN-z11']
        obslabel = 'z10-12_Donnan'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = T_data['dPhi']*1e3,fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'z10-12_Harikane'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e4, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['Phi_-']*1e4,T_data['Phi_+']*1e4),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'z13_Harikane'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e1, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['Phi_-']*1e1,T_data['Phi_+']*1e1),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'z13_Naidu'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], pow(10,T_data['logPhi'])*1e9, xerr = T_data['dmag']/2,
            alpha=.5,yerr = ( pow(10,T_data['logPhi'])*1e9 - pow(10,T_data['logPhi']-T_data['dPhi_-'])*1e9,
            pow(10,T_data['logPhi']+T_data['dPhi_+'])*1e9 - pow(10,T_data['logPhi'])*1e9 ),
            fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'z10_Morishita'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], pow(10,T_data['logPhi'])*1e9, xerr = T_data['dmag']/2,
            alpha=.5,yerr = ( pow(10,T_data['logPhi'])*1e9 - pow(10,T_data['logPhi']-T_data['logPhi_-'])*1e9,
            pow(10,T_data['logPhi']+T_data['logPhi_+'])*1e9 - pow(10,T_data['logPhi'])*1e9 ),
            fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)
        obslabel = 'z11_GN-z11'; iobs = obslabels.index(obslabel); print(obslabel)
        T_data = ascii.read('../z_10-13_LF_data/'+obslabel+'.dat',guess=False,delimiter=' ',)
        axl.errorbar(T_data['M_uv'], T_data['Phi']*1e3, xerr = T_data['dmag']/2,
            alpha=.5,yerr = (T_data['Phi_-']*1e3,T_data['Phi_+']*1e3),fmt=marks[iobs],ms=7.5,color=markcs[iobs],elinewidth=2,capsize=4)




    # # QLF add detection limit
    # area same, Roman_deep deeper than Euclid_deep
    for ilabel in range(3):
        label = labels[ilabel]
        deep_mag = M_absolute(Depth[label],z)
        Vc_cover = Vc(Area[label],z,1)/1e9
        detection = [np.linspace(-40,deep_mag,num=50),np.linspace(1./Vc_cover,1e10,num=50)]
        # horizontal
        axl.plot(detection[0],1./Vc_cover*np.ones(50),styles[0],linewidth=1,color=colors[ilabel],zorder=4)
        # vertical
        axl.plot(deep_mag*np.ones(50),detection[1],styles[0],linewidth=1,color=colors[ilabel],zorder=4)
        if ilabel ==0:
            axl.text(deep_mag-.2, 1e6, 'LSST Deep', fontsize=18,color=colors[0],zorder=4)
        elif ilabel ==1:
            axl.text(deep_mag-.2, 3e5, 'Roman Wide', fontsize=18,color=colors[1],zorder=4)
        else:
            axl.text(deep_mag-.2, 1e6, 'Euclid Wide', fontsize=18,color=colors[2],zorder=4)
        if ilabel==0 and z==6:
            print(label+' depth mag={:.1f}, Vc={:.1e}, 1/Vc={:.1e}'.format(deep_mag,Vc_cover,1./Vc_cover))
            print('faintest abundance: ',np.max(y_best[TL['M1450']<deep_mag]))
            print('integrate QLF: ',np.nansum(y_best[TL['M1450']<deep_mag]*dmag))
            print('Vc_cover: ',Vc_cover)
        N_detect[label].append(np.nansum(y_best[TL['M1450']<deep_mag]*dmag*Vc_cover))

    axl.set_xlim(np.max(TL['M1450']),np.min(TL['M1450']))
    axl.set_ylim(5e-3,1e2)
    axl.set_xlim(-18,-30)
    axl.set_ylim(1e-3,1e3)
    # axl.grid(1)
    from matplotlib.ticker import MultipleLocator
    ml = MultipleLocator(10)
    # if z==9:

    xleft = -18; xright = -30
    logybottom = -6; logytop = 7
    axl.set_xlim(-18,-30)
    axl.set_ylim(10**logybottom,10**logytop)

    axl.set_yscale('log')
    
    ## set y ticks
    y_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = logytop-logybottom+1)
    axl.yaxis.set_major_locator(y_major)
    y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    axl.yaxis.set_minor_locator(y_minor)
    axl.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    ## set x ticks
    axl.set_xticks(np.arange(xleft,xright-1,-1))

    # axl.set_xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
    # axl.set_ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
    axl.tick_params(labelsize=fstick)
    # plt.tight_layout()
    figl.savefig(prex+'LF_spreadz%d.png'%z,dpi=300,bbox_inches='tight')

axm.set_xlim(1e7,1e10); axm.set_xscale('log')
axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fstick)
figm.savefig(prex+'MF_spread.png',dpi=300,bbox_inches='tight')

# minimum integrant, may be 0
print('rhoM_min',rhoM_min)
print('nL_min',nL_min)

print('N_detect:',N_detect)