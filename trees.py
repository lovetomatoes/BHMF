# plot inflow gas Mdot for trees 
from PYmodule import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N_tree = 100
f_b = 0.16

for i in range(N_tree):
    treefile = '/Users/wli/treefiles1e11/tree_'+str(i)+'mer'
    T = ascii.read(treefile, guess=False, delimiter=' ')
    t = T['t_Myr']
    Mh = T['Mh_Ms']
    M_dot = f_b*(Mh[1:]-Mh[:-1])/(t[1:]-t[:-1])/1e6
    z = (T['z'][1:]+T['z'][:-1])/2
    plt.plot(z,M_dot)
plt.yscale('log')
plt.savefig('../Minflow_trees.png')

exit(0)
def linear_z(xa, ya, m, x):
    # linear interpolation
    for i in range(m):
        ms = i
        if x-xa[i]>=0.: # !!!!!!!! decreasing x as redshift
            break
    if (ms==0):
        ms=1; #xa[1],xa[0] extrapolate leftwards
    y1 = ya[ms-1]
    y2 = ya[ms]
    t=(x-xa[ms-1])/(xa[ms]-xa[ms-1])
    y=(1.-t)*y1+t*y2

    # log-log linear interpolation
    logxa = np.log(xa)
    for i in range(m):
        ms = i
        if np.log(x)-logxa[i]>=0.: # !!!!!!!! decreasing x as redshift
            break
    if (ms==0):
        ms=1; #xa[1],xa[0] extrapolate leftwards
    y1 = np.log(ya[ms-1])
    y2 = np.log(ya[ms])
    t=(np.log(x)-logxa[ms-1])/(logxa[ms]-logxa[ms-1])
    y=(1.-t)*y1+t*y2
    y = np.exp(y)
    return y

# #-------------------- median Mh for given one z 
def Mh_median_z(one_z):
    Mhs = []
    for itr in range(10):
        treefile = '../treefiles/tree_'+str(itr)+'mer'
        T = ascii.read(treefile,guess=False)
        Mhs.append(linear_z(T['z'],T['Mh_Ms'],len(T),one_z))
    Mhs = np.array(Mhs)
    return np.median(Mhs)


plt.figure(figsize=(10,8),dpi=400)
for itr in range(N_tree):
    treefile = '../treefiles/tree_'+str(itr)+'mer'
    T = ascii.read(treefile,guess=False)
    plt.plot(T['z'],T['Mh_Ms'],c='grey',label='_',zorder=1,alpha=.5,rasterized=True)

zs = np.linspace(5,50,num=20,endpoint=True)
plt.plot(zs, Mh_Tv(1.e3,zs),':',  c='black',label='_')
plt.plot(zs, Mh_Tv(1.e4,zs),':',  c='black',label='_')
plt.plot(zs, Mh_Tv(1.e5,zs),':',  c='black',label='_')
plt.plot(zs, Mh_Tv(1.e6,zs),':',  c='black',label='_')
plt.text(10,1.e5,s=r'$T_{\rm vir} = 10^3 \rm K$',fontsize=fstxt)
plt.text(10,5.e6,s=r'$T_{\rm vir} = 10^4 \rm K$',fontsize=fstxt)
plt.text(40,3.e8,s=r'$T_{\rm vir} = 10^5 \rm K$',fontsize=fstxt)
plt.text(40,9.e9,s=r'$T_{\rm vir} = 10^6 \rm K$',fontsize=fstxt)

zs = np.arange(50,5,step=-1)
Mh_medians=[]
for z in zs:
    Mh_medians.append(Mh_median_z(z))
Mh_medians = np.array(Mh_medians)

index =  np.where(Mh_medians>1.e5)
plt.plot(zs[index], Mh_medians[index],'--',linewidth=3,c='C3',zorder=2, label='median')
plt.xlabel(r'$\rm z$',fontsize=fslabel);
plt.ylabel(r'$\rm M_h(M_\odot)$',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.yscale('log')
plt.legend(loc='best',fontsize=fslegend)
plt.savefig('./figs/trees'+str(N_tree)+'median.pdf')