from PYmodule import *

dlogx = .0001
logx_min = -2.; logx_max = 1.+dlogx

def P_tot(a):
    logls = np.arange(logx_min,logx_max,dlogx)
    ls = pow(10., logls)
    return np.sum(pow(ls,a) * np.exp(-ls) * dlogx)

def P_left(a,l):
    # integration of dP ~ x^a exp(-x) dlogx; normalized by integral over (-2,1)
    x = pow(10., np.arange(logx_min,logx_max,dlogx))
    xx, ll = np.meshgrid(x, l)
    yy = xx<ll
    zz = yy*np.logical_and(logx_min<=np.log10(ll),np.log10(ll)<=logx_max) * pow(xx,a)*np.exp(-xx) * dlogx
    # zz = (xx<ll)*np.logical_and(logx_min<=np.log10(ll),np.log10(ll)<=logx_max) * pow(xx,a)*np.exp(-xx) * dlogx
    return np.sum(zz,axis=1)

def P_left_norm(a,l):
    return P_left(a,l)/P_tot(a)

""" 
def P_tot(a):
    logx = logx_min
    P_tot = 0
    while logx<1:
        x = pow(10., logx)
        P_tot += pow(x,a)*np.exp(-x) * dlogx
        logx += dlogx
    # print('P_tot',P_tot)
    return P_tot

def P_left(a,l):
    logl = np.log10(l)
    if logx_min <= logl <= logx_max:
        logx = logx_min
        P = 0
        while logx<logl:
            x = pow(10., logx)
            P += pow(x,a)*np.exp(-x) * dlogx
            logx += dlogx
        return P
    else:
        return np.nan
 """