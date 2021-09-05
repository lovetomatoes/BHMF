import os 
from PYmodule import *
for flambda in [0.16,0.17]:
    for len_T in [6,7]:
        z = 4
        sigma_fit = .09
        fname_old = './figs/z'+str(int(z))+'f'+str(int(flambda*100))+'s_fit'+str(int(sigma_fit*100))+'N'+str(int(len_T)-4)+'.png'
        fname_new = './figs/z'+str(int(z))+'f'+str(int(flambda*100))+'s0'+'N'+str(int(len_T)-4)+'.png'
        os.system('mv '+fname_old+' '+fname_new)
