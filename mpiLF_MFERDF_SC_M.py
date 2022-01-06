from PYmodule import *
import sys 
sys.path.append('/gpfs/share/home/1801110214')
from mpi4py import MPI
# produce LF(z=6) with MF & parameterized P(λ); λ has mass dependency

z = int(6)
tz = t_from_z(z)
z6datapre = '../z6/data2/'
# MF setting; Willot+ 2010
abin_mf =  np.logspace(2,12,num=100)
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
M_left = abin_mf[:-1]; M_right = abin_mf[1:]
dlog10M = np.log10(abin_mf[1]/abin_mf[0]); print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1
dn_MBH = MF(M_BH)*dlog10M # number density per M_BH bin

# LF data bins same w/ Matsu18
abin_lf = np.linspace(-29,-18,num=30)
dmag = abin_lf[1]-abin_lf[0] #; print(dmag)
L_left = abin_lf[:-1]; L_right = abin_lf[1:]
M1450  = (L_left+L_right)/2.
N_lf = len(M1450)

f_range = np.arange(.1, 1.1, .1)
d_range = np.arange(.1, 1.1, .1)
d_range = [.1]
l_range = np.append( [.01,.05], np.arange(.1,2.,.1))
a_range = np.append( [.01,.05], np.arange(.1, 3., .1)) # a>0 total P convergent
print(len(l_range)*len(a_range) )

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
assert len(f_range)*len(d_range) == size

i = 0
Chi2_min = 1e10; find_min = False; LFname_min = ''

f_duty = f_range[ rank%len(d_range) ]
d_fit = d_range[ rank%len(d_range) ]

for l_cut in l_range:
    for a in a_range:
        i = i+1
        # continue
    # # --------- Luminosity Function ---------
        Phi = np.zeros(N_lf)
        Phi_csv = 0.
        for ibin in range(N_lf):
            # #----------- Schechter lbd -----------
            x0 = kernelS_M1450_M(L_right[ibin], M_BH, l_cut, d_fit)
            x1 = kernelS_M1450_M(L_left[ibin],  M_BH, l_cut, d_fit)
            dP_M1450 = special.gammainc(a,x1) - special.gammainc(a,x0)

            dPhi = np.nansum(dn_MBH*dP_M1450)
            Phi_csv += dPhi
            Phi[ibin] = dPhi/dmag*f_duty
            Phi_Matsu = LF_M1450(M1450) # both Phi & Phi_Matsu in Mpc^-3 mag^-1
            Chi2 = np.nansum(pow(np.log(Phi) - np.log(Phi_Matsu),2))/N_lf

            T = Table(
                [M1450,Phi_Matsu,Phi,Chi2*np.ones(N_lf)],
                names=('M1450','Phi_Matsu','Phi','Chi2')
            )
            LFname = z6datapre+'L_ME_SC_M_'+ \
                    'f%.1f'%f_duty+ \
                    'l%.1e'%l_cut+ \
                    'a%.3f'%a
            # ascii.write(T, LFname,
            #             formats={'M1450':'6.2f','Phi_Matsu':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            #             overwrite=True)
            if np.nanmin([Chi2, Chi2_min]) == Chi2:
                    find_min = True
                    Chi2_min = Chi2
                    LFname_min = LFname

LFnames = np.array(comm.gather(LFname_min, root=0))
Chi2_mins = np.array(comm.gather(Chi2_min, root=0))
if rank == 0:
    print('Chi2_mins',Chi2_mins)
    print('Chi2_min=\n',np.nanmin(Chi2_mins), 'file:\n',LFnames[np.nanargmin(Chi2_mins)])