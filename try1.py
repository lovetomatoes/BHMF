from PYmodule import *
# N1 = 1; N2 = 2; N3 = 3; N4 = 4; N5 = 5
# a = np.ones((N1,N2,N3,N4,N5))

print((t_from_z(6)-t_from_z(15))/Myr)
print((t_from_z(6)-t_from_z(5))/Myr)
print((t_from_z(5)-t_from_z(4))/Myr)
print(M1450_Lbol(1e13*Lsun))
print(Lsun)

print('z: 50 to 17.58    : %3.2f Myr', (t_from_z(17.58)-t_from_z(50))/Myr)
print('z: 17.58 to 6     : %3.2f Myr', (t_from_z(6)-t_from_z(17.58))/Myr)
print('z: 6 to 4    : %3.2f Myr', (t_from_z(4)-t_from_z(6))/Myr)

# i = 0
# b = []
# for i1 in range(N1):
#     for i2 in range(N2):
#         for i3 in range(N3):
#             for i4 in range(N4):
#                 for i5 in range(N5):
#                     i += 1
#                     b.append(i-1)
#                     # print('i%(N2*N3);',(i-1)%(N3*N2), 'i2=',i2)
#                     # print('i%(N5*N4*N2*N3);',(i-1)%(N5*N4*N3*N2*N1), 'i1=',i1)
#                     # print((i-1)%(N5), 'i5=',i5)
#                     # print((i-1)//N5%N4, 'i4=',i4)
#                     # print((i-1)//N5//N4%N3, 'i3=',i3)
#                     # print((i-1)//N5//N4//N3%N2, 'i2=',i2)
#                     print((i-1)//N5//N4//N3//N2%N1, 'i1=',i1)


N1 = 2; N2 = 3; N3 = 4
a = np.ones((N1,N2,N3))

i = 0
b = []
N_tot = N1*N2*N3
p1 = np.zeros(N_tot); p2 = np.zeros(N_tot); p3 = np.zeros(N_tot)
for i1 in range(N1):
    for i2 in range(N2):
        for i3 in range(N3):
            p1[i] = i1
            p2[i] = i2
            p3[i] = i3
            b.append(i)
            i += 1

print(a.shape)
T = Table([p1,p2,p3,b],names=['p1','p2','p3','b'])
a = np.ones(4)
b[:] = a
# b = a
# b = a.copy()
a += 1
print(b)
print('1e9, 0.1 Eddington ratio: M1450=%.1f'%M1450_Lbol(L_M(1e9,.1)))
print('M1450=-29, edd=0.1, mass:%.1e'%M_L(Lbol_M1450(-29.),.1))
print('M1450=-29, edd=1., mass:%.1e'%M_L(Lbol_M1450(-29.),1.))
print('M1450=-29, edd=10, mass:%.1e'%M_L(Lbol_M1450(-29.),10.))

a = special.gamma([0.5, 1, 5])
a = special.gammainc(1,1)
print(a,1-1./np.exp(1))
# ascii.write(T,'../p1p2p3b',formats={'p1':'5.1e','p2':'5.1e','p3':'5.1e','b':'5.1e'},overwrite=True)