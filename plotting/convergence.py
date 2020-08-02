# FPRA- NSE
# authors: Max Bräuer and Jannes Wulff
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm
import numpy as np
from numpy import linalg as LA
from skimage.transform import resize

uf1="final_states_and_videos/global_convergence/50x50/u_822.txt"
vf1="final_states_and_videos/global_convergence/50x50/v_822.txt"
uf2="final_states_and_videos/global_convergence/100x100/u_2018.txt"
vf2="final_states_and_videos/global_convergence/100x100/v_2018.txt"
uf3="final_states_and_videos/global_convergence/200x200/u_6072.txt"
vf3="final_states_and_videos/global_convergence/200x200/v_6072.txt"

sz=52
u1= np.loadtxt(uf1,delimiter="\n").reshape((sz,sz))[1:-1,1:-1]
v1= np.loadtxt(vf1,delimiter="\n").reshape((sz,sz))[1:-1,1:-1]
u2= np.loadtxt(uf2,delimiter="\n").reshape((2*sz-2,2*sz-2))[1:-1,1:-1]
v2= np.loadtxt(vf2,delimiter="\n").reshape((2*sz-2,2*sz-2))[1:-1,1:-1]
u3= np.loadtxt(uf3,delimiter="\n").reshape((4*sz-6,4*sz-6))[1:-1,1:-1]
v3= np.loadtxt(vf3,delimiter="\n").reshape((4*sz-6,4*sz-6))[1:-1,1:-1]

u2=u2[::2,::2]
v2=v2[::2,::2]
u3=u3[::4,::4]
v3=v3[::4,::4]
# print(u3-u2)
# print(u2-u1)
# print(u2-u1)
fig = plt.figure(1,[10,10])
ax = fig.gca()
index1=np.arange(0,sz-2,1)
index2=np.arange(0,sz-2,1)
index1,index2=np.meshgrid(index1,index2)
print(u3[index1,index2])
test=np.log2((u3[index1,index2]-u2[index1,index2])/(u2[index1,index2]-u1[index1,index2]))
# n=LA.norm(u2-u1)
# z=LA.norm(u3-u2)
cmap = mpl.cm.get_cmap('magma')
norm = BoundaryNorm(np.linspace(0, 4, 10), cmap.N);
uvplot = ax.imshow(np.abs(test), norm = norm, origin = "lower", cmap=cmap)
bar = fig.colorbar(uvplot)
plt.show()
print(1/50*np.sqrt(np.nansum(np.square(test))))
# print(n)
# print(z)
"""
print(len(u3))
print(len(u2))
print(len(u1))
"""
