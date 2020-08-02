# FPRA- NSE
# authors: Max Bräuer and Jannes Wulff
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm
import numpy as np
import re
import glob
from skimage.transform import resize
def get_data_abs(i):
    u=np.loadtxt(ufilelist[i],delimiter="\n").reshape((num,num))
    v=np.loadtxt(vfilelist[i],delimiter="\n").reshape((num,num))
    return np.sqrt(np.power(u,2)+np.power(v,2))
def get_data_vec(i):
    # print(i)
    u=np.loadtxt(ufilelist[i],delimiter="\n").reshape((num,num))
    v=np.loadtxt(vfilelist[i],delimiter="\n").reshape((num,num))
    return u, v

# LaTeX style output
# plt.rc('text', usetex=True)
plt.rc('font', family='sans',size=20)

# for mp4 output
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=30, metadata=dict(artist='MB'), bitrate=1800)

# sort the files naturally
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)
"""
# Read written data in output
filelist =glob.glob("output/*.txt")
dre = re.compile(r'(\d+)')
filelist.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
uvend=int(len(filelist)/3)
ufilelist=filelist[int(len(filelist)/3):2*int(len(filelist)/3)]
vfilelist=filelist[2*int(len(filelist)/3)::]
"""
ufilelist=["final_states_and_videos/Re10-10000_final_uv_data/Re10000_250x250/u_7878.txt"]
vfilelist=["final_states_and_videos/Re10-10000_final_uv_data/Re10000_250x250/v_7878.txt"]
uvend=1
# Extract important information out of .txt files
fn=len(ufilelist)
num = int(np.sqrt(len(np.loadtxt(ufilelist[0],delimiter="\n"))))
velo_abs = np.zeros((num,num))
# border_u_x=25
# border_l_x=5
# border_u_y=236
# border_l_y=216

border_u_y=240
border_l_y=210
border_u_x=30
border_l_x=0


velo_abs=get_data_abs(uvend-1)
# uvplot.set_data(velo_abs[border_l:border_u-1,border_l:border_u-1])
u_res,v_res=get_data_vec(uvend-1)
# velo_vec.set_UVC(u_res[border_l:border_u-1,border_l:border_u-1],v_res[border_l:border_u-1,border_l:border_u-1])
# print(u_res[0:20-1,0:20-1])
# Some important plot objects
fig = plt.figure(1,[10,10])
ax = fig.gca(xlim=(border_l_x,border_u_x),ylim=(border_l_y,border_u_y))
# Colormap + Colorbar
factor=2
sz=int(border_u_x-border_l_x)/factor
X = np.arange(border_l_x, border_u_x, factor)
Y =  np.arange(border_l_y, border_u_y, factor)
U, V = np.meshgrid(X, Y)
udata=resize(u_res[border_l_y:border_u_y,border_l_x:border_u_x],(sz,sz))
vdata=resize(v_res[border_l_y:border_u_y,border_l_x:border_u_x],(sz,sz))
# velo_vec = ax.quiver(Y,X, X, Y, scale=2.0, color="white")
# velo_vec.set_UVC(udata,vdata)
cmap = mpl.cm.get_cmap('magma')
norm = BoundaryNorm(np.linspace(0, 0.2, 20), cmap.N);
uvplot = ax.imshow(velo_abs[border_l_y:border_u_y,border_l_x:border_u_x], norm = norm, origin = "lower", cmap=cmap,extent=[border_l_x,border_u_x,border_l_y,border_u_y])
bar = fig.colorbar(uvplot)

velo_vec = ax.quiver(X,Y,udata, vdata, scale=0.5,color="white")

# # TEST:get current data points
# def get_data_abs_test(i):
#     u=np.power(X/100,2)*np.power(Y/100,2)
#     v=u*np.sin(np.pi/10*i)
#     return np.sqrt(np.power(u,2)+np.power(v,2))
# def get_data_vec_test(i):
#     u=np.power(X/100,2)*np.power(X/100,2)*np.power(Y/100,2)
#     v=u*np.sin(np.pi/10*i)
#     return resize(u, (10, 10)), resize(v, (10, 10))

# get current data points
# initial plot configuration
def init():
    velo_abs=get_data_abs(0)
    uvplot.set_data(velo_abs)
    u_res,v_res=get_data_vec(0)
    velo_vec.set_UVC(u_res,v_res)
    return uvplot, velo_vec
#animation function
def animate(i):
    velo_abs =get_data_abs(i)
    uvplot.set_data(velo_abs)
    u_res,v_res=get_data_vec(i)
    velo_vec.set_UVC(u_res,v_res)

    return uvplot, velo_vec

# Static Plot of last data
"""
# Animation Plot
ani = animation.FuncAnimation(
    fig, animate,frames=fn, init_func=init, interval=10, blit=True)
"""
# plt.ylim(border_l_x,border_u_x)
# plt.xlim(border_l_y,border_u_y)

plt.xlabel(r'$x$ direction')
plt.ylabel(r"$y$ direction")
plt.title(r"Velocities in the box")
plt.savefig("eddi_ol.pdf")
plt.show()
# To save the clip
# ani.save('im.mp4', writer=writer)
