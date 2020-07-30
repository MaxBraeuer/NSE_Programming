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

# LaTeX style output
plt.rc('text', usetex=True)
plt.rc('font', family='sans',size=20)

# for mp4 output
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=30, metadata=dict(artist='MB'), bitrate=1800)

# sort the files naturally
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# Read written data in output
filelist =glob.glob("output/*.txt")
dre = re.compile(r'(\d+)')
filelist.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
uvend=int(len(filelist)/3)
ufilelist=filelist[int(len(filelist)/3):2*int(len(filelist)/3)]
vfilelist=filelist[2*int(len(filelist)/3)::]

# Extract important information out of .txt files
fn=len(ufilelist)
num = int(np.sqrt(len(np.loadtxt(ufilelist[0],delimiter="\n"))))
velo_abs = np.zeros((num,num))

# Some important plot objects
fig = plt.figure(1, [10, 10])
ax = fig.gca()
# Colormap + Colorbar
cmap = mpl.cm.get_cmap('magma')
norm = BoundaryNorm(np.linspace(0, 1, 10), cmap.N);
uvplot = ax.imshow(velo_abs, norm = norm, origin = "lower", cmap=cmap)
bar = fig.colorbar(uvplot)
X, Y = np.mgrid[0:100:10j, 0:100:10j]
print(X)
velo_vec = ax.quiver(Y, X, Y, X, scale=2.0, color="white")

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
def get_data_abs(i):
    u=np.loadtxt(ufilelist[i],delimiter="\n").reshape((num,num))
    v=np.loadtxt(vfilelist[i],delimiter="\n").reshape((num,num))
    return np.sqrt(np.power(u,2)+np.power(v,2))
def get_data_vec(i):
    print(i)
    u=np.loadtxt(ufilelist[i],delimiter="\n").reshape((num,num))
    v=np.loadtxt(vfilelist[i],delimiter="\n").reshape((num,num))
    return resize(u, (10, 10)), resize(v, (10, 10))
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
velo_abs=get_data_abs(uvend-1)
uvplot.set_data(velo_abs)
u_res,v_res=get_data_vec(uvend-1)
velo_vec.set_UVC(u_res,v_res)
"""
# Animation Plot
ani = animation.FuncAnimation(
    fig, animate,frames=fn, init_func=init, interval=10, blit=True)
"""
plt.xlabel(r'$x$ direction')
plt.ylabel(r"$y$ direction")
plt.title(r"Velocities in the box")
plt.savefig("uv_10000.pdf")
plt.show()
# To save the clip
# ani.save('im.mp4', writer=writer)
