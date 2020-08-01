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


# Re 10000
#u_ghia = [1.00000, 0.47221, 0.47783, 0.48070, 0.47804, 0.34635, 0.20673, 0.08344, 0.03111, -0.07540, -0.23186, -0.32709, -0.38000, -0.41657, -0.42537, -0.42735, 0.00000]
#y_ghia = [1.00000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]

#v_ghia = [0.00000, -0.54302, -0.52987, -0.49099, -0.45863, -0.41496, -0.36737, -0.30719, 0.00831, 0.27224, 0.28003, 0.35070, 0.41487, 0.43124, 0.43733, 0.43983, 0.0000]
#x_ghia = [1.000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8074, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]

# Re 1000
#u_ghia = [1.0000, 0.47221, 0.47783, 0.48070, 0.4780, 0.34635, 0.20673, 0.08344, 0.03111, -0.07540, -0.23186, -0.32709, -0.38000, -0.41657, -0.42537, -0.42735, 0.00000]
#y_ghia = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.000]

#v_ghia = [0.00000, -0.54302, -0.52987, -0.49099, -0.45863, -0.41496, -0.36737, -0.30719, 0.00831, 0.27224, 0.28003, 0.35070, 0.41487, 0.43124, 0.43733, 0.43983, 0.0000]
#x_ghia = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]

# Re 100
# u_ghia = [1.0000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000]
# y_ghia = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
#
# v_ghia = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000]
# x_ghia = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]



# # LaTeX style output
# plt.rc('text', usetex=True)
# plt.rc('font', family='sans',size=20)

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


# Static Plot of last data
velo_abs=get_data_abs(uvend-1)
uvplot.set_data(velo_abs)
u_res,v_res=get_data_vec(uvend-1)
velo_vec.set_UVC(u_res,v_res)

plt.xlabel(r'$x$ direction')
plt.ylabel(r"$y$ direction")
plt.title(r"Velocities in the box")
plt.savefig("uv_10000.pdf")
plt.show()
