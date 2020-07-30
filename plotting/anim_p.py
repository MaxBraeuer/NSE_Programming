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

# LaTeX style output
plt.rc('text', usetex=True)
plt.rc('font', family='sans',size=20)

# for mp4 output
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='MB'), bitrate=1800)

# sort the files naturally
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# Read written data in output
filelist =glob.glob("output/*.txt")
dre = re.compile(r'(\d+)')
filelist.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
pend = int(len(filelist)/3)
pfilelist=filelist[0:pend]

# Extract important information out of .txt files
fn=len(pfilelist)
num = int(np.sqrt(len(np.loadtxt(pfilelist[0],delimiter="\n"))))
print(num)
pressure = np.zeros((num,num))


# Some important plot objects
fig = plt.figure(1, [10, 10])
ax = fig.gca()
# Colormap + Colorbar
cmap = mpl.cm.get_cmap('magma')
norm = BoundaryNorm(np.linspace(-0.1, 0.1, 10), cmap.N);
presplot = ax.imshow(pressure, norm = norm, origin = "lower", cmap=cmap)
bar = fig.colorbar(presplot)
# get current data points
def get_data(i):
    return np.loadtxt(pfilelist[i],delimiter="\n").reshape((num,num))
# initial plot configuration
def init():
    pressure=get_data(0)
    presplot.set_data(pressure)
    return presplot,
#animation function
def animate(i):
    print(i)
    pressure =get_data(i)
    presplot.set_data(pressure)
    return presplot,

# Static Values in the end
pressure =get_data(pend-1)
presplot.set_data(pressure)

# Animation Plot
"""
ani = animation.FuncAnimation(
    fig, animate,frames=fn, init_func=init, interval=50, blit=True)
"""
plt.xlabel(r'$x$ direction')
plt.ylabel(r"$y$ direction")
plt.title(r"Pressure in the box")
plt.savefig("p_10000.pdf")
plt.show()
# To save the clip
# ani.save('im.mp4', writer=writer)
