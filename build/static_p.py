# FPRA- NSE
# authors: Max Bräuer and Jannes Wulff
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm
import numpy as np
import re
import glob

# LaTeX style output
plt.rc('text', usetex=True)
plt.rc('font', family='sans',size= 20)

# sort the files naturally
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# Read written data in output
filelist =glob.glob("output/*.txt")
dre = re.compile(r'(\d+)')
filelist.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
pfilelist=filelist[0:int(len(filelist)/3)]

# Extract important information out of .txt files
fn=len(pfilelist)
num = int(np.sqrt(len(np.loadtxt(pfilelist[0],delimiter="\n"))))
pressure = np.zeros((num,num))
pressure =np.loadtxt(pfilelist[1],delimiter="\n").reshape((num,num))
# Some important plot objects
fig = plt.figure(1, [10, 10])
ax = fig.gca()
# Colormap + Colorbar
cmap = mpl.cm.get_cmap('magma')
norm = BoundaryNorm(np.linspace(0, 10, 11), cmap.N);
presplot = ax.imshow(pressure, norm = norm, origin = "lower", cmap=cmap)
bar = fig.colorbar(presplot)
presplot.set_data(pressure)
"""

# check with analytical solution
pressure = np.zeros((num,num))
x = np.linspace(-0.01, 1.01, 102)
y = np.linspace(-0.01, 1.01, 102)
X, Y = np.meshgrid(x, y)
pressure=np.sin(np.pi*Y)*np.sinh(np.pi*X)
presplot.set_data(pressure)

"""
plt.xlabel(r'$x$ Richtung')
plt.ylabel(r"$y$ Richtung")
plt.title(r"Druckprofil im Inneren der Box.")
# plt.show()
plt.savefig("1,0_sol.pdf")
