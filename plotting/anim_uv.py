# FPRA- NSE
# authors: Max Bräuer and Jannes Wulff
import numpy as np
import re
import glob

# sort the files naturally
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# Read written data in output
filelist =glob.glob("output/*.txt")
dre = re.compile(r'(\d+)')
filelist.sort(key=lambda l: [int(s) if s.isdigit() else s.lower() for s in re.split(dre, l)])
ufilelist=filelist[int(len(filelist)/3):2*int(len(filelist)/3)]
vfilelist=filelist[2*int(len(filelist)/3)::]

# Extract important information out of .txt files
fn=len(ufilelist)
num = int(np.sqrt(len(np.loadtxt(ufilelist[0],delimiter="\n"))))
velo_abs = np.zeros((num,num))
