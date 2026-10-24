import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times"],
    "text.usetex": True,
    "font.size": 8,
    "pgf.rcfonts": False
})


plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

a = [1,2,3,4,5]
b = [0,1,2,1,0]

plt.figure(figsize=(4,3), dpi=150)

plt.plot(a,b, label='TEST123')
plt.legend()
plt.tight_layout()
plt.show()
