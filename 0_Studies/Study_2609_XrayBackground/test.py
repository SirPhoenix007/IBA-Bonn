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

length_spectrum = 250

x = np.linspace(0,length_spectrum,length_spectrum)
H = -40*x**2 + 0.001*x**4

a, b, c = 1,1,1
plt.plot(x,H)
for i in range(length_spectrum):
    H_i = H[i]
    y = []
    for j in range(length_spectrum):
        factor = 1 - ((i-j)**2/(a**2*(b+c*i)))
        y.append(H_i*factor)
    plt.plot(x,y, color='black', lw=0.2)

plt.ylim(np.min(H)-5 , np.max(H)+5)
plt.show()
