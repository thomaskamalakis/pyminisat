from satnet import satsim
from matplotlib import rc
import matplotlib.pyplot as plt

rc('text', usetex=True)
rc('font', size=14)
rc('lines',linewidth=2.0)

s = satsim()
s.simulate()
s.plot_rmin(close_all=True)
s.plot_results()
plt.show()
