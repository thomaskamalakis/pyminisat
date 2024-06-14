from satnet import satsim
from matplotlib import rc
import matplotlib.pyplot as plt
import time

rc('text', usetex=True)
rc('font', size=14)
rc('lines',linewidth=2.0)

t1 = time.time()
s = satsim()
s.simulate()
t2 = time.time()
print(t2-t1)
s.plot_rmin(close_all=True)
s.plot_results()
plt.show()
