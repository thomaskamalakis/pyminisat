from satnet import satsim
from matplotlib import rc
import matplotlib.pyplot as plt

rc('text', usetex=True)
rc('font', size=14)
rc('lines',linewidth=2.0)

s = satsim()
s.optimize_F_mp()

s.plot_rmin(close_all=True)

orb_i=1
sat_i=2
orb_j=2
s.calc_grid_dists()
s.calc_dist_power()
plt.close('all')
subset1 = ['U0', 'U1', 'U2', 'S1']
strs1 = {
    'U0' : '-ko',
    'U1' : '-s',
    'U2' : '-d',
    'S1' : '--r'
    }
subset2 = ['L0', 'L1', 'L2', 'S2']
strs2 = {
    'L0' : '-ko',
    'L1' : '-s',
    'L2' : '-d',
    'S2' : '--r'}
s.plot_dists(subset=subset1,
             plot_strs=strs1)

s.plot_dists(subset=subset2,
             plot_strs=strs2)

s.plot_cg(subset=subset1,
             plot_strs=strs1)

s.plot_cg(subset=subset2,
             plot_strs=strs2)

s.plot_PS(subset=subset1,
             plot_strs=strs1)

s.plot_PS(subset=subset2,
             plot_strs=strs2)

s.plot_Preq(subset=subset1,
             plot_strs=strs1)

s.plot_Preq(subset=subset2,
             plot_strs=strs2)
