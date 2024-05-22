from linkmodel import link, receiver, transmitter
import matplotlib.pyplot as plt
from constellation import constellation
import numpy as np
import pickle
from matplotlib import rc
import os
DEFAULTS = {
    'LMdB' : 3,
    'alt' : 550e3,
    'Nsat' : 66,
    'Norb' : 24,
    'F' : 13,
    'i' : 53 * np.pi / 180,
    'target_BER' : 1e-12,
    'l' : 1.55e-6,
    'eta' : 0.8,
    'redB' : 20,
    'Rb' : 10e9,
    'RL' : 100,
    'FndB' : 3,
    'TK' : 300,
    'DR' : 80e-3,
    'nR' : 0.8,
    'sR' : 1e-6,
    'sT' : 1e-6,
    'theta_T' : None,
    'PT' : None,
    'nT' : 0.8,
    'Nt' : 10000,
    'filename' : 'sat.data',    # Default filename to save data
    'sat' : 1,                  # Satellite assumed for F optimizations and minimal distance estimations
    'orb' : 1,                   # Orbit assumed for F optimizations and minimal distance estimations
    'no_cores' : os.cpu_count()
    }

CONSTELLATION_P = ['alt','Nsat','Norb','F','i', 'orb', 'sat']
RECEIVER_P = ['target_BER','l','eta','redB','Rb','RL','FndB','TK','DR','nR','sR']
TRANSMITTER_P = ['sT','theta_T','PT','nT']


class satsim:
    
    def set_param(self, k, params):
        v = params.get(k, DEFAULTS[k])
        setattr(self, k, v)
        
    def init_graph(self):
        rc('text', usetex=True)
        rc('font', size=16)
        rc('lines',linewidth=2.0)
    def __init__(self, **kwargs):
        
        # initialize transmitter
        params_t = {k: kwargs.get(k, DEFAULTS[k]) for k in TRANSMITTER_P}       
        self.transmitter = transmitter( **params_t )
        
        # initialize receiver
        params_r = {k: kwargs.get(k, DEFAULTS[k]) for k in RECEIVER_P}       
        self.receiver = receiver( **params_r )
        
        # initialize constellation
        params_c = {k: kwargs.get(k, DEFAULTS[k]) for k in CONSTELLATION_P}       
        self.constellation = constellation( **params_c )
        
        # time axis
        print(kwargs)

        self.set_param('Nt', kwargs)
        self.t = np.linspace(0, self.constellation.T, self.Nt)
        
        # Link margin
        self.set_param('LMdB', kwargs)

        # Satellite assumed for optimizations
        self.set_param('sat', kwargs)
        self.set_param('orb', kwargs)
        
        # Filename used to save data
        self.set_param('filename', kwargs)
        
        # Number of CPUs to be used in multiprocessing
        self.set_param('no_cores', kwargs)
        
        # Initialize anomalies
        self.constellation.calc_init_anomalies()
        
        # Calculate receiver sensitivity
        self.receiver.calc_Psens()
                    
    def optimize_F(self, plot_result = False):
        self.constellation.optimize_F(orb = self.orb, sat = self.sat)

    def optimize_F_mp(self, plot_result = False):
        self.constellation.optimize_F_mp(no_cores=self.no_cores)

    def save(self):
            
        with open(self.filename, 'wb') as f:
            pickle.dump(self, f)
    
    def plot_rmin(self, close_all=False):
        if close_all:
            plt.close('all')
            
        plt.figure()
        plt.plot(self.constellation.rmin/1e3, '-o')
        plt.xlabel('$F$')
        plt.ylabel(r'$r_\mathrm{min}\;\mathrm{[Km]}$')
        plt.grid()
    
    def calc_grid_dists(self):
        self.dists = self.constellation.grid_dists()
    
    def calc_dist_power(self):

        self.link = link(tx = self.transmitter,
                         rx = self.receiver,
                         LMdB = self.LMdB
                         )
        self.link.rx.calc_Psens()        
        self.PT_req_dBm = {}
        self.LPS_dB = {}
        self.L_dB = {}
        
        for k in self.dists:            
            self.PT_req_dBm[k] = np.zeros(self.t.shape)
            self.LPS_dB[k] = np.zeros(self.t.shape)
            self.L_dB[k] = np.zeros(self.t.shape)
            
            for i, t in enumerate(self.t):
                self.link.set_R( self.dists[k][i] )
                self.link.calc_required_PT()
                self.PT_req_dBm[k][i] = self.link.PT_req_dBm
                self.LPS_dB[k][i] = self.link.LPSdB
                self.L_dB[k][i] = self.link.LdB
                
            
    def simulate(self):
        self.calc_grid_dists()
        
    def plot_dists(self, subset = None, plot_strs=None,markevery=200):
        if not subset:
            subset = self.constellation.sat_grid(self.orb, self.sat).keys()
        
        fig = plt.figure()
        ax = fig.gca()
        
        for k in subset:
            if plot_strs:
                ax.plot(self.t, self.dists[k]/1e3, plot_strs[k], label = k, markevery=markevery)
            else:
                ax.plot(self.t, self.dists[k]/1e3, label = k, markevery=markevery)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid()
        
        plt.ylabel(r'$r_\mathrm{min}\;\mathrm{[Km]}$')
        plt.xlabel(r'$t\;\mathrm{[s]}$')
        
        plt.tight_layout()
      
    def plot_cg(self,subset = None, plot_strs=None,markevery=200):
        if not subset:
            subset = self.constellation.sat_grid(self.orb, self.sat).keys()
        
        plt.figure()
        for k in subset:
            if plot_strs:
                plt.plot(self.t, self.L_dB[k], plot_strs[k], label = k, markevery=markevery)
            else:
                plt.plot(self.t, self.L_dB[k], label = k, markevery=markevery)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid()
        
        plt.ylabel(r'$L_\mathrm{ch}\;\mathrm{[dB]}$')
        plt.xlabel(r'$t\;\mathrm{[s]}$')
        
        plt.tight_layout()
      
    def plot_PS(self,subset = None, plot_strs=None,markevery=200):
        if not subset:
            subset = self.constellation.sat_grid(self.orb, self.sat).keys()
        
        plt.figure()
        for k in subset:
            if plot_strs:
                plt.plot(self.t, self.LPS_dB[k], plot_strs[k], label = k, markevery=markevery)
            else:
                plt.plot(self.t, self.LPS_dB[k], label = k, markevery=markevery)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid()
        
        plt.ylabel(r'$L_\mathrm{PS}\;\mathrm{[dB]}$')
        plt.xlabel(r'$t\;\mathrm{[s]}$')
        
        plt.tight_layout()
        
    def plot_Preq(self,subset = None, plot_strs=None,markevery=200):
        if not subset:
            subset = self.constellation.sat_grid(self.orb, self.sat).keys()
        
        plt.figure()
        for k in subset:
            if plot_strs:
                plt.plot(self.t, self.PT_req_dBm[k], plot_strs[k], label = k, markevery=markevery)
            else:
                plt.plot(self.t, self.PT_req_dBm[k], label = k, markevery=markevery)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.grid()
        
        plt.ylabel(r'$P_\mathrm{T}\;\mathrm{[dBm]}$')
        plt.xlabel(r'$t\;\mathrm{[s]}$')
        
        plt.tight_layout()
        
        
    def __repr__(self):
        return self.constellation.__repr__() + '\n' + \
               self.transmitter.__repr__() + '\n' +  \
               self.receiver.__repr__() + '\n'
               
    def __str__(self):
        return self.__repr__()
               
     