#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:50:57 2024

@author: thomas
"""
from constants import constants
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
from multiprocessing import Pool
import os

class Arrow3D(FancyArrowPatch):
    """
    A class used to draw a 3D arrow in the constellation visualizations
    """
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)

setattr(Axes3D, 'arrow3D', _arrow3D)
co = constants()

class constellation:
    """
    The constellation class can be used to determine the satellite positions
    for each moment of time
    """
    def __init__(self,
                 alt = 550e3,               # satellite altitude [m]
                 Nsat = 66,                 # number of satellites per orbit
                 Norb = 24,                 # number of orbits
                 F = 13,                    # F-parameter to be used initially
                 i = 53 * np.pi / 180,      # inclination [rad]
                 orb = 0,                   # default satellite orbit to be used in the estimations
                 sat = 0,                   # satellite within this orbit to be used in the estimations
                 Nt = 10000):               # number of points in the time axis assumed
        self.alt = alt
        self.r = self.alt + co.Rearth
        self.T = 2 * np.pi * self.r ** (3/2) / np.sqrt(co.mu)
        self.h = 2 * np.pi * self.r ** 2 / self.T
        self.Nsat = Nsat
        self.Norb = Norb
        self.N = self.Nsat * self.Norb
        self.F = F
        self.i = i
        self.sat = sat
        self.orb = orb
       
        k = np.arange( self.Norb )
        self.k = k        
        self.Wk = 2 * np.pi * k / self.Norb
        self.Txx = np.cos( self.Wk )
        self.Tyx = np.sin( self.Wk )
        self.Txy = -np.sin(self.Wk) * np.cos(self.i)
        self.Tyy = np.cos(self.Wk) * np.cos(self.i)
        self.Tzy = -np.sin(i)
        m = np.arange( self.Nsat )    
        self.m = m
        
        self.mm, self.kk = np.meshgrid(m, k)
        self.x = np.zeros(self.mm.shape)
        self.y = np.zeros(self.mm.shape)
        self.z = np.zeros(self.mm.shape)
        
        self.calc_init_anomalies()
        self.init_time(Nt = Nt)
        
    # Initialize time axis
    def init_time(self, Nt = 10000):
        self.Nt = Nt
        self.t = np.linspace(0, self.T, Nt)
        
    # Calculate the initial values of the satellite anomalies    
    def calc_init_anomalies(self):
        self.theta_km0 = 2 * np.pi * (self.mm) / self.Nsat + 2 * np.pi * (self.kk) * self.F / self.N       
    
    # Calculate the anomalies at time t
    def calc_anomalies(self, t):
        self.theta_km = self.theta_km0 + 2 * np.pi / self.T * t
    
    # Calculate satellite position with respect to perifocal plane
    def calc_pos_alt(self):
        coeff = self.h ** 2 / co.mu
        self.ux = coeff * np.cos(self.theta_km)
        self.uy = coeff * np.sin(self.theta_km)
    
    # Calculate satellite position with respect to equitorial plan
    def calc_pos_geo(self):
        for k in range(self.Norb):
            self.x[k, :] = self.Txx[k] * self.ux[k, :] + self.Txy[k] * self.uy[k,:]
            self.y[k, :] = self.Tyx[k] * self.ux[k, :] + self.Tyy[k] * self.uy[k,:]
            self.z[k, :] = self.Tzy * self.uy[k, :]
    
    # Calculate satellite positions at time t
    def calc_pos(self, t):
        self.calc_anomalies(t)
        self.calc_pos_alt()
        self.calc_pos_geo()
    
    # Draw equatorial plane
    def plot_eq_plane(self, ax, N = 100, plot_str='.', show_pos_x = True):
        
        theta = np.linspace(0, 2*np.pi, N)
        x = self.r * np.cos(theta)
        y = self.r * np.sin(theta)
        z = self.r * np.zeros(theta.size)
        
        if show_pos_x:
            i = np.where( (x > 0) )
            x = x[i]
            y = y[i]
            z = z[i]
            
        ax.plot3D(x,y,z,plot_str) 
    
    # Initialize a 3D figure for constellation plotting 
    def init_fig(self, include_axis = False):
        ax = plt.axes(projection='3d')
        s = self.r
        
        if include_axis:
            ax.arrow3D(0,0,0,2 * s,0,0,mutation_scale=20,
                       arrowstyle="-|>",
                       linestyle='dashed')
            
            ax.arrow3D(0,0,0,0,2 * s,0,mutation_scale=20,
                       arrowstyle="-|>",
                       linestyle='dashed')
            
            ax.arrow3D(0,0,0,0,0,2 * s,mutation_scale=20,
                       arrowstyle="-|>",
                       linestyle='dashed')
            
            ax.text(2.2*s,0,0,'$x$')
            ax.text(0,2.2 * s,0,'$y$')
            ax.text(0,0,2.2 * s,'$z$')
        ax.view_init(20, 30, 0)
        ax.axis('off')
        return ax
    
    # Plot a satellite orbit in the 3D coordinate system
    def plot_orb(self, ax, k = 0, plot_str='-o'):
        ax.plot3D(self.x[k, :],self.y[k, :],self.z[k, :],plot_str)
        
    # Plot a single satellite in the 3D coordinate system
    def plot_sat(self, ax, k = 0, m = 0, plot_str='o'):
        ax.plot3D(self.x[k, m],self.y[k, m],self.z[k, m],plot_str) 
        
    # Calculate the distance between all satellites in orbit i and j
    def calc_inter_dist(self, i, j):
        xi = self.x[i, :]
        yi = self.y[i, :]
        zi = self.z[i, :]
        
        xj = self.x[j, :]
        yj = self.y[j, :]
        zj = self.z[j, :]

        xxi, xxj = np.meshgrid(xi, xj)
        yyi, yyj = np.meshgrid(yi, yj)
        zzi, zzj = np.meshgrid(zi, zj)
        
        r = (xxi - xxj) ** 2 + (yyi - yyj) ** 2 + (zzi - zzj) ** 2
        
        return np.sqrt(r)
    
    # Calculate distance between satellite (orb_i, sat_i) and (orb_j, sat_j)
    def calc_sat_dist(self, orb_i, sat_i, orb_j, sat_j):
        xi = self.x[orb_i, sat_i]
        yi = self.y[orb_i, sat_i]
        zi = self.z[orb_i, sat_i]
        xj = self.x[orb_j, sat_j]
        yj = self.y[orb_j, sat_j]
        zj = self.z[orb_j, sat_j]

        r = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
        
        return np.sqrt(r)
    
    # Calculate shortest distance between satellite (orb_i) and (sat_i) and all other satellites     
    def calc_nearest_dist(self, orb_i, sat_i):
        Dx = self.x - self.x[orb_i, sat_i]
        Dy = self.y - self.y[orb_i, sat_i]
        Dz = self.z - self.z[orb_i, sat_i]
        r = np.sqrt(Dx ** 2 + Dy ** 2 + Dz ** 2)
        r[orb_i, sat_i] = np.inf
        rmin = np.min(r)        
        return rmin
    
    # Estimate shortest distance between satellite (orb_i, sat_i) and the satellites of orbit orb_j    
    def calc_sat_orb_dist(self, orb_i, sat_i, orb_j):
        r = self.calc_sat_dist(orb_i, sat_i, orb_j, self.m)
        sat_j = np.argsort(r)
        return r[sat_j], sat_j

    # Optmize the constellation parameter F 
    def optimize_F(self):
        orb = self.orb
        sat = self.sat
        self.rmint = np.zeros([self.Norb, self.Nt])
        
        for F in range(self.Norb):
            self.F = F
            self.calc_init_anomalies()
            print(F)
            for i, t in enumerate(self.t):
                self.calc_pos(t)
                self.rmint[F, i] = self.calc_nearest_dist(orb, sat)
            
            self.rmin = np.min(self.rmint,axis=1)
        self.F = np.argmin(self.rmin)
        self.calc_init_anomalies()
    
    # Estimate the distances of all satellites from satellite (self.orb, self.sat)
    def calc_r(self, F):
        self.F = F
        self.calc_init_anomalies()
        r = np.zeros(self.t.size)
        for i, t in enumerate(self.t):
            self.calc_pos(t)
            r[i] = self.calc_nearest_dist(self.orb, self.sat)
        return r
    
    # Similar to optimize_F but this time uses Python multiprocessing
    def optimize_F_mp(self, no_cores = None):
      
        if not no_cores:
            no_cores = os.cpu_count()
            
        pool = Pool(no_cores)
        F = np.arange(self.Norb)        

        self.rmint=np.array(pool.map(self.calc_r, F))
        self.rmin = np.min(self.rmint,axis=1)        
        self.F = np.argmax(self.rmin)
        self.calc_init_anomalies            
    
    # next satellite in the same orbit
    def next_sat(self, sat):
        return sat+1 if sat < self.Nsat - 1 else 0
 
    # previous satellite in the same orbit
    def prev_sat(self, sat):
        return sat-1 if sat > 0 else self.Nsat
    
    # next orbit
    def next_orb(self, orb):
        return orb+1 if orb < self.Norb - 1 else 0
    
    # previous orbit    
    def prev_orb(self, orb):
        return orb-1 if orb > 0 else self.Norb
    
    # Nearby satellite 3x3 grid formed by the two closest intra-orbit satellites
    # and the closest three upper and three lower orbit satellites
    def sat_grid(self, orb, sat):
        self.calc_pos(0)
        next_orb = self.next_orb(orb)
        prev_orb = self.prev_orb(orb)
        _ , i_upper = self.calc_sat_orb_dist(orb, sat, next_orb)[:3]
        _ , i_lower = self.calc_sat_orb_dist(orb, sat, prev_orb)[:3]
        _ , i_same = self.calc_sat_orb_dist(orb, sat, orb)[:3]
        
        return {
           'U0' : ( next_orb, i_upper[0] ),
           'U1' : ( next_orb, i_upper[1] ),
           'U2' : ( next_orb, i_upper[2] ),
           'L0' : ( prev_orb, i_lower[0] ),
           'L1' : ( prev_orb, i_lower[1] ),
           'L2' : ( prev_orb, i_lower[2] ),
           'S1' : ( orb, i_same[1] ),
           'S2' : ( orb, i_same[2] )
        }
    
    # Calculate satellite distances inside the assumed grid
    def grid_dists(self):
        orb = self.orb
        sat = self.sat
        g = self.sat_grid(orb, sat)
        
        d = {}
        for k in g:
            d[k] = np.zeros(self.Nt)
            
        for i, t in enumerate(self.t):
            self.calc_pos(t)
            for k, (orb_j, sat_j) in g.items():
                d[k][i] = self.calc_sat_dist(orb, sat, orb_j, sat_j)
        self.grid_d = d
        return d        
    
    # Required link distances to maintain connectivity with the grid
    # def req_dists(self, Nt = 10000, 
    #                     orb_i = 0, 
    #                     sat_i = 0, 
    #                     orb_j = 1, 
    #                     neighbours = 3):
    #     self.req_d = np.zeros([Nt, neighbours])
    #     self.req_i = np.zeros([Nt, neighbours])
        
    #     for it, t in enumerate(self.t):
    #         self.calc_pos(t)
    #         r, i = self.calc_sat_orb_dist(orb_i, sat_i, orb_j)
    #         self.req_d[it, :] = r[0:neighbours]
    #         self.req_i[it, :] = i[0:neighbours]
    
    def __repr__(self):
        st = """Constellation with parameters:
- Altitude : {alt} Km
- Satellites per orbit : {Nsat}
- Number of orbits: {Norb}
- Inclination : {i} deg
- F parameter : {F}""".format(alt = self.alt/1e3,
                   Nsat = self.Nsat,
                   Norb = self.Norb,
                   i = self.i * 180 / np.pi,
                   F = self.F)
        return st
    
    def __str__(self):
        return self.__repr__()
        
        
            
        
        
            
        
        
        

        