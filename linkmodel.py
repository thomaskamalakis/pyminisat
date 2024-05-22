import numpy as np
from scipy.special import erfinv, erfc
from constants import constants as co
from utils import Q, invQ, to_linear, to_dB, to_dBm
"""
Here we implement the transmitter, receiver and link classes
"""


class receiver:
    """
    Basic receiver class
    """
    def __init__(self, 
                 target_BER = 1e-12,    # target BER
                 l = 1.55e-6,           # opweration wavelength [m]
                 eta = 0.8,             # internal quantum efficiency
                 redB = 20,             # assumed extinction ratio [dB]
                 Rb = 10e9,             # data rate [b/s]
                 RL = 100,              # load resistor [Ohm]
                 FndB = 3,              # noise figure [dB]
                 TK = 300,              # absolute temperature [K]
                 DR = 80e-3,            # receiver aperture diameter [m]
                 nR = 0.8,              # receiver efficiency
                 sR = 1e-6):            # receiver pointing error [rad]
        
        self.target_BER = target_BER
        self.l = l
        self.eta = eta
        self.redB = redB
        self.Rb = Rb
        self.RL = RL
        self.FndB = FndB
        self.TK = TK
        self.DR = DR
        self.sR = sR
        self.nR = nR
        self.GR = (self.DR * np.pi / self.l) ** 2
        self.LR = np.exp(-self.GR * self.sR ** 2)
        self.GRdB = to_dB(self.GR)
        self.LRdB = to_dB(self.LR)
        
    # Estimate receiver sensitivity
    def calc_Psens(self):         
        # Calculate some auxiliary parameters
        self.Be = self.Rb/2
        self.target_g = invQ( self.target_BER ) 
        self.Fn  = to_linear(self.FndB)
        self.re = to_linear(self.redB)
        self.f = co.c /  self.l
        self.R = self.eta * co.qe / (co.hP * self.f)
        
        # Noise variance-related parameters:
        self.sth2 = 4 * co.kB * self.Fn * self.TK * self.Be / self.RL
        self.a0 = 4 * self.R * co.qe * self.Be / (1 + self.re)
        self.a1 = 4 * self.R * co.qe * self.Be * self.re / (1 + self.re)
        self.d = self.target_g / (2 * self.R) * (self.re + 1) / (self.re - 1)
        
        # Estimation of the receiver sensitivty
        self.Psens = self.d ** 2 * (self.a0 + self.a1) + 2 * self.d * np.sqrt( self.d ** 2 * self.a0 * self.a1 + self.sth2)
        self.Psens_dBm = to_dBm(self.Psens)   

    # Estimate receiver sensitivity from incident optical power Pavg [W]
    def calc_BER(self, Pavg):
        self.Pavg = Pavg        
        self.P0 = 2 / (self.re + 1) * Pavg
        self.P1 = self.re * self.P0
        self.sh1_2 = 2 * co.qe * self.R * self.P1 * self.Be
        self.sh0_2 = 2 * co.qe * self.R * self.P0 * self.Be
        self.s1 = np.sqrt(self.sth2 + self.sh1_2)
        self.s0 = np.sqrt(self.sth2 + self.sh0_2)
        self.g = self.R * (self.P1 - self.P0) / (self.s1 + self.s0)
        
    def __repr__(self):
        st = """Receiver with parameters:
- target BER : {target_BER}
- wavelength : {l} mum
- quantum efficiency : {eta}
- assumed extinction ration : {redB} dB
- data rate : {Rb} Gb/s
- load resistor : {RL} Ohm
- noise figure : {F} dB
- temperature : {TK} K
- receiver diameter : {DR} mm
- pointing errors : {sR} murad
- receiver efficiency : {nR}
- receiver gain : {GRdB} dB
- pointing loss : {LRdB} dB""".format(
              target_BER = self.target_BER,
              l = self.l / 1e-6,
              eta = self.eta,
              redB = self.redB,
              Rb = self.Rb / 1e9,
              RL = self.RL,
              F = self.FndB,
              TK = self.TK,
              DR = self.DR/1e-3,
              sR = self.sR/1e-6,
              nR = self.nR,
              GRdB = self.GRdB,
              LRdB = self.LRdB)
        return st
        
        def __str__(self):
            return self.__repr__()

class transmitter:
    """
    Basic transmitter class
    """
    def __init__(self,
                 sT = 1e-6,         # transmitter pointing errors [rad]
                 theta_T = None,    # beam divergence [rad]
                 PT = None,         # transmitter optical power [W]
                 nT = 0.8):         # transmitter efficiency
        
        self.sT = sT
        # If beam divergence is not specified choose optimal divergence
        if not theta_T:
            theta_T = 4 * self.sT
        self.theta_T = theta_T
        self.PT = PT
        self.nT = nT
        self.GT = 16 / self.theta_T ** 2
        self.LT = np.exp(-self.GT * self.sT ** 2)
        self.GTdB = to_dB(self.GT)
        self.LTdB = to_dB(self.LT)
        
    def __repr__(self):
        st = """Transmitter with parameters:
- pointing errors : {sT} murad
- beam divergence : {theta_T} murad
- transmitter efficiency : {nT}
- transmitter gain : {GTdB} dB
- trasnmitter loss : {LTdB} dB""".format(
              sT = self.sT/1e-6,
              nT = self.nT,
              theta_T = self.theta_T/1e-6,
              GTdB = self.GTdB,
              LTdB = self.LTdB)
        return st
        
        def __str__(self):
            return self.__repr__()
class link:
    """
    The link class can be used to estimate the link budget for the intersatellite links
    """
    def __init__(self, 
                 tx = transmitter(),        # transmitter object
                 rx = receiver(),           # receiver object
                 R = 1000e3,                # tx/rx distance [m]
                 LMdB = 3):                 # link margin [dB]
    
        self.tx = tx
        self.rx = rx
        self.l = self.rx.l
        self.R = R
        self.LMdB = LMdB
        
    # Calculate free space loss and channel gain
    def calc_cg(self):
        self.LPS = (self.l / (4 * np.pi * self.R) ) ** 2 
        self.LPSdB = to_dB(self.LPS)
        self.L = self.tx.nT * self.rx.nR * self.tx.GT * self.rx.GR * self.LPS * self.tx.LT * self.rx.LR
        self.LdB = to_dB(self.L)
    
    # set link distance
    def set_R(self, R):
        self.R = R
    # Required transmission power to maintain the link               
    def calc_required_PT(self):
        if not self.rx.Psens:
            self.rx.calc_Psens()
        
        PR = self.rx.Psens
        self.calc_cg()
        self.PT_req = PR / self.L
        self.PT_req_dBm = to_dBm(self.PT_req)
        

        
        
            



