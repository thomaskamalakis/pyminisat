import numpy as np
from scipy.special import erfinv, erfc
from constants import constants as co

def invQ(y):
    return np.sqrt(2)*erfinv(1-2*y)

def Q(x):
    return 0.5 * erfc(x / np.sqrt(2))

def to_l(x):
    return 10 ** (x/10)

def to_dB(x):
    return 10 * np.log10(x)
    
def to_dBm(x):
    return 10 * np.log10(x/1e-3)

class receiver:
    
    def __init__(self, 
                 target_BER = 1e-12,
                 l = 1.55e-6,
                 eta = 0.8,
                 redB = 20,
                 Rb = 10e9,
                 RL = 100,
                 FndB = 3,
                 TK = 300,
                 DR = 80e-3,
                 nR = 0.8,
                 sR = 1e-6):
        
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
        
        
    def calc_Psens(self):         
        # Calculate some auxiliary parameters
        self.Be = self.Rb/2
        self.target_g = invQ( self.target_BER ) 
        self.Fn  = to_l(self.FndB)
        self.re = to_l(self.redB)
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
    
    def __init__(self,
                 sT = 1e-6,
                 theta_T = None,
                 PT = None,
                 nT = 0.8):
        
        self.sT = sT
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

    def __init__(self, 
                 tx = transmitter(),
                 rx = receiver(),
                 R = 1000e3,
                 LMdB = 3):
    
        self.tx = tx
        self.rx = rx
        self.l = self.rx.l
        self.R = R
        self.LMdB = LMdB
        

    def calc_cg(self):
        self.LPS = (self.l / (4 * np.pi * self.R) ) ** 2 
        self.LPSdB = to_dB(self.LPS)
        self.L = self.tx.nT * self.rx.nR * self.tx.GT * self.rx.GR * self.LPS * self.tx.LT * self.rx.LR
        self.LdB = to_dB(self.L)
    
    def set_R(self, R):
        self.R = R
        
    def set_PT(self,PT):
        self.PT = PT
    
    def set_PR(self,PR):
        self.PR = PR
                   
    def calc_required_PT(self):
        if not self.rx.Psens:
            self.rx.calc_Psens()
        
        PR = self.rx.Psens
        self.calc_cg()
        self.PT_req = PR / self.L
        self.PT_req_dBm = to_dBm(self.PT_req)
        

        
        
            



