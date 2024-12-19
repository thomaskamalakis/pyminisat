# Pyminisat: A Python package for simulating optical intersatellite links in a minisat constellation

## What is it all about?
This repo contains the minisat `Python` package which can be used to simulate optical intersatellite links in a minisat (Starlink-like) constellation. The code is writter by Mr. Michalis Gioulis and Prof. Thomas Kamalakis both affiliated with Harokopio University of Athens, Department of Informatics and Telematics. The code is available under the MIT open-source licence.

## Requirements and installation
You will need to install `numpy`, `scipy` and `matplotlib` in order to get this to work. To install this in a virtual environment you will first need to install `python-pip`. On Debian-based distributions simply issue:
```
sudo apt install python3-pip
```
Then install the `virtualenv` package using
```
pip3 install venv
```
or if you are on an externally managed environment, use:
```
sudo apt install python3-venv
```
In a folder of your choosing, clone the repo using
```
git clone https://github.com/thomaskamalakis/pyminisat
```
Now, create a virtual environment in the folder and activate it
```
python3 -m venv venv
source venv/bin/activate
```
You can now install the required packages
```
pip3 install numpy scipy matplotlib
```

## Usage
You can use the `Python` script `net_example.py` to run your first simulation. Just issue:
```
python3 net_example.py
```
If all goes well you should see some plots popping out which describe the required optical power for intersatellite optical links in the standard Starlink architecture.

## The `satsim` class

The `satsim` class in `satnet.py` is the main class that can be used to setup the simulation. When initializing the class the user may change any of the default arguments using keyword arguments.

-`LMdB` is the link margin assumed in dB (defaults to 3dB),

-`alt` is the amplitude of the satellites in meters (defaults to 550,000 meters),

-`Nsat` is the number of satellites per orbit (defaults to 66),

-`Norb` is the number of orbits considered in the constellation (defaults to 24),

-`F` is the phasing factor of the constellation (defaults to 13),

-`i` is the inclination in radians (defaults to 0.925 rad, i.e. 53 degrees),

-`target_BER` is the target value of the bit error rate (BER) (defaults to 1e-12),

-`l` is the optical wavelength of the laser link in meters (defaults to 1.55e-6, i.e. 1550nm),

-`eta` is the internal quantum efficiency of the photodetector (defaults to 0.8 - 80%),

-`redB` is the extinction ratio of the laser source in dB (defaults to 20dB),

-`Rb` is the target data rate in b/s (defaults to 10^10 b/s, i.e. 10Gb/s),

-`RL` is the load resistor at the detector in Ohm (defaults to 100 Ohm),

-`FndB` is the noise figure of the detector in dB (defaults to 3dB),

-`TK` is the detector noise temperature is K (defaults to 300K),

-`DR` is the diameter of the receiver in m (defaults to 80e-3, i.e. 80mm),

-`nR` is the efficiency of the detector (defaults to 0.8, i.e. 80%),

-`nT` is the efficiency of the transmitter (defaults to 0.8, i.e. 80%),

-`sR` is the pointing error standard deviation of the detector in rad (defaults to 1e-6, i.e. 1μrad),

-`sT` is the pointing error standard deviation of the transmitter in rad (defaults to 1e-6, i.e. 1μrad),

-`Nt` is the number of points considered in the time axis within one orbital period (defaults to 10000),

-`filename` is the filename for saving the simulation data and results (defaults to `sat.data`),

-`sat` and `orb` is the index of the satellite and the orbit assumed in the calculations. Since the constellation is symmetric, it does not matter which satellite and orbit we specify. Both default to a value of 1, i.e. the satellite (1,1) is assumed in the calculations.

-`no_cores` is the number of cores used (default to the number of cores reported by `os.cpu_count`)

## Initialization of the `satsim` class

To initialize the `satsim` class with the default parameters just use:

```
from satnet import satsim
s = satsim()
```

To change some of the simulation parameters use the keyword arguments as discussed previously. For example,

```
from satnet import satsim
s = satsim(Nsat = 132, Norb = 48)
```

initializes a simulation object where a 48 satellites of 132 satellites per each is assumed.

## Running the simulations

To run the simulation you simply use the `simulate` method of the `satsim` class. 
```
from satnet import satsim
s = satsim(Nsat = 132, Norb = 48)
s.simulate()
```
This method estimates the optimal phasing parameter `F` for the constellation at hand and then estimates the distances between satellites and the required transmission power at every time instance in the time axis. 






## More documentation, please!

We have submitted this work for publication in a Journal so hopefully after it has been accepted we will provide a more detailed implementation.


