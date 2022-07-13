#!/usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0, '..')
import matplotlib.pyplot as plt

from inputparameters import InputParameters
from eulerian_radhydro import EulerianRadHydro


### INPUTS ###
inp = InputParameters()
inp.mode = 'radhydro'
inp.limiter = 'double minmod'

# Geometry Specifications
inp.geometry = 'slab'   # geometry type
inp.rL = -0.25  # left boundary coordinate [cm]
inp.rR =  0.25   # right boundary coordinate [cm]
inp.N  = 500  # number of cells

# Material Parameters
inp.C_v      = 0.14472799784454     # specific heat [jerks / (cm3 keV)]
inp.gamma    = 5/3                  # thermodynamic constant [cm3 / g]
inp.kappa    = [577.35, 0, 1, 0]    # opacity constants [g / cm2]
inp.kappa_s  = 0                    # scattering opacity [g / cm2]
inp.a        = 0.01372              # radiation constant [jerks / (cm3 kev4)]
inp.c        = 299.792              # speed of light [cm / sh]


# Time Specifications
inp.cfl         = 0.3  # cfl factor
inp.T_start     = 0.0   # start time [sh]
inp.T_final     = 5.0   # final time [sh]
inp.relErFactor = 0.2 # Relative change in Er per time step
inp.maxTimeStep = 5.0e-2

# Initial Conditions
inp.rho = lambda r: 1. * (r < 0) + 3.00185103 * (r >= 0)
inp.u   = lambda r: 0.380431331 * (r < 0) + 0.126732249 * (r >= 0)
inp.T   = lambda r: 0.1 * (r < 0) + 0.366260705 * (r >= 0)
inp.Er  = lambda r: 1.37201720e-6 * (r < 0) + 2.46899872e-4 * (r >= 0)

# Boundary Conditions
inp.bc_L_hydro_type = "transmissive"
inp.bc_R_hydro_type = "transmissive"
inp.bc_L_rad_type   = "reflective"
inp.bc_R_rad_type   = "reflective"


### RUN PROBLEM ###
inp.checkInputs()
rh = EulerianRadHydro(inp, DEBUG=False)
rh.run(cycle_stop=None, print_freq=1)

### PLOTS ###
rh.fields.plotFields(["T", "Er"], ['b-o', 'r-o'], [-0.02, 0.02])
plt.legend(fontsize=12)