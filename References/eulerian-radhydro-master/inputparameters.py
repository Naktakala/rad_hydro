#!/usr/bin/env python3
import numpy as np
from sys import exit

def isScalar(value):
    return isinstance(value, float) or isinstance(value, int)

class InputParameters:
    #########################################################################
    # Description:
    #   Input parameter handler class. All quantities defined in an input 
    #   deck. Inputs checked with checkInputs function after the problem is 
    #   instanciated. 
    #########################################################################
    def __init__(self):
        # Options
        self.mode       = 'radhydro'
        self.limiter    = 'double minmod'
        
        # Geometry Specifications
        self.geometry   = 'slab'
        self.rL         = None
        self.rR         = None
        self.N          = None

        # Physical Constants
        self.gamma   = 5./3                 # thermodynamic constant for EOS
        self.C_v     = 1.                   # specific heat (jerks / cc-eV)
        self.kappa   = [577.35, 0, 1, 0]    # constants for opacity function of temp
        self.kappa_s = 0                    # scattering opacity (g / cm2)
        self.a       = 0.01372              # radiation constant (jerks / cc-keV4)
        self.c       = 299.792              # speed of light (cm )

        # Boundary Conditions
        self.bc_L_hydro_type    = None  # left boundary condition type
        self.bc_R_hydro_type    = None  # right boundary condition type
        self.bc_L_rad_type      = None  # left boundary condition type
        self.bc_R_rad_type      = None  # right boundary condition type
        self.bc_L_rad_val       = None  # left boundary value
        self.bc_R_rad_val       = None  # right boundary value

        # Initial Conditions
        self.rho    = None
        self.u      = None
        self.e      = None
        self.Er     = None

        # Iteration Parameters
        self.cfl         = None
        self.T_start     = None
        self.T_final     = None
        self.relErFactor  = None
        self.maxTimeStep = None
        self.dt          = None


    def checkInputs(self):
        # Options checks
        if (self.mode not in ['radhydro', 'hydro', 'rad']):
            exit("Need to be in radhydro, hydro, or rad mode.")
        if (self.limiter not in ['double minmod', 'minmod', 'vanLeer', 'none']):
            exit("Limiter not available. Choose 'double minmod', 'minmod', 'vanLeer', or 'none'.")
        # Geometry Checks
        if (self.geometry not in ['slab', 'cylinder', 'sphere']):
            exit("Geometry type {} not supported.".format(self.geometry))

        # Physical Constant Checks
        if (not isScalar(self.gamma) or self.gamma <= 1.0):
            exit("Need to specify gamma as a scalar larger than 1.0")

        # Initial Conditions
        if not callable(self.rho):
            exit("Need to specify rho initial condition as a function.")
        if not callable(self.u):
            exit("Need to specify u initial condition as a function.")
        if not callable(self.T):
            exit("Need to specify t initial condition as a function.")
        if not callable(self.Er):
            exit("Need to specify Er initial condition as a function.")

        # Boundary Conditions
        if self.bc_L_hydro_type not in ['fixed', 'transmissive', 'periodic']:
            exit("Need to specify fixed, transmissive, or periodic for left hydro bdry cond.")
        if self.bc_R_hydro_type not in ['fixed', 'transmissive', 'periodic']:
            exit("Need to specify fixed, transmissive, or periodic for right hydro bdry cond.")
        if self.bc_L_rad_type not in ["reflective", "source"]:
            exit("Need to specify reflective or source for left rad bdry cond.")
        if self.bc_R_rad_type not in ["reflective", "source"]:
            exit("Need to specify reflective or source for right rad bdry cond.")

        # Iteration parameters
        if (not isScalar(self.cfl) or self.cfl <= 0):
            exit("Need to specify a positive, scalar CFL factor.")
        if (not isScalar(self.T_start) or self.T_start < 0):
            exit("Need to specify a positive, scalar start time.")
        if (not isScalar(self.T_final) or self.T_final <= self.T_start):
            exit("Need to specify a scalar end time that is later than the start time.")
        if (not isScalar(self.maxTimeStep) or self.maxTimeStep > self.T_final):
            exit("Need to specify a scalar, positive maximum time step that is smaller than the final time.")
        if (not isScalar(self.dt) and self.dt != None):
            exit("Need to specify a fixed time step, or no time step.")
        
