#!/usr/bin/env python3

import numpy as np
import copy

class MomentumCoupling:
    #########################################################################
    # Description:
    #   Add the contributions of the radiation into the material momentum 
    #   equation.
    #
    # Parameters:
    #   eh =  Eulerian rad-hydro problem object
    #   predictor = Boolean to determine whether in predictor step or not
    #   DEBUG = Boolean for verbose output
    #########################################################################
    def __init__(self, rh, DEBUG=False):
        self.rh = rh
        self.geo = rh.geo
        self.fields = rh.fields
        self.mat = rh.mat
        self.DEBUG = DEBUG
        
    ##### FUNCTIONS #####
    
            
    #########################################################################
    # Description:
    #   Take the previously computed intermittent material momentum and add 
    #   in the radiation momentum deposition to the material momentum. This uses
    #   lagged edge radiation fluxes derived via continuity of flux.
    #########################################################################
    def addRadMomentumDeposition(self, dt, predictor):
        # Update opacities and diffusion coefficient
        if predictor:
            self.mat.computeKappa_t(self.fields.T_old)
            self.mat.computeD_edge(self.fields.rho_old, self.mat.kappa_t)
            self.fields.updateEdgeRadiationEnergy(predictor)
        if not predictor:
            self.mat.computeKappa_t(self.fields.T_n)
            self.mat.computeD_edge(self.fields.rho_n, self.mat.kappa_t)
            self.fields.updateEdgeRadiationEnergy(predictor)
            
        # Shorthand 
        A = self.geo.A
        V = self.geo.V
        
        # Get correct field variables
        if predictor:
            M = self.fields.M_n
            rho = self.fields.rho_n
            u = self.fields.u_n
            dt *= 0.5
        else:
            M = self.fields.M
            rho = self.fields.rho
            u = self.fields.u
            dt = dt
            
        # Edge radiation values
        Er_half = self.fields.Er_edge
        
        # Add radiation momentum deposition terms
        for i in range(self.geo.N):
            M[i] -= (1./3.)*(dt/V[i]) * (A[i+1]*Er_half[i+1] - A[i]*Er_half[i])
            u[i]  = M[i] / rho[i]
