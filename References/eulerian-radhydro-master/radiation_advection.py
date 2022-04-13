#!/usr/bin/env python3

import numpy as np
import copy

class RadiationAdvection:
    #########################################################################
    # Description:
    #   Compute intermittent predicted values of the radiation energy density
    #   considering only the advection of the material by solving only the 
    #   hyperbolic portion of the grey diffusion equation with forward Euler
    #   for a half time step with fluxes evaluated at cell edges obtained from
    #   generating slopes.
    #########################################################################
    def __init__(self, rh, DEBUG=False):
        self.rh = rh
        self.mat = rh.mat
        self.geo = rh.geo
        self.fields = rh.fields
        self.ld = rh.ld
        self.DEBUG = DEBUG
     
        
    ##### FUNCTIONS #####    
    
    
    #########################################################################
    # Description:
    #   Advect the radiation a half time step using lagged edge radiation advection
    #   edge fluxes. The boundary edge values are derived via continuity of flux, 
    #   and the inner edges via slope generation. In boundary cells, a one-sided
    #   slope to the interior is used to generate edge values.
    #########################################################################
    def advectRadiation(self, dt, predictor):
        # Shorthand
        A = self.geo.A
        V = self.geo.V
        Er_old = self.fields.Er_old
        
        # Update edge values
        self.ld.computeEdgeValues(False, predictor)
        Er_ld = self.fields.Er_ld
        u_ld = self.fields.U_ld[:,1,:] / self.fields.U_ld[:,0,:]
        
        # Construct fluxes and get update field
        if predictor:
            Er = self.fields.Er_n
            F_edge = self.fields.computeFluxes(Er_ld, u_ld)
            dt *= 0.5
        else:
            Er = self.fields.Er
            F_upwind = self.upwindFluxes(Er_ld, u_ld)
            dt = dt
                
        for i in range(self.geo.N):
            if predictor:
                F_L, F_R = F_edge[i,0,0], F_edge[i,0,1]
            else:
                F_L, F_R = F_upwind[i], F_upwind[i+1]
              
            Er[i] = Er_old[i] - 4./3.*dt/V[i] * (A[i+1]*F_R - A[i]*F_L)
        
        
    #########################################################################
    # Description:
    #   Upwind radiation advection fluxes based upon the edge values of velocity.
    #########################################################################
    def upwindFluxes(self, Er_ld, u_ld):
        # Init vector for storage
        F_upwind = np.zeros(self.geo.N+1)
                        
        for i in range(self.geo.N+1):
            # Assign left and and right values
            if i == 0:
                if self.rh.inp.bc_L_hydro_type == "fixed":
                    u_L, u_R = self.fields.u_bL, u_ld[i,0]
                elif self.rh.inp.bc_L_hydro_type == "transmissive":
                    u_L, u_R = u_ld[i,0], u_ld [i,0]
                elif self.rh.inp.bc_L_hydro_type == "periodic":
                    u_L, u_R = u_ld[-1,1], u_ld[i,0]
            elif 0 < i < self.geo.N:
                u_L, u_R = u_ld[i-1,1], u_ld[i,0]
            elif i == self.geo.N:
                if self.rh.inp.bc_R_hydro_type == "fixed":
                    u_L, u_R = u_ld[i-1,1], self.fields.u_bR
                elif self.rh.inp.bc_R_hydro_type == "transmissive":
                    u_L, u_R = u_ld[i-1,1], u_ld[i-1,1] 
                elif self.rh.inp.bc_R_hydro_type == "periodic":
                    u_L, u_R = u_ld[i-1,1], u_ld[0,0]
                    
            # Upwinding
            if u_L > 0 and u_R > 0:
                if i == 0:
                    F_upwind[i] = Er_ld[i,0,0] * u_L
                else:
                    F_upwind[i] = Er_ld[i-1,0,1] * u_L

            # If moving left across boundary
            elif u_L < 0 and u_R < 0:
                if i == self.geo.N:
                    F_upwind[i] = Er_ld[i-1,0,1] * u_R
                else:
                    F_upwind[i] = Er_ld[i,0,0] * u_R
                    
            # If both sides flow into the boundary
            elif (u_ld[i,1] <= 0 <= u_ld[i+1,0]):
                if i == 0:
                    F_upwind[i] = Er_ld[i,0,0] * (u_L + u_R)
                elif i == self.geo.N:
                    F_upwind[i] = Er_ld[i-1,0,1] * (u_L + u_R)
                else:
                    F_upwind[i] = Er_ld[i-1,0,1] * u_L + Er_ld[i,0,0] * u_R
                    
            # If both sides flow away from the boundary
            else:
                F_upwind[i] = 0.0
                
        return F_upwind
    
            
