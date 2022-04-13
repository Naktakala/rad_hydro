#!/usr/bin/env python3
import os

import numpy as np
import sys

class ReconstructData:
    def __init__(self, rh, DEGUG=False):
        self.rh = rh
        self.geo = rh.geo
        self.mat = rh.mat
        self.fields = rh.fields
        
        
    ##### FUNCTIONS #####
    
     
    #########################################################################
    # Description:
    #   Compute the edge values for hydro or radiation by reconstructing cell
    #   centered quantities with slopes. There are subtle differences in the 
    #   handling of hydro and radiation quantities which are discussed in the
    #   ``constructSlopes'' function here and the ``addGhostCells'' function in
    #   the fields class.
    #########################################################################
    def computeEdgeValues(self, hydro, predictor):
        if hydro:
            if predictor:
                U = self.fields.U_old
            else:
                U = self.fields.U_n
        else:
            if predictor:
                U = self.fields.Er_old
            else:
                U = self.fields.Er_n
                
        # Init storage for edge values. First value corresponds to i,L then i,R
        U_ld = np.zeros([U.shape[0], U.shape[1], 2]) 
        
        # Compute slopes
        slopes = self.constructSlopes(U)

        # Compute edge values
        U_ld[:,:,0] = U - 0.5*slopes
        U_ld[:,:,1] = U + 0.5*slopes
                      
        # If hydro
        if (U_ld.shape[1] == 3):
            # Compute edge pressures from edge conserved hydro variable edge values
            P_ld = np.zeros([U.shape[0], 2])
            rho_e = U_ld[:,2,:] - 0.5 * U_ld[:,1,:]**2 / U_ld[:,0,:]
            P_ld = self.mat.pressureEOS(rho_e)
            
            # Update fields 
            self.fields.U_ld = np.copy(U_ld)
            self.fields.P_ld = np.copy(P_ld)
            self.fields.U_slopes = np.copy(slopes)
        
        # If radiation
        else:
            self.fields.Er_ld = np.copy(U_ld)
            self.fields.Er_slopes = np.copy(slopes)
    
    
    #########################################################################
    # Description:
    #   Generate slopes by taking differences across cell interfaces for hydro 
    #   quantities or radiation. For hydro, 2 ghost cells are used, meaning a slope
    #   is constructed for each interior cell and one ghost cell on either side. 
    #   This is done to obtain a left edge value for boundary cells. For radiation,
    #   One sided slopes are taken in the boundary cells and centered slopes elsewhere.
    #   On the boundaries, a Marshack or reflective condition is enforced.
    #########################################################################
    def constructSlopes(self, U): 
        # Weights for slope weighting including ghost cell
        weights = np.zeros(U.shape[0])
            
        edge_slopes = np.zeros([U.shape[0], U.shape[1], 2])
        slopes = np.zeros([U.shape[0], U.shape[1]])
        
        # Compute left and right slopes, then average for cell i slope        
        for i in range(self.geo.N):
            # Left boundary
            if i == 0:
                if U.shape[1] == 3:
                    if self.rh.inp.bc_L_hydro_type == "fixed":
                        edge_slopes[i,:,0] = U[i] - self.fields.U_bL
                    elif self.rh.inp.bc_R_hydro_type == "transmissive":
                        edge_slopes[i,:,0] = 0.0
                    elif self.rh.inp.bc_R_hydro_type == "periodic":
                        edge_slopes[i,:,0] = U[i] - U[-1]
                elif U.shape[1] == 1:
                    if self.rh.inp.bc_L_rad_type == "reflective":
                        edge_slopes[i,:,0] = 0.0
                    elif self.rh.inp.bc_L_rad_type == "source":
                        edge_slopes[i,:,0] = U[i] - self.fields.Er_bL

            elif 0 < i < self.geo.N-1:
                edge_slopes[i,:,1] = U[i+1] - U[i]
                edge_slopes[i+1,:,0] = U[i+1] - U[i]
            
            # Right boundary
            elif i == self.geo.N:
                if U.shape[1] == 3:
                    if self.rh.inp.bc_R_hydro_type == "fixed":
                        edge_slopes[i,:,1] = self.fields.U_bL - U[i]
                    elif self.rh.inp.bc_R_hydro_type == "transmissive":
                        edge_slopes[i,:,1] = 0.0
                    elif self.rh.inp.bc_R_hydro_type == "periodic":
                        edge_slopes[i,:,0] = U[0] - U[i]
                    edge_slopes[i,:,0] = U[i] - self.fields.U_bL
                elif U.shape[1] == 1:
                    if self.rh.inp.bc_R_rad_type == "reflective":
                        edge_slopes[i,:,1] = 0.0
                    elif self.rh.inp.bc_L_rad_type == "source":
                        edge_slopes[i,:,1] = self.fields.Er_bR - U[i] 
            
            # Compute cell averaged slopes
            slopes[i] += 1./2 * (1 + weights[i]) * edge_slopes[i,:,0]
            slopes[i] += 1./2 * (1 - weights[i]) * edge_slopes[i,:,1]
            
            # Limit slopes
            if self.rh.inp.limiter == "double minmod":
                slopes[i] = self.DoubleMinmod(slopes[i], edge_slopes[i])
            elif self.rh.inp.limiter == "minmod":
                slopes[i] = self.Minmod(slopes[i], edge_slopes[i])
            elif self.rh.inp.limiter == "vanLeer":
                slopes[i] = self.vanLeer(slopes[i], edge_slopes[i], weights[i])
            elif self.rh.inp.limiter == "none":
                pass
        
        return slopes
    
    
    #########################################################################
    # Description:
    #   Double Minmod slope limiter. According to McClarren paper, preserves the 
    #   equilibrium diffusion limit.
    #########################################################################
    def DoubleMinmod(self, slope, edge_slopes):
        for i in range(len(slope)):
            if np.sign(slope[i]) == np.sign(edge_slopes[i,0]) == np.sign(edge_slopes[i,1]):
                a = np.abs(slope[i])
                b = 2 * np.abs(edge_slopes[i,0])
                c = 2 * np.abs(edge_slopes[i,1])
                slope[i] = min(a, b, c)
            else:
                slope[i] = 0.0
        return slope
    
    #########################################################################
    # Description:
    #   Minmod slope limiter.
    #########################################################################
    def Minmod(self, slope, edge_slopes):
        for i in range(len(slope)):
            if np.sign(slope[i]) == np.sign(edge_slopes[i,0]) == np.sign(edge_slopes[i,1]):
                a = np.abs(slope[i])
                b = 1.0 * np.abs(edge_slopes[i,0])
                c = 1.0 * np.abs(edge_slopes[i,1])
                slope[i] = min(a, b, c)
            else:
                slope[i] = 0.0
        return slope
    
    #########################################################################
    # Description:
    #   Van Leer slope limiter
    #########################################################################
    def vanLeer(self, slope, edge_slopes, weight):
        # Compute ratio of left to right slopes
        r = np.zeros(slope.shape)
        for i in range(len(slope)):
            if np.abs(edge_slopes[i,0]) == 0.0 and np.abs(edge_slopes[i,1]) == 0.0:
                r[i] = 1
            elif np.abs(edge_slopes[i,0]) != 0.0 and np.abs(edge_slopes[i,1]) == 0.0:
                r[i] = np.inf
            else:
                r[i] = edge_slopes[i,0]/edge_slopes[i,1]                        
                    
        xi = np.zeros(r.shape)
        for i in range(len(slope)):
            if r[i] <= 0:
                xi[i] = 0
            elif r[i] >= 0.0:
                if np.isinf(r[i]):
                    xi[i] = 0
                else:
                    xi[i] = min( 2*r[i]/(1+r[i]), 2/(1-weight + (1+weight)*r[i]) )
        return xi * slope
    