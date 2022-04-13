#!/usr/bin/env python3

from reconstruct_data import ReconstructData
import numpy as np
import sys

class RiemannSolver:
    #########################################################################
    # Description:
    #   Solve the Euler equations over a full time step using fluxes obtained
    #   by solving a Riemann problem at each cell edge with predicted values 
    #   at the half time step.
    #########################################################################
    def __init__(self, rh, DEBUG=False):
        # General setup #
        self.rh = rh
        self.geo = rh.geo
        self.mat = rh.mat
        self.fields = rh.fields
        self.ld = rh.ld
        self.DEBUG = DEBUG

        
    ##### FUNCTIONS #####
    
    #########################################################################
    # Description:
    #   Compute HLLC fluxes based upon the wave speeds for the local Riemann
    #   problem. Compute star region quantities as needed based upon wave 
    #   speed conditions. This function begins by reconstructing the cell
    #   centered data at the predicted half time step value, computing corresponding
    #   fluxes, and then modifying them based on the specifications of the HLLC
    #   approximate Riemann solver.
    #########################################################################      
    def computeHLLCFluxes(self, U_ld, P_ld):
        # Estimate wave speeds
        S = self.estimateWaveSpeeds(U_ld, P_ld)
        
        # Compute edge fluxes 
        F_edge = self.fields.computeFluxes(U_ld, P_ld)
        
        # Init hllc flux storage ordered left to right
        F_hllc = np.zeros([self.geo.N+1, 3])
        
        for i in range(self.geo.N+1):
            # if i==50:
            #     print(S[i,0], S[i,1], S[i,2])
            # If left wave speed > 0, take left flux
            if S[i,0] >= 0:
                if i == 0:
                    if self.rh.inp.bc_L_hydro_type == "fixed":
                        U_L = self.fields.U_bL
                        P_L = self.fields.P_bL
                    elif self.rh.inp.bc_L_hydro_type == "transmissive":
                        U_L = U_ld[i,:,0]
                        P_L = P_ld[i,0]
                    elif self.rh.inp.bc_L_hydro_type == "periodic":
                        U_L = U_ld[-1,:,1]
                        P_L = P_ld[-1,1]
                    F_L = self.fields.computeFluxes(U_L, P_L).flatten()
                else:
                    F_L = F_edge[i-1,:,1]
                    
                F_hllc[i] = F_L

                # if i==50:
                #     print("Left")
                
            # If left wave speed < 0 and star wave speed > 0, take left star flux
            elif (S[i,0] <= 0 and S[i,1] >= 0):
                if i == 0:
                    if self.rh.inp.bc_L_hydro_type == "fixed":
                        U_L = self.fields.U_bL
                        P_L = self.fields.P_bL
                    elif self.rh.inp.bc_L_hydro_type == "transmissive":
                        U_L = U_ld[i,:,0]
                        P_L = P_ld[i,0]
                    elif self.rh.inp.bc_L_hydro_type == "periodic":
                        U_L = U_ld[-1,:,1]
                        P_L = P_ld[-1,1]
                    F_L = self.fields.computeFluxes(U_L, P_L).flatten()
                else:
                    U_L = U_ld[i-1,:,1]
                    P_L = P_ld[i-1,1]
                    F_L = F_edge[i-1,:,1]

                U_L_s = self.computeStarRegionValues(U_L, P_L, S[i,0], S[i,1])
                
                F_hllc[i] = F_L + S[i,0] * (U_L_s - U_L)
                # if i==50:
                #     print("LeftStar")
                #     print("U_L_s",U_L_s)
                
            # If star wave speed < 0 and right wave speed > 0, take right star flux
            elif (S[i,1] <= 0 and S[i,2] >= 0):
                if i == self.geo.N:
                    if self.rh.inp.bc_R_hydro_type == "fixed":
                        U_R = self.fields.U_bR
                        P_R = self.fields.P_bR
                    elif self.rh.inp.bc_R_hydro_type == "transmissive":
                        U_R = U_ld[i,:,1]
                        P_R = P_ld[i,1]
                    elif self.rh.inp.bc_R_hydro_type == "periodic":
                        U_R = U_ld[0,:,0]
                        P_R = P_ld[0,0]
                    F_R = self.fields.computeFluxes(U_R, P_R).flatten()
                else:
                    U_R = U_ld[i,:,0]
                    P_R = P_ld[i,0]   
                    F_R = F_edge[i,:,0]
                    
                U_R_s = self.computeStarRegionValues(U_R, P_R, S[i,2], S[i,1])

                # if i==50:
                #     print("U_R_s",U_R_s)
                #     print("RightStar")
                
                F_hllc[i] = F_R + S[i,2] * (U_R_s - U_R)
                
            # If right wave speed < 0, take right flux
            elif S[i,2] <= 0:
                if i == self.geo.N:
                    if self.rh.inp.bc_R_hydro_type == "fixed":
                        U_R = self.fields.U_bR
                        P_R = self.fields.P_bR
                    elif self.rh.inp.bc_R_hydro_type == "transmissive":
                        U_R = U_ld[i,:,1]
                        P_R = P_ld[i,1]
                    elif self.rh.inp.bc_R_hydro_type == "periodic":
                        U_R = U_ld[0,:,0]
                        P_R = P_ld[0,0]
                    F_R = self.fields.computeFluxes(U_R, P_R).flatten()
                else:
                    F_R = F_edge[i,:,0]
                    
                F_hllc[i] = F_R
                # if i==50:
                #     print("Right")
                                
        return F_hllc
    
    
    #########################################################################
    # Description:
    #   Compute the conserved hydro variable in either left or right star 
    #   region for input to HLLC fluxes of conserved hydro quantities. This 
    #   function is used as necessary to compute values for Riemann problem 
    #   solutions lying in the middle states between waves and is applied 
    #   to single edge quantities.
    #########################################################################
    def computeStarRegionValues(self, U, P, S, S_s):
        # Init star region storage
        U_s = np.zeros(U.shape)
        
        # Define velocity and common coefficient
        u = U[1] / U[0]
        coef = U[0] * (S - u) / (S - S_s)
        
        # Star region density
        U_s[0] = coef * 1
        
        # Star region momentum
        U_s[1] = coef * S_s
        
        # Star region material energy
        U_s[2] = coef * (U[2]/U[0] + (S_s - u) * (S_s + P/(U[0]*(S-u))))
        
        return U_s
    
    
     #########################################################################
    # Description:
    #   Estimate wave speeds for the HLLC approximate Riemann problem. This
    #   function uses approximate edge sound speeds to compute left, right,
    #   and star region wave speeds for the Riemann problem on each cell edge.
    #########################################################################
    def estimateWaveSpeeds(self, U_ld, P_ld):
        # Shorthand
        u_ld = U_ld[:,1,:] / U_ld[:,0,:]
        
        # Compute sound speeds
        c_edge = self.mat.computeSoundSpeed(U_ld[:,0,:], P_ld)
        
        # Init wave speed vector ordered left to right
        S = np.zeros([self.geo.N+1, 3])
    
        # Compute wave speeds
        for i in range(self.geo.N+1):
            # Define left and right values for edge i
            if i == 0:
                if self.rh.inp.bc_L_hydro_type == "fixed":
                    U_L, P_L = self.fields.U_bL, self.fields.P_bL
                    u_L, c_L = self.fields.u_bL, self.mat.computeSoundSpeed(U_L[0], P_L)
                elif self.rh.inp.bc_L_hydro_type == "transmissive":
                    U_L, P_L = U_ld[i,:,0], P_ld[i,0]
                    u_L, c_L = u_ld[i,0], c_edge[i,0]
                elif self.rh.inp.bc_L_hydro_type == "periodic":
                    U_L, P_L = U_ld[-1,:,1], P_ld[-1,1]
                    u_L, c_L = u_ld[-1,1], c_edge[-1,1]
            if i < self.geo.N:
                U_R, P_R = U_ld[i,:,0], P_ld[i,0]
                u_R, c_R = u_ld[i,0], c_edge[i,0]
            if i > 0:
                U_L, P_L = U_ld[i-1,:,1], P_ld[i-1,1]
                u_L, c_L = u_ld[i-1,1], c_edge[i-1,1]
            if i == self.geo.N:
                if self.rh.inp.bc_R_hydro_type == "fixed":
                    U_R, P_R = self.fields.U_bR, self.fields.P_bR
                    u_R, c_R = self.fields.u_bR, self.mat.computeSoundSpeed(U_R[0], P_R)
                elif self.rh.inp.bc_R_hydro_type == "transmissive":
                    U_R, P_R = U_ld[i-1,:,1], P_ld[i-1,1]
                    u_L, c_L = u_ld[i-1,1], c_edge[i-1,1]
                elif self.rh.inp.bc_R_hydro_type == "periodic":
                    U_R, P_R = U_ld[0,0], P_ld[0,0]
                    u_R, c_R = u_ld[0,0], c_edge[0,0]
                
            # Left wave
            S[i,0] = min(u_L - c_L, u_R - c_R)
                    
            # Right wave
            S[i,2] = max(u_L + c_L, u_R + c_R)
                    
            # Star wave
            dP = P_L - P_R
            dF_M = U_L[1]*(S[i,0] - u_L) - \
                    U_R[1]*(S[i,2] - u_R)
            dF_rho = U_L[0]*(S[i,0] - u_L) - \
                        U_R[0]*(S[i,2] - u_R)
            S[i,1] = (dF_M - dP) / dF_rho
        return S