#!/usr/bin/env python3

import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
import copy

class EnergyCoupling:
    #########################################################################
    # Description:
    #   Couple the radiation energy and internal energy together by simultaneously
    #   solving the radiation energy density and the total material energy equations.
    #########################################################################
    def __init__(self, rh, DEBUG=False):
        self.rh     = rh
        self.fields = rh.fields
        self.geo    = rh.geo
        self.mat    = rh.mat
        self.DEBUG  = DEBUG
        self.systemMatrix = None
        self.rhs = None
    
        
    ##### FUNCTIONS #####
    
    def solveCoupledSystem(self, dt, predictor):
        # Update diffusion coefficients to average of densities and update opacities
        if predictor:
            self.mat.computeKappa_a(self.fields.T_old)
            rho_avg = 0.5 * (self.fields.rho_n + self.fields.rho_old)
            self.mat.computeD_edge(rho_avg, self.mat.kappa_t)
        else:
            self.mat.computeKappa_a(self.fields.T_n)
            self.mat.computeD_edge(self.fields.rho_n, self.mat.kappa_t)
            
        # Shorthand
        V = self.geo.V
        A = self.geo.A
        dr = self.geo.dr
        a = self.mat.a
        c = self.mat.c
        C_v = self.mat.C_v
        kappa_a = self.mat.kappa_a
        D_edge = self.mat.D_edge
        T_old = self.fields.T_old
        Er_old = self.fields.Er_old
        Er_edge = self.fields.Er_edge
        
        # Get correct field variables
        if predictor:
            rho = self.fields.rho_n
            u = self.fields.u_n
            u_p = self.fields.u_old
            e = self.fields.e_n
            P = self.fields.P_n
            T = self.fields.T_n
            Em = self.fields.Em_n
            Er = self.fields.Er_n
            dt *= 0.5
        else: 
            rho = self.fields.rho
            u = self.fields.u
            u_p = self.fields.u_n
            e = self.fields.e
            P = self.fields.P
            T = self.fields.T
            Em = self.fields.Em
            Er = self.fields.Er
            dt = dt
            
            
        # Diagonals of coupled equations 
        if self.rh.inp.mode == "radhydro":
            maindiag = np.zeros(2*self.geo.N)
            upperdiag = np.zeros(2*self.geo.N)
            lowerdiag = np.zeros(2*self.geo.N)
            diag_e_to_Er = np.zeros(2*self.geo.N)
            diag_Er_to_e = np.zeros(2*self.geo.N)
            self.rhs = np.zeros(2*self.geo.N)
        elif self.rh.inp.mode == "rad":
            maindiag = np.zeros(self.geo.N)
            upperdiag = np.zeros(self.geo.N)
            lowerdiag = np.zeros(self.geo.N)
            self.rhs = np.zeros(self.geo.N)
                
        
        # Loop through equations and cells
        for eqn in range(2):
            for i in range(self.geo.N):
                # Alt. index 
                j = i + self.geo.N
        
                ### RADIATION ENERGY EQUATIONS ###
                if eqn == 0:
                    
                    # Time derivative terms
                    maindiag[i] = 1.0
                    if self.rh.inp.mode == "radhydro":
                        self.rhs[i] = Er[i]  
                    elif self.rh.inp.mode == "rad": 
                        self.rhs[i] = Er_old[i]
                        self.rhs[i] += dt * rho[i] * kappa_a[i] * c * a * T_old[i]**4
                    
                    # Radiation absorption terms
                    maindiag[i] += 0.5 * dt * rho[i] * kappa_a[i] * c
                    self.rhs[i] -= 0.5 * dt * rho[i] * kappa_a[i] * c * Er_old[i]
                        
                    # Radiation current terms
                    # Left cell edge current
                    if i > 0:
                        coef_F_L = -2 * c / (dr[i]/D_edge[i,0] + dr[i-1]/D_edge[i-1,1]) * dt/V[i]
                        
                        # i'th unk component of Er in left edge currrent expression
                        maindiag[i] -= 0.5 * A[i] * coef_F_L
                        
                        # i-1'st unk component of Er in left edge current expression 
                        lowerdiag[i-1] =  0.5 * A[i] * coef_F_L
                        
                        # Add old left edge current
                        self.rhs[i] += 0.5 * A[i] * coef_F_L * (Er_old[i] - Er_old[i-1])
                    
                    # Right cell edge current
                    if i < self.geo.N-1:
                        coef_F_R = -2 * c / (dr[i+1]/D_edge[i+1,0] + dr[i]/D_edge[i,1]) * dt/V[i]
                        
                        # i'th unk component of Er in right edge current expression
                        maindiag[i] -= 0.5 * A[i+1] * coef_F_R
                        
                        # i+1'st unk component of Er in right edge current expression
                        upperdiag[i+1] = 0.5 * A[i+1] * coef_F_R
                        
                        # Add old right edge current
                        self.rhs[i] -= 0.5 * A[i+1] * coef_F_R * (Er_old[i+1] - Er_old[i])
                                            
                    # Boundary terms
                    if i == 0:
                        # Handle left boundary current (0 if reflective)
                        if self.rh.inp.bc_L_rad_type == "source":
                            coef_F_L = -2 * c / (4 + dr[0]/D_edge[0,0]) * dt/V[i]
                            
                            # 1st unk component of Er for left boundary current
                            maindiag[i] -= 0.5 * A[0] * coef_F_L
                            
                            # Add contributions from fixed source and old 1st Er component of boundary current
                            self.rhs[i] += 0.5 * A[0] * coef_F_L * (Er_old[0] - self.fields.Er_bL)
                            self.rhs[i] -= 0.5 * A[0] * coef_F_L * self.fields.Er_bL
                        
                    if i == self.geo.N-1:
                        # Handle right boundary current (0 if reflective)
                        if self.rh.inp.bc_R_rad_type == "source":
                            coef_F_R = -2 * c / (4 + dr[-1]/D_edge[-1,1]) * dt/V[i]
                            print(coef_F_R, Er_old[-1])
                            
                            
                            # N'th unk component of Er for left boundary current
                            maindiag[i] -= 0.5 * A[-1] * coef_F_R 
                            
                            # Add contributions from fixed source and old N'th component of boundary current
                            self.rhs[i] -= 0.5 * A[-1] * coef_F_R * (self.fields.Er_bR - Er_old[-1])
                            self.rhs[i] -= 0.5 * A[-1] * coef_F_R * self.fields.Er_bR
                    
                    # Material energy coupling
                    if self.rh.inp.mode == "radhydro":
                        # Material emmission term    
                        diag_e_to_Er[j] = -0.5 * dt * rho[i]*kappa_a[i]*a*c * 4*T[i]**3/C_v
                        self.rhs[i] += 0.5 * dt * rho[i]*kappa_a[i]*a*c * T_old[i]**4
                        self.rhs[i] += 0.5 * dt * rho[i]*kappa_a[i]*a*c * (T[i]**4 - 4*T[i]**3*e[i]/C_v)
                    
                        # Radiation kinetic energy deposition terms
                        self.rhs[i] += 1./3. * dt/V[i] * (A[i+1]*Er_edge[i+1] - A[i]*Er_edge[i]) * u_p[i] 
                
                
                ### MATERIAL ENERGY EQUATIONS ###
                elif eqn == 1 and self.rh.inp.mode == "radhydro":
                    
                    # Time derivative terms
                    maindiag[j] = 1.0
                    self.rhs[j] = Em[i]/rho[i] - 0.5*u[i]**2
                    
                    # Material emmission terms
                    maindiag[j] += 0.5 * dt * kappa_a[i]*a*c * 4*T[i]**3/C_v
                    self.rhs[j] -= 0.5 * dt * kappa_a[i]*a*c * T_old[i]**4
                    self.rhs[j] -= 0.5 * dt * kappa_a[i]*a*c * (T[i]**4 - 4*T[i]**3*e[i]/C_v)
                    
                    # Radiation absorption terms
                    diag_Er_to_e[i] = -0.5 * dt * kappa_a[i] * c
                    self.rhs[j] += 0.5 * dt * kappa_a[i] * c * Er_old[i]
                    
                    # Radiation kinetic energy deposition terms
                    self.rhs[j] -= 1./3. * dt/V[i] * (A[i+1]*Er_edge[i+1] - A[i]*Er_edge[i]) * u_p[i]/rho[i]
    
                elif eqn == 1 and self.rh.inp.mode == "rad":
                    break
                
                
        ### ASSEMBLE FULL MATRIX AND RHS ###
        if self.rh.inp.mode == "radhydro":
            data = np.array([maindiag, upperdiag, lowerdiag, diag_e_to_Er, diag_Er_to_e])
            diags = np.array([0, 1, -1, self.geo.N, -self.geo.N])
        elif self.rh.inp.mode == "rad":
            data = np.array([maindiag, upperdiag, lowerdiag])
            diags = np.array([0, 1, -1])
            
        self.systemMatrix = spdiags(data, diags, len(self.rhs), len(self.rhs), format="csr")
        x = spsolve(self.systemMatrix, self.rhs)
        

        ### UPDATE FIELD VARIABLES ###
        for field in range(2):
            for i in range(self.geo.N):
                if field == 0:
                    Er[i] = x[i]
                elif field == 1 and self.rh.inp.mode == "radhydro":
                    e[i] = x[i+self.geo.N]
                    Em[i] = rho[i] * (0.5*u[i]**2 + e[i])
                    P[i] = self.mat.pressureEOS(rho[i]*e[i])
                    T[i] = self.mat.temperatureEOS(e[i])
                    
                elif field == 1 and self.rh.inp.mode == "rad":
                    break