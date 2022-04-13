#!/usr/bin/env python3
import numpy as np

class Materials:
    #########################################################################
    # Description:
    #   Object containing all material parameters and fundamental constants.
    ######################################################################### 
    def __init__(self, rh):
        self.rh = rh
        self.inp = rh.inp
        
        # Specific energy density (constant)
        self.C_v = self.inp.C_v
        
        # Compressibility coefficient (constant)
        self.gamma = self.inp.gamma
        
        # Radiaition constant
        self.a = self.inp.a
        
        # Speed of light
        self.c = self.inp.c
        
        # Kappa function
        k1, k2, k3, n = self.inp.kappa
        kappa_func = lambda T: k1 / (k2 * T**n + k3)
        
        # Define kappa function on each cell
        self.kappa_funcs = []
        for i in range(self.rh.geo.N):
            self.kappa_funcs.append(kappa_func)
        
        # Absorption opacity at cell centers (variable)
        self.kappa_a = np.zeros(self.rh.geo.N)
        
        # Scattering opacity (constant)
        self.kappa_s = np.ones(self.rh.geo.N) * self.inp.kappa_s
        
        # Total opacity at cell edges (variable)
        self.kappa_t = np.zeros([self.rh.geo.N, 2])
        
        # Edge diffusion coefficients
        self.D_edge = np.zeros([self.rh.geo.N, 2])
        
        
    ##### FUNCTIONS #####
    
    
    #########################################################################
    # Description:
    #   Initialize material properties from fields
    #########################################################################
    def initFromFields(self):
        self.computeKappa_a(self.rh.fields.T_old)
        self.computeKappa_t(self.rh.fields.T_old)
        self.computeD_edge(self.rh.fields.rho_old, self.kappa_t)
    
    
    #########################################################################
    # Description:
    #   Equation of state for pressure as a function of mass density times 
    #   internal energy density.
    #########################################################################
    def pressureEOS(self, I):
        return (self.gamma - 1) * I
    
    
    #########################################################################
    # Description:
    #   Equation of state for temperature as a function of specific energy density  
    #########################################################################
    def temperatureEOS(self, e):
        return e / self.C_v
    
    
    #########################################################################
    # Description:
    #   Compute sound speed for a state given by a mass density and pressure
    #########################################################################
    def computeSoundSpeed(self, rho, P):
        return np.sqrt(self.gamma * P / rho)
    
    
     #########################################################################
    # Description:
    #   Compute the absorption opacity for given temperatures
    #########################################################################
    def computeKappa_a(self, T):
        for i in range(self.rh.geo.N):
            self.kappa_a[i] = self.kappa_funcs[i](T[i])
    
    
    #########################################################################
    # Description:
    #   Compute the total opacity at cell edges by computing edge temperatures,
    #   computing the absorption opacity based upon those, and then taking
    #   the sum of the temperature dependent absorption opacity and the
    #   constant scattering opacity.
    #########################################################################        
    def computeKappa_t(self, T):
        # Left boundary
        if self.inp.bc_L_rad_type == "source":
            T_bL = ((self.rh.fields.Er_bL / self.a + T[0]**4) / 2)**(1/4)
        else:
            T_bL = T[0]
        self.kappa_t[0,0] = self.kappa_funcs[0](T_bL) + self.kappa_s[0]
        
        # Interior 
        for i in range(self.rh.geo.N-1):
            T_edge = ((T[i]**4 + T[i+1]**4) / 2)**(1/4)
            self.kappa_t[i,1] = self.kappa_funcs[i](T_edge) + self.kappa_s[i]
            self.kappa_t[i+1,0] = self.kappa_funcs[i+1](T_edge) + self.kappa_s[i+1]
            
        # Right boundary
        if self.inp.bc_R_rad_type == "source":
            T_bR = ((self.rh.fields.Er_bR / self.a + T[-1]**4) / 2)**(1/4)
        else:
            T_bR = T[-1]
        self.kappa_t[-1,1] = self.kappa_funcs[-1](T_bR) + self.kappa_s[-1]
        
    
    
    #########################################################################
    # Description:
    #   Compute the diffusion coefficient at a cell edges from the left and
    #   right based upon mass densities and opacities derived from temperatures.
    #########################################################################
    def computeD_edge(self, rho, kappa_t):   
        # Left boundary
        self.D_edge[0,0] = 1. / (3 * rho[0] * kappa_t[0,0])
        
        # Interior
        for i in range(self.rh.geo.N-1):
            self.D_edge[i,1] = 1. / (3 * rho[i] * kappa_t[i,1])
            self.D_edge[i+1,0] = 1. / (3 * rho[i+1] * kappa_t[i+1,0])
            
        # Right boundary
        self.D_edge[-1,1] = 1. / (3 * rho[-1] * kappa_t[-1,1])        
    