#!/usr/bin/env python3
from reconstruct_data import ReconstructData
import matplotlib.pyplot as plt
import numpy as np
import sys

class Fields:
    #########################################################################
    # Description:
    #   Solution handling class. Here, the physical variables and conserved
    #   hydro quantities are initialized and stored. General functions 
    #   pertaining common to several different steps in this algorithm are 
    #   defined here for easy and consistent access.Functions pertaining to solving the  
    #   actual physics equations are located in their respective classes
    #   corresponding to the step of the algorithm.
    #########################################################################
    def __init__(self, rh):
        self.rh = rh
        self.inp = rh.inp
        self.mat = rh.mat
        self.geo = rh.geo
        
        # Init density storage
        self.rho = self.initializeAtCenters(self.inp.rho)
        self.rho_n = np.copy(self.rho)
        self.rho_old = np.copy(self.rho)

        # Init velocity storage
        self.u = self.initializeAtCenters(self.inp.u)
        self.u_n= np.copy(self.u)
        self.u_old = np.copy(self.u)
        
        # Init temperature storage
        self.T = self.initializeAtCenters(self.inp.T)
        self.T_n = np.copy(self.T)
        self.T_old = np.copy(self.T)
        
        # Init specific internal energy storage
        self.e = self.inp.C_v * self.T
        self.e_n = np.copy(self.e)
        self.e_old = np.copy(self.e)
        
        # Init pressure storage
        self.P = self.rh.mat.pressureEOS(self.rho * self.e)
        self.P_n = np.copy(self.P)
        self.P_old = np.copy(self.P)
        
        # Init momentum storage
        self.M = self.rho * self.u
        self.M_n = np.copy(self.M)
        self.M_old = np.copy(self.M)
    
        # Total material energy
        self.Em = self.rho * (0.5 * self.u**2 + self.e)
        self.Em_n = np.copy(self.Em) 
        self.Em_old = np.copy(self.Em)
        
        # Hydro vectors
        self.U_n = np.hstack((self.rho_n, self.M_n, self.Em_n))
        self.U_old = np.copy(self.U_n)
        
        # Hydro boudary conditions
        if self.inp.bc_L_hydro_type == "fixed":
            self.rho_bL = self.inp.rho(self.geo.rL)
            self.u_bL = self.inp.u(self.geo.rL)
            self.T_bL = self.inp.T(self.geo.rL)
            self.e_bL = self.inp.C_v * self.T_bL
            self.P_bL = self.mat.pressureEOS(self.rho_bL * self.e_bL)
            self.M_bL = self.rho_bL * self.u_bL
            self.Em_bL = self.rho_bL * (0.5 * self.u_bL**2 + self.e_bL)
            self.U_bL = np.array([self.rho_bL, self.M_bL, self.Em_bL])
        else:
            self.rho_bL = None
            self.u_bL = None
            self.T_bL = None
            self.e_bL = None
            self.P_bL = None
            self.M_bL = None
            self.Em_bL = None
            self.U_bL = None 
        
        if self.inp.bc_R_hydro_type == "fixed":
            self.rho_bR = self.inp.rho(self.geo.rR)
            self.u_bR = self.inp.u(self.geo.rR)
            self.T_bR = self.inp.T(self.geo.rR)
            self.e_bR = self.inp.C_v * self.T_bR
            self.P_bR = self.mat.pressureEOS(self.rho_bR * self.e_bR)
            self.M_bR = self.rho_bR * self.u_bR
            self.Em_bR = self.rho_bR * (0.5 * self.u_bR**2 + self.e_bR)
            self.cs_bR = self.mat.computeSoundSpeed(self.rho_bR, self.P_bR)
            self.U_bR = np.array([self.rho_bR, self.M_bR, self.Em_bR])
        else:
            self.rho_bR = None
            self.u_bR = None
            self.T_bR = None
            self.e_bR = None
            self.P_bR = None
            self.M_bR = None
            self.Em_bR = None
            self.U_bR = None
        
        
        # Init radiation energy density
        self.Er = self.initializeAtCenters(self.inp.Er)
        self.Er_n = np.copy(self.Er)
        self.Er_old = np.copy(self.Er)
        
        # Set radiation boundary conditions
        if self.inp.bc_L_rad_type == "source":
            if self.inp.bc_L_rad_val == None:
                self.Er_bL = self.inp.Er(self.inp.rL)
            else:
                self.Er_bL = self.inp.bc_L_rad_val
        else:
            self.Er_bL = None

        if self.inp.bc_R_rad_type == "source":
            if self.inp.bc_R_rad_val == None:
                self.Er_bR = self.inp.Er(self.inp.rR)
            else:
                self.Er_bR = self.inp.bc_R_rad_val
        else:
            self.Er_bR = None
        
        # Hydro LD values --- Constructed after init
        self.U_ld, self.P_ld = np.zeros((self.geo.N,3,2)), np.zeros((self.geo.N,2))
        
        # Radiation LD values --- Constructed after init
        self.Er_ld = np.zeros((self.geo.N,1,2))
        
        # Slopes
        self.U_slopes = np.zeros((self.geo.N,3,2))
        self.Er_slopes = np.zeros((self.geo.N,2))
        
        # Radiation edge values -- Constructed after init
        self.Er_edge = np.zeros(self.geo.N+1)
        
        # Energy conservation book keeping
        self.material_energy = []
        self.material_advection = []
        self.radiation_energy = []
        self.radiation_leakage = []
        self.radiation_advection = []
        self.total_energy = []
        
        # Set initial energy in system
        Emat, Erad = 0, 0
        for i in range(self.geo.N):
            Emat += self.geo.V[i] * self.Em_old[i]
            Erad += self.geo.V[i] * self.Er_old[i]
        print("Emat Erad", Emat, Erad)
        self.material_energy.append(Emat)
        self.material_advection.append(0)
        self.radiation_energy.append(Erad)
        self.radiation_leakage.append(0)
        self.radiation_advection.append(0)
        self.total_energy.append(Emat + Erad)
        
        self.total_material_advection = 0
        self.total_radiation_advection = 0                          
        self.total_radiation_leakage = 0
        
        
    ##### FUNCTIONS #####


    #########################################################################
    # Description:
    #   Initialize solution variables at cell centers using specified ICs
    #########################################################################
    def initializeAtCenters(self, function):
        values = np.zeros([self.geo.N, 1])
        if function is not None:
            for i in range(self.geo.N):
                values[i] = function(self.geo.r[i])
        return values
    
    #########################################################################
    # Description:
    #   Move most recent field variables to old field variables.
    #########################################################################
    def stepSolutions(self):
        self.rho_old = np.copy(self.rho)
        self.u_old   = np.copy(self.u)
        self.e_old   = np.copy(self.e)
        self.P_old   = np.copy(self.P)
        self.T_old   = np.copy(self.T)
        self.M_old   = np.copy(self.M)
        self.Em_old  = np.copy(self.Em)
        self.Er_old  = np.copy(self.Er)
        
        
    #########################################################################
    # Description:
    #   Form a vector containing all hydro quantities at cell centers of the mesh.
    #########################################################################
    def updateConservedHydroVector(self, predictor=True):
        if predictor:
            self.U_old = np.hstack((self.rho_old, self.M_old, self.Em_old))
        else:
            self.U_n = np.hstack((self.rho_n, self.M_old, self.Em_n))    
            
    
    #########################################################################
    # Description:
    #   Compute edge radiation values using continuity of current across cell
    #   interfaces. This yields a diffusion coefficient weighted harmonic
    #   mean at each interface. This by definition is only performed for 
    #   interior cell interfaces
    #########################################################################    
    def updateEdgeRadiationEnergy(self, predictor):
        # Shorthand
        dr = self.geo.dr
        D_edge = self.mat.D_edge
        Er_edge = self.Er_edge
        
        # Query correct variables
        if predictor:
            Er = self.Er_old
            rho = self.rho
            T = self.T_old
        else:
            Er = self.Er_n
            rho = self.rho_n
            T = self.T_n
            
        # Left boundary
        if self.inp.bc_L_rad_type == "source":
            kappa_t_L = self.mat.kappa_funcs[0](T[0]) + self.mat.kappa_s[0]
            coef_L = 3 * rho[0] * kappa_t_L * dr[0]
            Er_edge[0] = (3*coef_L*self.Er_bL + 4*Er[0]) / (coef_L + 4)  
        else:
            Er_edge[0] = Er[0]
        
        # Interior edges
        for i in range(1, self.geo.N):
            coef_L = D_edge[i-1,1] / dr[i-1]
            coef_R = D_edge[i,0] / dr[i]
            Er_edge[i] = (coef_R*Er[i] + coef_L*Er[i-1])/(coef_L + coef_R)
           
        # Right boundary
        if self.inp.bc_R_rad_type == "source":
            kappa_t_R = self.mat.kappa_funcs[-1](T[-1]) + self.mat.kappa_s[-1]
            coef_R = 3 * rho[-1] * kappa_t_R * dr[-1]
            Er_edge[-1] = (3*coef_R*self.Er_bR + 4*Er[-1]) / (coef_L + 4)
        else:
            Er_edge[-1] = Er[-1]
            
    
    #########################################################################
    # Description:
    #   Compute edge fluxes from reconstructed edge data for either hydro or
    #   radiation energy. This function generated edge values for the conserved
    #   quantities (hydro/radiation), computes the auxillary variable 
    #   (pressure/velocity), and then computes fluxes from that.
    #########################################################################
    def computeFluxes(self, U_ld, X_ld):
        U_ld = np.atleast_3d(U_ld)
        X_ld = np.atleast_2d(X_ld)
        
        # Init edge fluxes
        F_edge = np.zeros(U_ld.shape)
        
        # Hydro fluxes
        if (U_ld.shape[1] == 3):
            # mass flux
            F_edge[:,0,:] = U_ld[:,1,:]
            # momentum flux
            F_edge[:,1,:] = U_ld[:,1,:]**2/U_ld[:,0,:] + X_ld
            # energy flux
            F_edge[:,2,:] = (U_ld[:,2,:] + X_ld) * U_ld[:,1,:]/U_ld[:,0,:]
            
        # Radiation flux
        else:
            # Radiation energy flux
            F_edge[:,0,:] = U_ld[:,0,:] * X_ld 
            
        return F_edge 
            
    
    #########################################################################
    # Description:
    #   Energy conservation checker
    #########################################################################
    def energyConservationCheck(self):
        # Shorthand
        A = self.geo.A
        V = self.geo.V
        dt = self.rh.dt
        dr = self.geo.dr
        c = self.mat.c
        
        # Shorthand field variables
        Em = self.Em
        Em_np = self.Em_n
        Er = self.Er
        Er_np = self.Er_n
        Er_old = self.Er_old
        Er_n = 0.5 * (Er + Er_old)
        T_np = self.T_n
        rho_np = self.rho_n
        P_np = self.P_n
        u_np = self.u_n
        D_edge = self.mat.D_edge
        
        material_energy = self.material_energy
        material_advection = self.material_advection
        radiation_energy = self.radiation_energy
        radiation_leakage = self.radiation_leakage
        radiation_advection = self.radiation_advection
        total_energy = self.total_energy
        
        # Total material and radiation energy
        Emat, Erad = 0, 0
        for i in range(self.geo.N):
            Emat += V[i] * Em[i] 
            Erad += V[i] * Er[i]
        print("Er0 ErNm1 Em0 EmNm1",Er[0],Er[self.geo.N-1],Em[0],Em[self.geo.N-1])
        material_energy.append(Emat)
        radiation_energy.append(Erad)
        
        # Advection leakage and work energy
        if self.inp.bc_L_hydro_type == "transmissive":
            F_L_madv = (Em_np[0] + P_np[0]) * u_np[0]
            if self.inp.mode == "radhydro":
                F_L_radv = 4./3. * Er_np[0] * u_np[0]
            else:
                F_L_radv = 0.0
        else:
            F_L_madv, F_L_radv = 0.0, 0.0
            
        if self.inp.bc_R_hydro_type == "transmissive":
            F_R_madv = (Em_np[-1] + P_np[-1]) * u_np[-1]
            if self.inp.mode == "radhydro":
                F_R_radv = 4./3. * Er_np[-1] * u_np[-1]
            else:
                F_R_radv = 0.0
        else: 
            F_R_madv, F_R_radv = 0.0, 0.0
        mat_adv = (A[-1] * F_R_madv - A[0] * F_L_madv) * dt
        rad_adv = (A[-1] * F_R_radv - A[0] * F_L_radv) * dt
        material_advection.append(mat_adv)
        radiation_advection.append(rad_adv)
        self.total_material_advection += mat_adv
        self.total_radiation_advection += rad_adv
        
        # Radiation leakage
        if self.inp.bc_L_rad_type == "source":
            coef_F_L = -2 * c / (4 + dr[0] / D_edge[0,0])
            F_L_rad = coef_F_L * (Er_n[0] - self.Er_bL)
        else:
            F_L_rad = 0.0
        if self.inp.bc_R_rad_type == "source":
            coef_F_R = -2 * c / (4 * dr[-1] / D_edge[-1,1])
            F_R_rad = coef_F_R * (self.Er_bR - Er_n[-1])
        else:
            F_R_rad = 0.0
        rad_leakage = (A[-1] * F_R_rad - A[0] * F_L_rad) * dt
        radiation_leakage.append(rad_leakage)
        self.total_radiation_leakage += rad_leakage
        
        # Total energy
        total = Emat + Erad + mat_adv +  rad_adv + rad_leakage
        total_energy.append(total)
        
        # Compute energy balance
        dEmat = Emat - material_energy[0]
        dErad = Erad - radiation_energy[0]
        total_mat_adv = self.total_material_advection
        total_rad_adv = self.total_radiation_advection
        total_rad_leak = self.total_radiation_leakage

        print("dEmat          ",dEmat         ,"\n",
              "dErad          ",dErad         ,"\n",
              "total_mat_adv  ",total_mat_adv ,"\n",
              "total_rad_adv  ",total_rad_adv ,"\n",
              "total_rad_leak ",total_rad_leak,"\n",
              "mat_adv        ",mat_adv       ,"\n",
              "rad_adv        ",rad_adv       ,"\n",
              "rad_leakage    ",rad_leakage   ,"\n",
              )
        
        return dEmat + dErad + total_mat_adv + total_rad_adv + total_rad_leak
        

    #########################################################################
    # Description:
    #   This function plots the desired field variable specified.
    #########################################################################
    def plotFields(self, variables, styles=[], xlims=[], ylims=[],savename="",title=""):
        fig = plt.figure()
        plt.rc("text", usetex=True)
        for i, var in enumerate(variables):
            vals = getattr(self, var)
            if len(variables) == 1 and "Er" in var:
                label = "Radiation Energy"
                plt.ylabel("Radiation Energy (jerks/cm$^3$)", fontsize=16)
            elif "Er" in var:
                vals = (vals / self.inp.a)**(1/4)
                label = "Radiation Temperature"
            elif "T" in var:
                label = "Temperature"
                plt.ylabel("Temperature (keV)", fontsize=16)
            elif "rho" in var:
                label = "Mass Density"
                plt.ylabel("Mass Density (g/cm$^3$)", fontsize=16)
            elif "u" in var:
                label = "Velocity"
                plt.ylabel("Velocity (cm/sh)", fontsize=16)
            elif "e" in var:
                label = "Internal Energy"
                plt.ylabel("Internal Energy (jerks/g)", fontsize=16)
                
            plt.plot(self.geo.r, vals, styles[i], ms=2.5, label=label)
        
        plt.xlabel("r (cm)",fontsize=16)
        if xlims == []:
            plt.xlim([self.geo.rL, self.geo.rR])
        else:
            plt.xlim(xlims)
        if ylims != []:
            plt.ylim(ylims)
        plt.grid()
        plt.legend(fontsize=10)
        plt.title(title)
        if savename != "":
            plt.savefig(savename)
        # plt.show()
        
        