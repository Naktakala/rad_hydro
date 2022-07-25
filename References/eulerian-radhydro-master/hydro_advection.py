#!/usr/bin/env python3


from riemann_solver import RiemannSolver
import numpy as np
import sys

class HydroAdvection:
    #########################################################################
    # Description:
    #   Compute the predicted values of the conserved hydro quantities 
    #   considering only hydro effects by solving an the Euler equations
    #   using forward Euler discretization in time with fluxes evaluated at 
    #   the cell edges obtained from the cell edge values obtained from 
    #   generating slopes.
    #########################################################################
    def __init__(self, rh,  DEBUG=False):
        # General setup
        self.rh = rh
        self.geo = rh.geo
        self.mat = rh.mat
        self.fields = rh.fields
        self.ld = rh.ld
        self.DEBUG = DEBUG
        
    
    ##### FUNCTIONS #####
    
    
    #########################################################################
    # Description:
    #   Solve the advection problem over a half time step with lagged edge fluxes.
    #########################################################################
    def advectHydro(self, dt, predictor):
        # Update edge values
        self.fields.updateConservedHydroVector(predictor)
        self.ld.computeEdgeValues(True, predictor)
    
        # Shorthand
        A = self.geo.A
        V = self.geo.V
        U_old = self.fields.U_old
    
        # Get field variables to update
        if predictor:
            rho = self.fields.rho_n
            M = self.fields.M_n
            Em = self.fields.Em_n
            u = self.fields.u_n
            e = self.fields.e_n
            P = self.fields.P_n
            T = self.fields.T_n
        else:
            rho = self.fields.rho
            M = self.fields.M
            Em = self.fields.Em
            u = self.fields.u
            e = self.fields.e
            P = self.fields.P
            T = self.fields.T
                         
        # Compute hydro fluxes and set time step based on predictor/corrector
        U_ld, P_ld = self.fields.U_ld, self.fields.P_ld
        if predictor:
            F_edge = self.fields.computeFluxes(U_ld, P_ld)
            dt *= 0.5
        else:
            F_hllc = self.rh.rs.computeHLLCFluxes(U_ld, P_ld)
            dt = dt
        
        # Advect hydro and update field variable
        U = np.zeros(U_old.shape)        
        for i in range(self.geo.N):
            if predictor:
                F_L, F_R = F_edge[i,:,0], F_edge[i,:,1]
            else:
                F_L, F_R = F_hllc[i], F_hllc[i+1]
                    
            U[i] = U_old[i] - dt/V[i] * (A[i+1]*F_R - A[i]*F_L)

                            
            # Update field variables
            rho[i] = U[i,0]
            M[i] = U[i,1]
            Em[i] = U[i,2]
            u[i] = M[i]/rho[i]
            e[i] = Em[i]/rho[i] - 0.5 * u[i]**2
            P[i] = self.mat.pressureEOS(rho[i]*e[i])
            T[i] = self.mat.temperatureEOS(e[i])

        # if not predictor:
        #     c=0
        #     for Uvec in U:
        #         print(c,Uvec)
        #         if c==50:
        #             print(dt)
        #             print("U_cm1L ", U_ld[c-1,:,0])
        #             print("U_cm1R ", U_ld[c-1,:,1])
        #
        #             print("U_cL ", U_ld[c,:,0])
        #             print("U_cR ", U_ld[c,:,1])
        #
        #             print("U_cp1L ", U_ld[c+1,:,0])
        #             print("U_cp1R ", U_ld[c+1,:,1])
        #
        #             print("F_L ", -F_hllc[c])
        #             print("F_R ",  F_hllc[c+1])
        #         c+=1
        #     sys.exit()

        # Check for negativities
        # if predictor:
        #     if any([irho < 0 for irho in rho]):
        #         print("Negative mass density encountered in predictor step.\n")
        #         sys.exit(1)
        #     elif any([ie < 0 for ie in e]):
        #         print("Negative internal energy density encountered in predictor step.\n")
        #         sys.exit(1)
        # if not predictor:
        #     if any([irho < 0 for irho in rho]):
        #         print("Negative mass density encountered in corrector step.\n")
        #         sys.exit(1)
        #     elif any([ie < 0 for ie in e]):
        #         print("Negative internal energy density encountered in corrector step.\n")
        #         sys.exit(1)