#!/usr/bin/env python3
from geometry import Geometry
from fields import Fields
from materials import Materials
from reconstruct_data import ReconstructData
from riemann_solver import RiemannSolver
from hydro_advection import HydroAdvection
from radiation_advection import RadiationAdvection
from momentum_coupling import MomentumCoupling
from energy_coupling import EnergyCoupling



import numpy as np

class EulerianRadHydro:
    #########################################################################
    # Description:
    #   Interface class to run Eulerian frame radiation-hydrodynamics 
    #   problems. This class is inited by simply checking inputs, creating 
    #   the mesh, generating initial conditions, and seting up storage for
    #   information about the progress of the simulation
    #########################################################################
    def __init__(self, inp, DEBUG=False):
        self.inp = inp
        self.inp.checkInputs()
        self.DEBUG = DEBUG
        
        # Init spatial mesh and variables
        self.geo    = Geometry(self)
        self.mat    = Materials(self)
        self.fields = Fields(self)
        self.ld     = ReconstructData(self)
        self.rs     = RiemannSolver(self)
        self.mat.initFromFields()
        
        # Init physics 
        self.hydro_advection = HydroAdvection(self)
        self.radiation_advection = RadiationAdvection(self)
        self.momentum_coupling = MomentumCoupling(self)
        self.energy_coupling = EnergyCoupling(self)

        # Init time parameters
        self.time      = 0
        self.time_list = [inp.T_start]
        self.cycle_num = 0
        self.T_final   = inp.T_final
        self.dt        = 0
        self.dt_list   = []
        self.output_times = []
        
        
        
    #########################################################################
    # Description:
    #   Run the full Eulerian rad-hydro problem
    #########################################################################
    def run(self, print_freq=10, cycle_stop=None):
        print(len(self.output_times))
        # Compute first time step size
        self.computeTimeStep()
        
        while self.time <= self.T_final + 1e-8:
            self.cycle_num += 1
            
            self.runTimeStep()
            energy_diff = self.fields.energyConservationCheck()
            
            if self.cycle_num % print_freq == 0:
                print("==========")
                print("Cycle Number: %i" %self.cycle_num)
                print("Sim. Time: %1.4e" %self.time)
                print("Time Step: %1.4e" %self.dt_list[-1])
                print("Energy conservation check: %1.4e" %energy_diff)
                print("==========\n")

            t = 0
            dt = self.dt_list[-1]
            # print("A",len(self.output_times))
            for output_time in self.output_times:
                t += 1
                if (self.time - dt) < output_time and self.time >= output_time:
                    savename = "Test3a_t" + "{:03d}".format(t) + ".png"
                    self.fields.plotFields(["T", "Er"], 
                                           ['b-o', 'r-o'], 
                                           [-0.25,0.25],
                                           [0.095,0.45],savename,
                                           "time={:.4f} ".format(self.time)+ \
                                            "e_balance={:+.4e} ".format(energy_diff[0]))
            
            if cycle_stop != None:
                if self.cycle_num >= cycle_stop:
                    break
                
            self.computeTimeStep()
            self.fields.stepSolutions()


    #########################################################################
    # Description:                                                          
    #   Run a single time step of the MUSCL-Hancock scheme                   
    #########################################################################    
    def runTimeStep(self): 
        for predictor in [True, False]:                       
            if self.inp.mode == 'hydro':
                # Hydro advection
                self.hydro_advection.advectHydro(self.dt, predictor)
            
            elif self.inp.mode == 'radhydro':
                # Hydro advection
                self.hydro_advection.advectHydro(self.dt, predictor)
                
                # Radiation advection
                self.radiation_advection.advectRadiation(self.dt, predictor)
                
                # Couple momentum
                self.momentum_coupling.addRadMomentumDeposition(self.dt, predictor)
                
                # Couple energy
                self.energy_coupling.solveCoupledSystem(self.dt, predictor)
                
            elif self.inp.mode == 'rad':
                self.energy_coupling.solveCoupledSystem(self.dt, predictor)
                
            
        
    #########################################################################
    # Description:
    #   Compute the maximum stable time step based upon solutions, cell 
    #   sound speeds, and specified CFL limits
    #########################################################################
    def computeTimeStep(self):
        if self.inp.dt is not None:
            self.dt = self.inp.dt
            self.dt_list.append(self.dt)
            self.time += self.dt
            self.time_list.append(self.time)
        
        elif self.inp.dt is None:        
            # Shorthand
            u = self.fields.u_old
            P = self.fields.P_old
            rho = self.fields.rho_old
            dr = self.geo.dr
            cfl = self.inp.cfl
            
            c_s = self.mat.computeSoundSpeed(rho, P)
            
            # Return limiting time step
            self.dt = min(self.inp.maxTimeStep, 2*cfl*np.min(dr / (np.abs(u) + np.abs(c_s))))
            self.dt_list.append(self.dt)
            self.time += self.dt
            self.time_list.append(self.time)
            
            