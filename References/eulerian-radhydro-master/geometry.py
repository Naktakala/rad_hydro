#!/usr/bin/env python3
import numpy as np

class Geometry:
    def __init__(self, rh):
        #########################################################################
        # Description:
        #   Spatial mesh handler class. Store all spatial domain and mesh 
        #   information including cell edge and center locations, cell edge areas
        #   and cell volumes.
        #########################################################################
        self.rh = rh
        self.inp = rh.inp
        self.geometry = self.inp.geometry

        # Geometry parameters
        self.N = self.inp.N
        self.rL = self.inp.rL
        self.rR = self.inp.rR

        # Radius at cell edges (user defined)
        self.r_half = np.linspace(self.rL, self.rR, self.N+1)

        # Cell widths
        self.dr = np.zeros(self.N)
        for i in range(self.N):
            self.dr[i] = self.r_half[i+1] - self.r_half[i]

        # Cell center coordinates
        self.r = np.zeros(self.N)
        for i in range(self.N):
            self.r[i] = self.r_half[i] + self.dr[i]/2

        # Areas (defined at edges)
        self.A = np.zeros(self.N+1)
        if (self.geometry == 'slab'):
            self.A.fill(1.0)
        elif (self.geometry == 'cylinder'):
            for i in range(self.N+1):
                self.A[i] = 2 * np.pi * self.r_half[i]
        elif (self.geometry == 'sphere'):
            for i in range(self.N+1):
                self.A[i] = 4 * np.pi * self.r_half[i]**2

        # Volumes (defined over cells)
        self.V = np.zeros(self.N)
        if (self.geometry == 'slab'):
            for i in range(self.N):
                self.V[i] = self.r_half[i+1] - self.r_half[i]
        elif (self.geometry == 'cylinder'):
            for i in range(self.N):
                self.V[i] = np.pi * (self.r_half[i+1]**2 - self.r_half[i]**2)
        elif (self.geometry == 'sphere'):
            for i in range(self.N):
                self.V[i] = 4/3 * np.pi * (self.r_half[i+1]**3 - self.r_half[i]**3)
        