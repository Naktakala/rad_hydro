import numpy as np

def computeShockConditions(M0, rho0, T0, LOUD=False):
    # Define constants
    Cv     = 1.   # specific heat
    gamma  = 5./3  # thermodynamic constant
    
    # Define preshock conditions
    e0    = Cv * T0   # presock specific energy density
    cs0   = np.sqrt(e0 * gamma*(gamma-1)) # preshock sound speed (cm/sh)
    v0    = M0 * cs0  # preshock velocity (cm/sh)
    
    
    # Define postshock conditions
    v1    = v0*(gamma-1)/(gamma+1) + 2*cs0**2/(v0*(gamma+1)) # postshock velocity
    rho1  = rho0 * v0/v1 # postshock density
    e1    = e0 + (v0**2 - v1**2)/(2*gamma) # postshock specific energy density
    T1    = e1/Cv # postshock temperature
    cs1   = np.sqrt(e1*gamma*(gamma-1)) # postshock sound speed
    M1    = v1/cs1  # postshock Mach number
    
    # Print preshock and postshock states
    if LOUD:
        print('\n')
        print('Preshock Velocity = %.3e' %v0, '(cm/sh) = ', M0, '(M)')
        print('Postshock Velocity = %.3e' %v1, '(cm/sh) = ', M1, '(M)')
        print('\n')
        print('Preshock Temperature = %.3e' %T0, '(keV)')
        print('Postshock Temperature = %.3e' %T1, '(keV)')
        print('\n')
        print('Preshock Density = %.3e' %rho0, '(g/cm^3)')
        print('Postshock Density = %.3e' %rho1, '(g/cm^3)')
    
    # Check jump conditions
    rmass = (rho0*v0) / (rho1*v1)
    rmom  = (rho0*v0**2 + rho0*e0*(gamma-1)) / (rho1*v1**2 + rho1*e1*(gamma-1))
    rten  = (0.5*rho0*v0**3 + rho0*e0*v0*gamma) / (0.5*rho1*v1**3 + rho1*e1*v1*gamma)
    
    if LOUD:
        print('\n')
        print('Mass Flux Jump Ratio: ',         rmass)
        print('Momentum Flux Jump Ratio: ',     rmom)
        print('Total Energy Flux Jump Ratio: ', rten)
    
    # Transform to frame where preshock velocity is zero
    ss = -v0    # shock speed
    v0 = v0 + ss    # preshock velocity in dynamic frame
    v1 = v1 + ss    # postshock velocity in dynamic frame
    M0 = (v0-ss)/cs0    # preshock Mach number in dynamic frame
    M1 = (v1-ss)/cs1    # postshock Mach number in dynamic frame
    
    if LOUD:
        print('\n')
        print('Dynamic Preshock Velocity = %.3e' %v0, '(cm/sh) = ', M0, '(M)' )
        print('Dynamic Postshock Velocity = %.3e' %v1, '(cm/sh) = %.3e' %M1, '(M)')
    
    # Compute ratios of delta flux to delta solution
    ratio_mass = (rho0*v0 - rho1*v1) / (rho0-rho1)
    
    flux0      = rho0*v0**2 + rho0*e0*(gamma-1)
    flux1      = rho1*v1**2 + rho1*e1*(gamma-1)
    ratio_mom  = (flux0 - flux1) / (rho0*v0 - rho1*v1)
    
    flux0      = 0.5*rho0*v0**3 + rho0*e0*v0*gamma
    flux1      = 0.5*rho1*v1**3 + rho1*e1*v1*gamma
    u0         = 0.5*rho0*v0**2 + rho0*e0
    u1         = 0.5*rho1*v1**2 + rho1*e1
    ratio_ten  = (flux0 - flux1) / (u0 - u1)
    
    # Print flux to solution ratios in dynamic frame
    if LOUD:
        print('\n')
        print('The following verification parameters should all equal the shock speed:')
        print('Dynamic Mass Flux-to-Solution Ratio = %.3e' %ratio_mass)
        print('Dynamic Momentum Flux-to-Solution Ratio = %.3e' %ratio_mom)
        print('Dynamic Total Energy Flux-to-Solution Ratio = %.3e' %ratio_ten)
        
    return v0-ss, rho1, v1-ss, T1

