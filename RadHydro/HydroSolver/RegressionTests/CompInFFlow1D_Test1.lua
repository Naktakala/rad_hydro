-- 1D Compressible Inviscid test
-- Test:
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
            "Expected "..tostring(num_procs)..
            ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end


--############################################### Setup mesh
chiMeshHandlerCreate()

mesh={}
N=100
L=1.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)
chiVolumeMesherSetupOrthogonalBoundaries()

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE, "Gamma")
chiPhysicsMaterialSetProperty(materials[1],"Gamma",SINGLE_VALUE,1.4)

--############################################### Setup Physics
phys1 = chiCreateCompInFFlowSolver();
chiSolverAddRegion(phys1,region1)
chiSolverSetBasicOption(phys1,"maximum_dt"   ,1.0e-3)
chiSolverSetBasicOption(phys1,"CFL"          ,0.2)
chiSolverSetBasicOption(phys1,"max_timesteps",200)

vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)

chiCompInFFlowSetFieldInitialValue(phys1,vol_L,"rho",1.0)
chiCompInFFlowSetFieldInitialValue(phys1,vol_R,"rho",0.125)
--chiCompInFFlowSetFieldInitialValue(phys1,vol_R,"rho",1.0)

chiCompInFFlowSetFieldInitialValue(phys1,vol_L,"u",0.0)
chiCompInFFlowSetFieldInitialValue(phys1,vol_R,"u",0.0)

chiCompInFFlowSetFieldInitialValue(phys1,vol_L,"p",1.0)
chiCompInFFlowSetFieldInitialValue(phys1,vol_R,"p",0.1)
--chiCompInFFlowSetFieldInitialValue(phys1,vol_R,"p",1.0)

--e=p/(gamma-1)/rho  e_L = 1/(0.4*1) = 2.5
--                   e_R = 0.1/(0.4*0.125) = 2

-- E=0.5*rho*u^2 + rho*e  E_L = 2.5
--                        E_R = 0.25

chiSolverInitialize(phys1)
chiSolverExecute(phys1)