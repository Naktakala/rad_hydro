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
chiSolverSetBasicOption(phys1,"maximum_dt",2.0e-8)
chiSolverSetBasicOption(phys1,"CFL"       ,0.9)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)