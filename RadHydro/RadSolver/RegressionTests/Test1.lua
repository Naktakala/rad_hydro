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
N=10
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

--############################################### Setup Physics
phys1 = chiCreateRadTranDOSolver();

chiSolverSetBasicOption(phys1, "num_angles", 10)
chiSolverSetBasicOption(phys1, "scattering_order", 1)
chiSolverSetBasicOption(phys1, "num_groups", 1)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)
