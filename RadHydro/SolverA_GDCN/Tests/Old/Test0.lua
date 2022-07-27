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
N=1000
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

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE, "gamma")
chiPhysicsMaterialSetProperty(materials[1],"gamma",SINGLE_VALUE,1.4)

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE, "Cv")
chiPhysicsMaterialSetProperty(materials[1],"Cv",SINGLE_VALUE,1.0)

function MaterialKappaSFunction(T, mat_id)
   return 0.0
end
function MaterialKappaAFunction(T, mat_id)
    kappa_1 = 577.35
    kappa_2 = 0.0
    kappa_3 = 1.0
    n_exponent = 0.0;

    --return kappa_1/( kappa_2 * T^n_exponent + kappa_3 )
    return 0.0
end

--############################################### Setup Physics
solver_name = "RadHydroSolverA"
phys1 = chiCreateSolverA(solver_name);

chiSolverSetBasicOption(phys1, "maximum_dt"   , 1.0e-2)
chiSolverSetBasicOption(phys1, "CFL"          , 0.3)
chiSolverSetBasicOption(phys1, "max_timesteps", 2000)
chiSolverSetBasicOption(phys1, "max_time"     , 0.2)

--vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
--vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)
--
vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"rho",1.0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"rho",1.0)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"w",1.0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"w",1.0)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"p",1.0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"p",1.0)

zmax = 4
zmin = 5
TRANSMISSIVE = 0
FIXED = 1
bctype = FIXED
chiRadHydroSetBCSetting(phys1,zmin,TRANSMISSIVE)
chiRadHydroSetBCSetting(phys1,zmax,TRANSMISSIVE)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

ff_list = {}

table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-rho") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-u") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-v") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-w") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-p") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-e") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-temperature") )
table.insert(ff_list, chiGetFieldFunctionHandleByName(solver_name.."-radE") )

chiExportMultiFieldFunctionToVTK(ff_list, "ZResults")

--############################################### Line plot
--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT, 0.0,0.0,0.0+dx/2)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0,L-dx/2)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, N)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,ff_list[1]) --rho
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,ff_list[4]) --w
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,ff_list[5]) --p
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,ff_list[6]) --e

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)

arrays = chiFFInterpolationGetValue(cline)

--############################################### Write to file
x = mesh
rho = arrays[1]
u   = arrays[2]
p   = arrays[3]
e   = arrays[4]

ofile = io.open("Test1_output.txt", "w")
io.output(ofile)
io.write("   i           x    density     pressure    velocity    energy  \n")
for k,v in pairs(arrays[1]) do
    io.write(string.format("%4d  %10.2e  %10.2e  %10.2e  %10.2e  %10.2e\n",
            k, 0.5*(x[k]+x[k+1]), rho[k], p[k], u[k], e[k]))
end
io.close(ofile)

os.execute("mv Test1_output.txt RadHydro/SolverA_GDCN/Tests/")
os.execute("python3 RadHydro/SolverA_GDCN/Tests/Test1.py")
os.execute("mv Test1_output.png RadHydro/SolverA_GDCN/Tests/")
