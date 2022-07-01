-- 1D Compressible Inviscid test
-- Test:
num_procs = 1

function script_path()
    local str = debug.getinfo(2, "S").source:sub(2)
    return str:match("(.*/)")
end
dofile(script_path().."Z_Utils.lua")



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
L=0.5
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
gamma = 5/3
Cv = 0.14472799784454
chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE, "gamma")
chiPhysicsMaterialSetProperty(materials[1],"gamma",SINGLE_VALUE,gamma)

chiPhysicsMaterialAddProperty(materials[1],SCALAR_VALUE, "Cv")
chiPhysicsMaterialSetProperty(materials[1],"Cv",SINGLE_VALUE,Cv)

function MaterialKappaSFunction(T, mat_id)
   return 0.0
end
function MaterialKappaAFunction(T, mat_id)
    kappa_1 = 577.35
    kappa_2 = 0.0
    kappa_3 = 1.0
    n_exponent = 0.0;

    return kappa_1/( kappa_2 * T^n_exponent + kappa_3 )
    --return 0.0
end

--############################################### Setup Physics
solver_name = "RadHydroSolverA"
phys1 = chiCreateSolverA(solver_name);

chiSolverSetBasicOption(phys1, "maximum_dt"   , 5e-2)
chiSolverSetBasicOption(phys1, "CFL"          , 0.3)
chiSolverSetBasicOption(phys1, "max_timesteps", 10000)
chiSolverSetBasicOption(phys1, "max_time"     , 5.0)

--vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
--vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)
--
vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)

rho0 = 1.0
T0 = 0.1
radE0 = a_const * T0^4

e0 = InternalEGiven_T_Cv(T0,Cv)

cs0 = SoundSpeedGiven_e_gamma(e0, gamma)
u0 = 1.2 * cs0

u1 = (u0^2 *(gamma - 1) + 2*cs0^2)/(u0*(gamma+1))
rho1 = rho0*u0/u1

e1 = (1/2/gamma)*(u0^2 - u1^2) + e0

cs1 = SoundSpeedGiven_e_gamma(e1, gamma)

T1 = e1/Cv

radE1 = a_const * T1^4

print(string.format(" u0    %8.5g",u0   )..
      string.format(" cs0   %8.5g",cs0  )..
      string.format(" rho0  %8.5g",rho0 )..
      string.format(" e0    %8.5g",e0   )..
      string.format(" T0    %8.5g",T0   )..
      string.format(" radE0 %8.5g",radE0))
print(string.format(" u1    %8.5g",u1   )..
      string.format(" cs1   %8.5g",cs1  )..
      string.format(" rho1  %8.5g",rho1 )..
      string.format(" e1    %8.5g",e1   )..
      string.format(" T1    %8.5g",T1   )..
      string.format(" radE1 %8.5g",radE1))
os.exit()

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"rho",rho0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"rho",rho1)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"w",u0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"w",u1)

--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"p",1.0)
--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"p",0.1)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"temperature",T0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"temperature",T1)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"radE",radE0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"radE",radE1)

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
