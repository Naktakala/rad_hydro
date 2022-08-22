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
N1 = 20/4
N2 = 400/4
L=0.5
L1 = (L/2)*(1-0.08*2)
L2 = (L/2)*(0.08*2)
dx1 = L1/N1
dx2 = L2/N2

print("L1dx1",L1,dx1)
print("L2dx2",L2,dx2)
xmin = 0.0
i=1
for k=0,N1 do
    mesh[i] = xmin + k*dx1
    i = i + 1
end
print("Last mesh point: ", i-1, mesh[i-1])
for k=1,2*N2 do
    mesh[i] = L1 + k*dx2
    i = i + 1
end
print("Last mesh point: ", i-1, mesh[i-1])
for k=1,N1 do
    mesh[i] = L1 + 2*L2 + k*dx1
    i = i + 1
end
print("Last mesh point: ", i-1, mesh[i-1])
N = 2*N1 + 2*N2

--mesh={}
--N=1000
--
--L=0.5
--xmin = 0.0
--dx = L/N
--for i=1,(N+1) do
--    k=i-1
--    mesh[i] = xmin + k*dx
--end
chiMeshCreateUnpartitioned1DOrthoMesh(mesh)
chiVolumeMesherExecute();
--os.exit()

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
solver_name = "RadHydroSolverC"
phys1 = chiCreateSolverC(solver_name);

chiSolverSetBasicOption(phys1, "maximum_dt"    , 5e-2)
chiSolverSetBasicOption(phys1, "CFL"           , 0.3)
chiSolverSetBasicOption(phys1, "max_timesteps" , -10)
total_time = 2.0
chiSolverSetBasicOption(phys1, "max_time"      , total_time)
time_vals = ""
Nt = 100/4
dt = total_time/Nt
for i=1,Nt do
    time_vals = time_vals..string.format("%.3g ", dt*i)
end
chiSolverSetBasicOption(phys1, "export_times"  , time_vals)
chiSolverSetBasicOption(phys1, "output_prefix" , "XTest2_")

--vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
--vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)
--
vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)
--0.0137201720
rho0 = 1.0
T0 = 0.1
radE0 = a_const * T0^4

e0 = InternalEGiven_T_Cv(T0,Cv)

cs0 = SoundSpeedGiven_e_gamma(e0, gamma)
u0 = 1.2 * cs0

rho1,T1,u1 = chiRadHydroMakePostShockConditionsRH(Cv, gamma, rho0, T0, u0,
        3.00185103, 3.66260705e-01, 1.26565579e-01)
--rho1,T1,u1 = chiRadHydroMakePostShockConditionsHydroOnly(Cv, gamma, rho0, T0, u0,
--        3.00185103, 3.66260705e-01, 1.26565579e-01)

e1 = Cv*T1

radE1 = a_const * T1^4

--u1 = u1 - 0.00016667 - 0.00004167 --400
--u1 = u1 - 0.00050001 --100

--[0]   rho  3.00187244e+00 u    1.26731956e-01 T    3.66258825e-01 e    5.30079064e-02 radE 2.46894802e-04 E    1.83229493e-01


--rho1=3.00187622e+00
--u1= 1.26730048e-01
--T1=3.66262920e-01
--e1=5.30084991e-02
--radE1=2.46905844e-04

--rho1=3.00185103e+00
--u1= 1.26732249e-01
--T1=3.66260705e-01
--e1=1.83229115e-02
--radE1=2.46899872e-06
--radE1=4.11859125e-06

print(string.format(" rho0  %8.8e",rho0 )..
      string.format(" u0    %8.8e",u0   )..
      string.format(" T0    %8.8e",T0   )..
      string.format(" e0    %8.8e",e0   )..
      string.format(" radE0 %8.8e",radE0))
print(string.format(" rho1  %8.8e",rho1 )..
      string.format(" u1    %8.8e",u1   )..
      string.format(" T1    %8.8e",T1   )..
      string.format(" e1    %8.8e",e1   )..
      string.format(" radE1 %8.8e",radE1))
--os.exit()

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"rho",rho0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"rho",rho1)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"w",u0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"w",u1)

--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"p",1.0)
--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"p",0.1)

chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"temperature",T0)
chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"temperature",T1)

--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_L,"radE",radE0)
--chiRadHydroSolverSetScalarFieldWithLV(phys1,vol_R,"radE",radE1)

zmax = 4
zmin = 5
TRANSMISSIVE = 0
FIXED = 1
bctype = TRANSMISSIVE
chiRadHydroSetBCSetting(phys1,zmin,bctype,rho0,0,0,u0,e0,p0,radE0)
chiRadHydroSetBCSetting(phys1,zmax,bctype,rho1,0,0,u1,e1,p1,radE1)

chiSolverInitialize(phys1)
--os.exit()
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




