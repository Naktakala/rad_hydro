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
chiSolverSetBasicOption(phys1,"maximum_dt"   ,1.0e-2)
chiSolverSetBasicOption(phys1,"CFL"          ,0.3)
chiSolverSetBasicOption(phys1,"max_timesteps",2000)
chiSolverSetBasicOption(phys1,"max_time"     ,0.2)

--vol_R = chiLogicalVolumeCreate(RPP,-10,10,-10,10,0,L/2)
--vol_L = chiLogicalVolumeCreate(RPP,-10,10,-10,10,L/2,L)
--
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

ff_list = {}

table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-rho") )
table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-u") )
table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-v") )
table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-w") )
table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-p") )
table.insert(ff_list, chiGetFieldFunctionHandleByName("CompInFFlow-e") )

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

ofile = io.open("CompInFFlow1D_Test1_output.txt", "w")
io.output(ofile)
io.write("   i           x    density     pressure    velocity    energy  \n")
for k,v in pairs(arrays[1]) do
    io.write(string.format("%4d  %10.2e  %10.2e  %10.2e  %10.2e  %10.2e\n",
                           k, 0.5*(x[k]+x[k+1]), rho[k], p[k], u[k], e[k]))
end
io.close(ofile)

os.execute("mv CompInFFlow1D_Test1_output.txt RadHydro/HydroSolver/RegressionTests/")
os.execute("python3 RadHydro/HydroSolver/RegressionTests/CompInFFlow1D_Test1.py")
os.execute("mv CompInFFlow1D_Test1_output.png RadHydro/HydroSolver/RegressionTests/")
