import numpy as np 
import matplotlib.pyplot as plt 

imh = 0
ctr = 1
iph = 2

bc_type = 0

class Face:
    def __init__(self, area : float, 
                       centroid : np.array, 
                       normal : np.array, 
                       neighbor = -1) -> None:
        self.area = area 
        self.centroid = centroid
        self.normal = normal
        self.neighbor = neighbor

    def IsOnBoundary(self) -> bool:
        if self.neighbor < 0:
            return True
        
        return False
        

class Cell:
    def __init__(self, z_center : float, dz : float, neighbors) -> None:
        x_c = np.array([0.0, 0.0, z_center])
        dx = np.array([0.0, 0.0, dz])
        self.centroid = x_c
        self.volume = dz
        
        self.faces = [Face(1.0, x_c-dx/2, np.array([0.0,0.0,-1.0]), neighbors[0]),
                      Face(1.0, x_c+dx/2, np.array([0.0,0.0,+1.0]), neighbors[1])]

def ComputeDfdf(D_i : float, D_n : float,
             x_ci : np.array, x_cn : np.array,
             x_f : np.array, n_f : np.array) -> float:
    dx_f = x_cn - x_ci
    dx_m = x_f - x_ci 
    dx_p = x_cn - x_f

    s_f   = dx_f/np.linalg.norm(dx_f)**2
    s_f_m = dx_m/np.linalg.norm(dx_m)**2
    s_f_p = dx_p/np.linalg.norm(dx_p)**2

    d_f = np.dot(n_f, s_f)
    d_m = np.dot(n_f, s_f_m)
    d_p = np.dot(n_f, s_f_p)

    D_f = (1.0/d_f) * ((D_i*d_m * D_n*d_p)/(D_i*d_m + D_n*d_p))

    return D_f, d_f

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

# ===================================== Problem specs
N_c = 20
x_min = -1.0
x_max = 1.0

D = np.zeros(N_c) + 3/2
sigma = np.zeros(N_c) + 0.1
q = np.zeros(N_c) + 1.0
BC = [["vacuum", 0.0],
      ["vacuum", 0.0]]

# ===================================== Init data structures
dx = (x_max - x_min)/N_c
x = np.zeros(N_c)
cells = []
for i in range(0, N_c):
    x[i] = i*dx + 0.5*dx 
    if i == 0: # Left bndry
        cells.append(Cell(x[i], dx, [-1,1]))
    elif i == (N_c-1):
        cells.append(Cell(x[i], dx, [N_c-2,-2]))
    else:
        cells.append(Cell(x[i], dx, [i-1,i+1]))

# ===================================== Assemble system
b = np.zeros(N_c)
A = np.zeros([N_c, N_c])

for i in range(0, N_c):
    cell_i : Cell = cells[i]
    V_i  = cell_i.volume
    x_ci = cell_i.centroid

    A[i,i]  = sigma[i]*V_i
    b[i]    = q[i]*V_i
    
    for face_item in cell_i.faces:
        face : Face = face_item
        x_f = face.centroid
        n_f = face.normal
        A_f = face.area
        
        if not face.IsOnBoundary():
            n = face.neighbor
            cell_n : Cell = cells[n]
            x_cn = cell_n.centroid

            D_f, d_f = ComputeDfdf(D[i],D[n],x_ci,x_cn,x_f,n_f)

            A[i,i] += D_f*A_f*d_f
            A[i,n] -= D_f*A_f*d_f

        else:
            bc = BC[abs(face.neighbor)-1]
            print(bc[bc_type], face.neighbor,abs(face.neighbor)-1)
            if bc[bc_type] == "dirichlet":
                phi_f = bc[1]
                dx_m = x_f - x_ci 
                s_f_m = dx_m/np.linalg.norm(dx_m)**2
                d_m = np.dot(n_f, s_f_m)

                A[i,i] += D[i]*A_f*d_m
                b[i]   += D[i]*A_f*d_m*phi_f
            elif bc[bc_type] == "vacuum":
                print("Yes")
                dx_m = x_f - x_ci 
                s_f_m = dx_m/np.linalg.norm(dx_m)**2
                d_m = np.dot(n_f, s_f_m)

                A[i,i] += D[i]*A_f*d_m/(1+2*D[i]*d_m)
            else:
                print("Unsupported BC-type: ", bc[bc_type])
                exit

print("A symmetric=", check_symmetric(A))

phi = np.linalg.solve(A, b)

plt.figure()
plt.plot(x, phi)
plt.savefig("Zphi.png")