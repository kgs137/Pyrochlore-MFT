import time
import os
from numpy import *
import numpy as np
from numba import jit
from copy import deepcopy
try:
    import cupy as cp
except:
    import numpy as cp
Nmu = 4
a = 1.0                                     # Taking unit lattice spacing # 0.354*10**-9 #[m] = 0.354 nm
kB = 1.0
# global eta
# global nearest_neighbors

# Bravais Lattice Vectors
l0 = array([0.0, 0.0, 0.0]);     l1 = array([0, 0.5*a, 0.5*a]);
l2 = array([0.5*a, 0.0, 0.5*a]); l3 = array([0.5*a, 0.5*a, 0.0]);
l = array([l0, l1, l2, l3])
# Basis Vectors
d0 = array([0.0, 0.0, 0.0]);      d1 = array([0.0, 0.25*a, 0.25*a]);
d2 = array([0.25*a, 0.0, 0.25*a]);    d3 = array([0.25*a, 0.25*a, 0.0]);
basis_vectors = array([d0, d1, d2, d3])
# Reciprocal Lattice Vectors
k0 = array([0, 0, 0]);             k1 = (2*pi)*cross(l2, l3)/dot(l1,cross(l2,l3));
k2 = (2*pi)*cross(l3, l1)/dot(l2,cross(l3,l1));         k3 = (2*pi)*cross(l1, l2)/dot(l3,cross(l1,l2));
# Local Spin Reference Vectors (Unit vectors pointing out of the tetrahedron along the local <111> direction, ui.uj = -1/3
u0 = (3**-0.5)*array([-1.0, -1.0, -1.0]);    u1 = (3**-0.5)*array([-1.0, 1.0, 1.0]);
u2 = (3**-0.5)*array([1.0, -1.0, 1.0]);   u3 = (3**-0.5)*array([1.0, 1.0, -1.0]);
local_axis = array([u0, u1, u2, u3])
local_axis_gpu = cp.asarray(local_axis, dtype=float32)

def check_eta(Nx, Ny, Nz,
              D = 0.0, J1 = 0.0, J2 = 0.0, J3a = 0.0, J3b = 0.0, R = -20.0,
              # D = 1.3224, J1 = 3.41, J2 = -0.14, J3a = -0.00466, J3b = 0.4,
              # D = 1.41; J1 = 3.72; J2 = 0.0; J3a = 0.0; J3b = 0.0               # + # H_s Hertog and Gingras 2000
              # D = 1.3224; J1 = 3.41; J2 = -0.14; J3a = 0.25; J3b = 0.025       # + # H_g
              # D = 1.3224; J1 = 3.41; J2 = -0.14; J3a = 0.03; J3b = 0.031       # + # H_g+ Bramwell 2018
              # D = 1.3224; J1 = 3.41; J2 = 0.006; J3a = -0.01416; J3b = 0.10216 # + # H_OP Samarakoon et al. 2020
              # D = 1.3224, J1 = 3.41, J2 = 0.0, J3a = -0.00466, J3b = 0.0439,       # + # H_OP Hallen 2022    <-----------
              # D = 1.3224, J1 = 3.41, J2 = 0.0, J3a = 0.0,     J3b = 0.4,          # + # H_J' Hallen 2022
              # D = 1.3224; J1 = 3.41; J2 = 0.0; J3a = -0.00466; J3b = 0.4           # + # H_para
              # D = 0.0; J1 = -5.7; J2 = 0.0; J3a = 0.0; J3b = 0.0                 # + # H_NN Hallen 2022
              # D = 1.3224; J1 = 0.0; J2 = 0.0; J3a = 0.0; J3b = 0.0               # + # Test
              verbose = False, dipole_mode="real", radius_of_summation=10):
    global eta
    global nearest_neighbors
    try:
        os.mkdir("interaction-matrices")
    except:
        0
    
    N = Nx*Ny*Nz*Nmu                            # 4*8*8*8 ~ 2000
    if verbose == True:
        print(str(N)+" spins.")
    ann = sqrt(2*0.25**2) * a                   # This is the spacing between Ti ions # 0.35355339 * a
    kB = 1.0 #[K/K]                             # Taking unit energy # kB = 8.617333*10**-5 #[eV/K]
    e = 2.71828182846
    alpha = 1.0
    r0 = 1.0

    kF = (2*np.pi)/(3.576)
    
    state_rnd  = zeros((Nx,Ny,Nz,Nmu), dtype=int8)
    state_2i2o  = ones((Nx,Ny,Nz,Nmu), dtype=int8)
    realspace  = zeros((Nx,Ny,Nz,Nmu,3), dtype=float)
    spin_orient = zeros((Nx,Ny,Nz,Nmu,3), dtype=float)
    lattice    = zeros((Nx,Ny,Nz,3), dtype=float)
    
    for x in arange(0, Nx):
        for y in arange(0, Ny):
            for z in arange(0, Nz):
                lattice[x, y, z, :]   = x*l[1] + y*l[2] + z*l[3]
                state_rnd[x, y, z, :] = (-1)**random.randint(100), (-1)**random.randint(100), (-1)**random.randint(100), (-1)**random.randint(100)
                for mu in arange(0,Nmu):
                    realspace[x, y, z, mu, :] = lattice[x, y, z, :] + basis_vectors[mu]
                    spin_orient[x, y, z, mu, :] = local_axis[mu] # used for plotting
                    state_2i2o[x,y,z,mu] = (-1)**mu
    L = max(realspace[:,:,:,:,:].reshape(-1))
    
    # Build Replica translation vectors for Ewald summation
    A = array([Nx*l1, Ny*l2, Nz*l3])
    B = array([Nx*k1, Ny*k2, Nz*k3])
    V = dot(A[0],cross(A[1],A[2]))
    replica_int = []
    replica_r = []
    replica_k = []
    replica_index = []
    replica_rnn = []
    
    for i in arange(-radius_of_summation, radius_of_summation+1):
        for j in arange(-radius_of_summation, radius_of_summation+1):
            for k in arange(-0, 0+1):
                if sqrt(i**2 + j**2 + k**2) <= radius_of_summation:
                    replica_int.append(array([i, j, k]))
                    replica_r.append(i*A[0] + j*A[1] + k*A[2])
                    replica_k.append(i*B[0] + j*B[1] + k*B[2])
                    if sqrt(i**2 + j**2 + k**2) <= sqrt(3): # Nearest replicas only
                        replica_index.append((2*(i*A[0] + j*A[1] + k*A[2])).astype(int))
                        replica_rnn.append(i*A[0] + j*A[1] + k*A[2])
    replica_int = array(replica_int); replica_r = array(replica_r); replica_k = array(replica_k); replica_index = array(replica_index); replica_rnn = array(replica_rnn);
    Nr = replica_int.shape[0]; #print(Nr);
    Nrnn = replica_rnn.shape[0]
    
    n0 = where((replica_int[:,0] == 0) & (replica_int[:,1] == 0) & (replica_int[:,2] == 0))[0][0]
    
    def Baux(r):
        return (1/r**3)*(scipy.special.erfc(alpha*r) + (2*alpha*r/sqrt(pi))*(exp(-2*(alpha**2)*(r**2))))
    
    def Caux(r):
        return (1/r**5)*(3*scipy.special.erfc(alpha*r) + (2*alpha*r/sqrt(pi))*(exp(-2*(alpha**2)*(r**2)))*(3+2*(alpha**2)*(r**2)))
    
    try:
        rij = load("interaction-matrices/rij_"+str(Nx)+str(Ny)+str(Nz)+".npy")
        rijdist = load("interaction-matrices/rijdist_"+str(Nx)+str(Ny)+str(Nz)+".npy")
    except:
        rijdist = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu)) # For checking distance distribution # zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu)) # Initializing as zeros may be a problem
        rij = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu,3))   # Identities are zero
    
        for xi in arange(0,Nx):
            print(xi/Nx)
            for yi in arange(0,Ny):
                for zi in arange(0,Nz):
                    for mui in arange(0,Nmu):
    
                        for xj in arange(0,Nx):
                            for yj in arange(0,Ny):
                                for zj in arange(0,Nz):
                                    for muj in arange(0,Nmu):
    
                                        if not (xi == xj and yi == yj and zi == zj and mui == muj):  #If not-self interaction
                                            ri = realspace[xi,yi,zi,mui,:]          #Unpack real space location vector
                                            rj = realspace[xj,yj,zj,muj,:]
    
                                            rij[xi,yi,zi,mui,xj,yj,zj,muj,:] = rj-ri
                                            rijdist[xi,yi,zi,mui,xj,yj,zj,muj] = dot(rj-ri,rj-ri)**0.5
                                            
        save("interaction-matrices/rijdist_"+str(Nx)+str(Ny)+str(Nz)+".npy", rijdist)
        save("interaction-matrices/rij_"+str(Nx)+str(Ny)+str(Nz)+".npy", rij)
    
    try:
        if radius_of_summation == 0:
            eta = (J1*load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy") +
                   J2*load("interaction-matrices/eta_J2_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy") +
                   J3a*load("interaction-matrices/eta_J3a_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy") +
                   J3b*load("interaction-matrices/eta_J3b_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy"))
            nearest_neighbors = where(load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+"_r0.npy") != 0, 1, 0)
        else:
            eta = (J1*load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy") +
                   J2*load("interaction-matrices/eta_J2_"+str(Nx)+str(Ny)+str(Nz)+".npy") +
                   J3a*load("interaction-matrices/eta_J3a_"+str(Nx)+str(Ny)+str(Nz)+".npy") +
                   J3b*load("interaction-matrices/eta_J3b_"+str(Nx)+str(Ny)+str(Nz)+".npy")) 
            nearest_neighbors = where(load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+".npy") != 0, 1, 0)
    except:
        nn = np.zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu), dtype = int8)
        nnn = np.zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu), dtype = int8)
        nnnna = np.zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu), dtype = int8)
        nnnnb = np.zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu), dtype = int8)
        eta_J1 = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))
        eta_J2 = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))
        eta_J3a = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))
        eta_J3b = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))
        for n in arange(0, Nrnn):
            print("n1", n/Nrnn)
            for xi in arange(0, Nx):
                for yi in arange(0, Ny):
                    for zi in arange(0, Nz):
                        for mui in arange(0, Nmu):
                            for xj in arange(0, Nx):
                                for yj in arange(0, Ny):
                                    for zj in arange(0, Nz):
                                        for muj in arange(0, Nmu):
                                            if ann == dot(rij[xi,yi,zi,mui,xj,yj,zj,muj]+replica_rnn[n], rij[xi,yi,zi,mui,xj,yj,zj,muj]+replica_rnn[n])**0.5:
                                                eta_J1[xi,yi,zi,mui, xj,yj,zj,muj] = dot(local_axis[mui], local_axis[muj])
                                                nn[xi,yi,zi,mui, xj,yj,zj,muj] = 1
        time0 = time.time()
        for xi, yi, zi, mui, xj, yj, zj, muj in list(zip(*where(nn==1))):
            for x2, y2, z2, mu2 in list(zip(*where(nn[xi,yi,zi,mui,:,:,:,:]==0))):
                for n2 in arange(0, Nrnn):
                    rplusn = dot(rij[xj,yj,zj,muj,x2,y2,z2,mu2]+replica_rnn[n2],rij[xj,yj,zj,muj,x2,y2,z2,mu2]+replica_rnn[n2])**0.5
                    if rplusn == ann and mui == mu2 and (xi!=x2 or yi!=y2 or zi!=z2):
                        eta_J3a[xi,yi,zi,mui, x2,y2,z2,mu2] = dot(local_axis[mui], local_axis[mu2])
                        nnnna[xi,yi,zi,mui, x2,y2,z2,mu2] = 1
                    elif rplusn == ann and mui != mu2:
                        eta_J2[xi,yi,zi,mui, x2,y2,z2,mu2] = dot(local_axis[mui], local_axis[mu2])
                        nnn[xi,yi,zi,mui, x2,y2,z2,mu2] = 1
    
        print(time.time() - time0)
        time0 = time.time()    
        for xi, yi, zi, mui, x2, y2, z2, mu2 in list(zip(*where(nnn==1))):
            for x3,y3,z3,mu3 in list(zip(*where((nnn[xi,yi,zi,mui,:,:,:,:]==0) & (nnnna[xi,yi,zi,mui,:,:,:,:]==0)))):
                for n3 in arange(0, Nrnn):
                    rplusn = dot(rij[x2,y2,z2,mu2,x3,y3,z3,mu3]+replica_rnn[n3], rij[x2,y2,z2,mu2,x3,y3,z3,mu3]+replica_rnn[n3])**0.5
                    if rplusn == ann and mui == mu3 and (xi!=x3 or yi!=y3 or zi!=z3):
                        eta_J3b[xi,yi,zi,mui, x3,y3,z3,mu3] = dot(local_axis[mui], local_axis[mu3])
                        nnnnb[xi,yi,zi,mui, x3,y3,z3,mu3] = 1
        print(time.time() - time0)
                                                           
        print("First NN: ", Nx*Ny*Nz*Nmu*6, sum(nn.reshape(-1)))
        print("Second NN: ", Nx*Ny*Nz*Nmu*6*2, sum(nnn.reshape(-1)))
        print("Third-a NN: ", Nx*Ny*Nz*Nmu*6, sum(nnnna.reshape(-1)))
        print("Third-b NN: ", Nx*Ny*Nz*Nmu*6, sum(nnnnb.reshape(-1)))
        
        if radius_of_summation == 0:
            save("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy", eta_J1)
            save("interaction-matrices/eta_J2_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy", eta_J2)
            save("interaction-matrices/eta_J3a_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy", eta_J3a)
            save("interaction-matrices/eta_J3b_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy", eta_J3b)
        else:
            save("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy", eta_J1)
            save("interaction-matrices/eta_J2_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy", eta_J2)
            save("interaction-matrices/eta_J3a_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy", eta_J3a)
            save("interaction-matrices/eta_J3b_"+str(Nx)+str(Ny)+str(Nz)+"_open.npy", eta_J3b)
        eta = J1*eta_J1 + J2*eta_J2 + J3a*eta_J3a + J3b*eta_J3b
        nearest_neighbors = where(load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+".npy") != 0, 1, 0)
    
    # Tests
    # 6 neighbors to a spin
    if sum(nearest_neighbors.reshape(-1))/(6*N) == 1:
        # expected neighbors in bulk and over boundary
        if verbose == True and nearest_neighbors[1,1,1,2,1,1,1,1] == 1 and nearest_neighbors[0,0,0,0,Nx-1,0,0,1] == 1:
            print("Nearest Neighbors: Good")
            
    # product of neighboring <111> axes should be -1/3
    if verbose == True and dot(local_axis[0], local_axis[3]) == min(load("interaction-matrices/eta_J1_"+str(Nx)+str(Ny)+str(Nz)+".npy").reshape(-1)):
        print("Tetrahedra Geometry: Good")
        
    try:
        if dipole_mode == "Ewald":
    #         eta += 2.0*D*(ann**3)*load("1interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+".npy")
            eta_D = load("interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy")
            eta_R = load("interaction-matrices/eta_R_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy")

        else:
    #         eta += 2.0*D*(ann**3)*load("1interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_real.npy")
            eta_D = load("interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy")
            eta_R = load("interaction-matrices/eta_R_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy")

    except:
        time0 = time.time()
        eta_D = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))
        eta_R = zeros((Nx,Ny,Nz,Nmu,Nx,Ny,Nz,Nmu))

        
        # complete the n0-th (central) replica loops first
        replica_kn = replica_k[n0]; replica_rn = replica_r[n0]
        for xi in arange(0, Nx):
            for yi in arange(0, Ny):
                for zi in arange(0, Nz):
                    for mui in arange(0, Nmu):
                        ax_mui = local_axis[mui]
                        eta_D[xi,yi,zi,mui,xi,yi,zi,mui] -= ((2*alpha**3)/(3*sqrt(pi)))*dot(local_axis[mui],local_axis[mui])
    
                        for xj in arange(0, Nx):
                            for yj in arange(0, Ny):
                                for zj in arange(0, Nz):
                                    for muj in arange(0, Nmu):
                                        rij_ = rij[xi,yi,zi,mui,xj,yj,zj,muj,:]
                                        rplusn = dot(rij_+replica_rn, rij_+replica_rn)**0.5
                                    # Ewald Real Space
    #                                     if not (xi==xj and yi==yj and zi==zj and mui==muj):
    #                                         eta_D[xi,yi,zi,mui,xj,yj,zj,muj] += ( (1/2)* (dot(ax_mui, local_axis[muj]) * 
    #                                                                                       (1/rplusn**3)*(math.erfc(alpha*rplusn) + (2*alpha*rplusn/sqrt(pi))*(exp(-2*(alpha**2)*(rplusn**2)))) - 
    #                                                                                       dot(ax_mui,(rij_+replica_rn)) * dot(local_axis[muj],(rij_+replica_rn)) * 
    #                                                                                       (1/rplusn**5)*(3*math.erfc(alpha*rplusn) + (2*alpha*rplusn/sqrt(pi))*(exp(-2*(alpha**2)*(rplusn**2)))*(3+2*(alpha**2)*(rplusn**2)))) )
                                    # Just Real Space
                                        if not (xi==xj and yi==yj and zi==zj and mui==muj):
                                            eta_D[xi,yi,zi,mui,xj,yj,zj,muj] += (dot(ax_mui, local_axis[muj])*(rplusn**-3)) - 3*(dot(ax_mui,(rij[xi,yi,zi,mui,xj,yj,zj,muj,:]+replica_rn)) * dot(local_axis[muj],(rij[xi,yi,zi,mui,xj,yj,zj,muj,:]+replica_rn)))*(rplusn**-5)
                                            eta_R[xi,yi,zi,mui,xj,yj,zj,muj] += dot(ax_mui,local_axis[muj])*(np.cos(2*kF*rplusn)*(rplusn**-3) - np.sin(2*kF*rplusn)*(rplusn**-4)*(2*kF)**-1)

    
        
        @jit
        def over_j(n,xi,yi,zi,mui,replica_kn,replica_rn,recip_coef): # found bug! replicas were not being passed.
            eta_D_temp = zeros((Nx,Ny,Nz,Nmu))
            eta_R_temp = zeros((Nx,Ny,Nz,Nmu))
            ax_mui = local_axis[mui]
            for xj in arange(0, Nx):
                for yj in arange(0, Ny):
                    for zj in arange(0, Nz):
                        for muj in arange(0, Nmu):
                            rij_ = rij[xi,yi,zi,mui,xj,yj,zj,muj,:]
                            rplusn = dot(rij_+replica_rn, rij_+replica_rn)**0.5
                        # Ewald Real Space
    # #                         eta_D[xi,yi,zi,mui,xj,yj,zj,muj] += (1/2)* (dot(ax_mui, local_axis[muj]) * Baux(rplusn) - dot(ax_mui,(rij_+replica_rn)) * dot(local_axis[muj],(rij_+replica_rn)) * Caux(rplusn))
    #                         eta_D_temp[xj,yj,zj,muj] += ( (1/2)* (dot(ax_mui, local_axis[muj]) * 
    #                                                                   (1/rplusn**3)*(math.erfc(alpha*rplusn) + (2*alpha*rplusn/sqrt(pi))*(exp(-2*(alpha**2)*(rplusn**2)))) - 
    #                                                                   dot(ax_mui,(rij_+replica_rn)) * dot(local_axis[muj],(rij_+replica_rn)) * 
    #                                                                   (1/rplusn**5)*(3*math.erfc(alpha*rplusn) + (2*alpha*rplusn/sqrt(pi))*(exp(-2*(alpha**2)*(rplusn**2)))*(3+2*(alpha**2)*(rplusn**2)))) )
                        # Ewald Reciprocal Space
    #                         eta_D_temp[xj,yj,zj,muj] += recip_coef * dot(ax_mui,replica_kn)*dot(local_axis[muj],replica_kn) * e**(cos(dot(replica_kn,rij_)/L))
    #                     # Just Real Space
                            eta_D_temp[xj,yj,zj,muj] += (dot(ax_mui, local_axis[muj])*(rplusn**-3)) - 3*(dot(ax_mui,(rij[xi,yi,zi,mui,xj,yj,zj,muj,:]+replica_rn)) * dot(local_axis[muj],(rij[xi,yi,zi,mui,xj,yj,zj,muj,:]+replica_rn)))*(rplusn**-5)
                            eta_R_temp[xj,yj,zj,muj] += dot(ax_mui,local_axis[muj])*(np.cos(2*kF*rplusn)*(rplusn**-3) - np.sin(2*kF*rplusn)*(rplusn**-4)*(2*kF)**-1)
            return eta_D_temp, eta_R_temp

    # @jit
    # def R_over_j(n,xi,yi,zi,mui,replica_rn):
    #     eta_R_temp = zeros((Nx,Ny,Nz,Nmu))
    #     ax_mui = local_axis[mui]
    #     for xj in arange(0, Nx):
    #         for yj in arange(0, Ny):
    #             for zj in arange(0, Nz):
    #                 for muj in arange(0, Nmu):
    #                     rij_ = rij[xi,yi,zi,mui,xj,yj,zj,muj,:]
    #                     rplusn = dot(rij_+replica_rn, rij_+replica_rn)**0.5
    #                     eta_R_temp[xj,yj,zj,muj] += dot(ax_mui,local_axis[muj])*(np.cos(2*kF*rplusn)*(rplusn**-3) - np.sin(2*kF*rplusn)*(rplusn**-4)*(2*kF)**-1)
    #     return eta_R_temp
    
        print(sum(D*ann**3*eta_D[1,1,1,2,:,:,:,:].reshape(-1)),sum(D*ann**3*eta_D[0,0,0,0,:,:,:,:].reshape(-1))) 
        print(sum(R*ann**3*eta_R[1,1,1,2,:,:,:,:].reshape(-1)),sum(R*ann**3*eta_R[0,0,0,0,:,:,:,:].reshape(-1)))    

    
        for n in arange(0, Nr):
            print("D and R: ", str(round(100*n/Nr,3))+"%")
            replica_kn = replica_k[n]
            replica_rn = replica_r[n]
            # Reciprocal Space
            if n != n0:
                recip_coef = 1.0#(1/(2*L**3)) * (4*pi/dot(replica_kn,replica_kn)) * e**(-(pi**2) *dot(replica_kn,replica_kn)/((L*alpha)**2))
                for xi in arange(0, Nx):
                    for yi in arange(0, Ny):
                        for zi in arange(0, Nz):
                            for mui in arange(0, Nmu):
                                D_, R_ = over_j(n,xi,yi,zi,mui,replica_kn,replica_rn,recip_coef)
                                eta_D[xi,yi,zi,mui, :,:,:,:] += D_
                                eta_R[xi,yi,zi,mui, :,:,:,:] += R_
        print("That took ", time.time() - time0)

    if dipole_mode == "Ewald":
        save("interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy", eta_D)
        eta += 2* D*ann**3*eta_D # A factor of two is missing from the Ewald sum
        eta += 2* R*ann**3*eta_R # A factor of two is missing from the Ewald sum
    else:
        save("interaction-matrices/eta_D_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy", eta_D)
        save("interaction-matrices/eta_R_"+str(Nx)+str(Ny)+str(Nz)+"_r"+str(radius_of_summation)+"_open.npy", eta_R)
        eta +=    D*ann**3*eta_D
        eta +=    R*ann**3*eta_R
    # if gpu_enabled == True:
    try:
        global eta_gpu
        eta_gpu = cp.asarray(eta, dtype=float32)
        global eta_gpu_mid # for subtracting away the self-energy in rates and energy
        eta_gpu_mid = (eta_gpu[
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 0,
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 0]
        + eta_gpu[
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 1,
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 1]
        + eta_gpu[
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 2,
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 2]
        + eta_gpu[
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 3,
        int(round(Nx/2 - 0.5)), int(round(Ny/2 - 0.5)), int(round(Nz/2 - 0.5)), 3])/4
    except:
        0
    # Tests
    # Total coupling for spins at bulk and boundary (Should be ~12 K)
    if verbose == True:
        print(sum(D*ann**3*eta_D[1,1,1,2,:,:,:,:].reshape(-1)),sum(D*ann**3*eta_D[0,0,0,0,:,:,:,:].reshape(-1))) 
        print(sum(R*ann**3*eta_R[1,1,1,2,:,:,:,:].reshape(-1)),sum(R*ann**3*eta_R[0,0,0,0,:,:,:,:].reshape(-1)))
    # return eta

@jit
def add_zeeman_gpu(h_field, Nx,Ny,Nz): # h_field gives the magnitude and direction of the external magnetic field
    global eta_gpu
    for yi in arange(0, Ny):
        for zi in arange(0, Nz):
            for mui in arange(0, Nmu):
                eta_gpu[:,yi,zi,mui, :,yi,zi,mui] += cp.dot(local_axis_gpu[mui], h_field)*cp.identity(Nx, dtype=float32)
    return 0

# def add_zeeman (cpu)

@jit
def calc_energy(state, Nx, Ny, Nz):#, eta_B_ = 0):
    energy = 0.0
    for xi in arange(0,Nx):
        for yi in arange(0,Ny):
            for zi in arange(0,Nz):
                for mui in arange(0,Nmu):
                    sigma_i = state[xi,yi,zi,mui]
                    for xj in arange(0,Nx):
                        for yj in arange(0,Ny):
                            for zj in arange(0,Nz):
                                for muj in arange(0,Nmu):
                                    if not (xi == xj and yi == yj and zi == zj and mui == muj):  #If not-self interaction
                                        if (xi+yi*Nx+zi*Nx*Ny+mui*Nx*Ny*Nz < xj+yj*Nx+zj*Nx*Ny+muj*Nx*Ny*Nz):   #Restricted against double-counting
                                            energy += eta[xi,yi,zi,mui, xj,yj,zj,muj]*sigma_i*state[xj,yj,zj,muj]
    return energy



@jit
def calc_energy_dif(xi, yi, zi, mui, state, Nx,Ny,Nz):
    energy = np.sum(np.multiply(eta[xi,yi,zi,mui,:,:,:,:], state[:,:,:,:]))*state[xi,yi,zi,mui]
    energy -= eta[xi,yi,zi,mui,xi,yi,zi,mui]     # Subtract self-interaction
    return -2*energy

# index conversions
def nutoi(nu, Nx,Ny,Nz):
    mu = int((nu) // (Nx*Ny*Nz))
    z = int((nu - mu*Nx*Ny*Nz) // (Nx*Ny))
    y = int((nu - z*Nx*Ny - mu*Nx*Ny*Nz) // (Nx))
    x = int((nu - y*Nx - z*Nx*Ny - mu*Nx*Ny*Nz) // (1))
    return array([x, y, z, mu])

def itonu(i,j,k,mu, Nx,Ny,Nz):
    return i + j*Nx + k*Nx*Ny + mu*Nx*Ny*Nz

@jit
def flatten_state(state, Nx,Ny,Nz):
    flat_state = np.zeros(Nx*Ny*Nz*Nmu)
    for j in range(Ny):
        for k in range(Nz):
            for mu in range(Nmu):
                flat_state[j*Nx+k*Nx*Ny+mu*Nx*Ny*Nz : j*Nx+k*Nx*Ny+mu*Nx*Ny*Nz + Nx] = state[:,j,k,mu]
    return flat_state

@jit
def calc_monopole_population(state, Nx,Ny,Nz):
    counts_by_category = np.zeros(5, dtype=int) # -4, -2, 0, +2, +4 or AI, 3I, 2I2O, 3O, AO
    vals = array([-4, -2, 0, +2, +4])
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                # Explicit Tetrahedra
                net = sum(state[i,j,k,:].reshape(-1))
                for n in range(5):
                    if net == vals[n]:
                        arg = n
                counts_by_category[arg] += 1
                
                # Implicit Tetrahedra
                net = -1*(state[i,j,k,0] + state[i-1,j,k,1] + state[i,j-1,k,2] + state[i,j,k-1,3])
                for n in range(5):
                    if net == vals[n]:
                        arg = n
                counts_by_category[arg] += 1
                
    return counts_by_category

@jit
def calc_labels(Nx,Ny,Nz):
    labels = zeros(Nx*Ny*Nz*Nmu, dtype=int)
    for i in arange(0, Nx):
        for j in arange(0, Ny):
            for k in arange(0, Nz):
                for mu in arange(0,Nmu):
                    labels[i+j*Nx+k*Nx*Ny+mu*Nx*Ny*Nz] = i+j*Nx+k*Nx*Ny+mu*Nx*Ny*Nz
    return labels


@jit
def calc_rates(state_temp, T_temp, Nx,Ny,Nz): # add eta_mid
    energy_dim = -2*np.multiply(state_temp[:,:,:,:], np.tensordot(eta[:,:,:,:, :,:,:,:], state_temp[:,:,:,:], [[4,5,6,7],[0,1,2,3]]))
    energy_dim -= -2*eta[1,1,1,1, 1,1,1,1]*np.ones((Nx,Ny,Nz,Nmu)) # Subtract self-interaction, a constant
    rates_dim = np.exp(-(energy_dim/(2*kB*T_temp)))
    return rates_dim
    
def calc_rates_gpu(state_temp, T_temp, Nx,Ny,Nz): # send the gpu version of eta_B
    global eta_gpu
    energy_dim = -2*cp.multiply(state_temp[:,:,:,:], cp.tensordot(eta_gpu[:,:,:,:, :,:,:,:], state_temp[:,:,:,:], [[4,5,6,7],[0,1,2,3]]))
    energy_dim -= -2*eta_gpu_mid*cp.ones((Nx,Ny,Nz,Nmu)) # Subtract self-interaction, a constant
    rates_dim = cp.exp(-(energy_dim/(2*kB*T_temp)))
    return rates_dim

def calc_energy_dif_gpu(xi, yi, zi, mui, state):
    global eta_gpu
    energy = cp.sum(cp.multiply(eta_gpu[xi,yi,zi,mui,:,:,:,:], state[:,:,:,:]))*state[xi,yi,zi,mui]
    energy -= eta_gpu_mid     # Subtract self-interaction
    return -2*energy  

@jit
def calc_effective_field(state, Nx,Ny,Nz):
    effective_field = np.zeros((Nx, Ny, Nz, Nmu, 3))
    for xi, yi, zi, mui, xj, yj, zj, muj in list(zip(*where(nearest_neighbors!=0))):
        effective_field[xi,yi,zi,mui, :] += state[xj,yj,zj,muj]*local_axis[muj]
    return effective_field
    
@jit
def rates_mask(rates, state, h = None, threshold = 0.0):
    h = np.cross(local_axis, h, 1, 4)
    h = sum(multiply(h, h), 4)**0.5
    return where(h > threshold, rates, 0)
    
@jit
def update_effective_field(effective_field, state, xj,yj,zj,muj):
    for xi, yi, zi, mui in list(zip(*where(nearest_neighbors[:,:,:,:, xj,yj,zj,muj]!=0))):
        effective_field[xi,yi,zi,mui, :] += state[xj,yj,zj,muj]*local_axis[muj]
    return effective_field

def loopstep(loc, state):
    # Initialize Loop
    nu1 = itonu(loc[0], loc[1], loc[2], loc[3])
    nu_previous = deepcopy(nu1)
    loc1 = deepcopy(loc)
    sigma1 = state[loc1[0],loc1[1],loc1[2],loc1[3]]
    loop = [ nu1 ]
    looping = True
    
    # Determine which tetrahedron to begin with; will alternate henceforth to ensure only two spins per tetrahedron are flipped.
    sign = (-1)**randint(0,2) # This is the truth value of xi==xj and yi==yj and zi==zj.
    while looping == True:
        try:
            # Gather valid neighboring spins
            valid_spins_ = []
            for nu in arange(0, N): # Run over lattice index.
                if nns[nu1, nu] != 0: # Check for nearest neighbor.
                    loc_ = nutoi(nu)
                    if ((sign > 0 and loc1[0]==loc_[0] and loc1[1]==loc_[1] and loc1[2]==loc_[2]) or 
                        (sign < 0 and loc1[0]!=loc_[0] and loc1[1]!=loc_[1] and loc1[2]!=loc_[2])): # Check tetrahedra type.
                        if sigma1 != state[loc_[0],loc_[1],loc_[2],loc_[3]]: # Check for opposite spin.
                            if nu != nu_previous: # Check not previous spin. Alternating sign might already do this.
                                valid_spins_.append(nu)
                else:
                    0
            # Select the next spin in the loop.
            nu_previous = deepcopy(nu1)
            nu1 = array(valid_spins_)[ randint(0,len(valid_spins_)) ]
            # Reassign spin value
            loc1 = nutoi(nu1)
            sigma1 = state[loc1[0],loc1[1],loc1[2],loc1[3]]
            # Alternate tetrahedron sign
            sign *= -1
            # Append the loop / Terminate the loop ... What is a more robust way to do this?
            if nu1 in loop:
                looping = False
            else:
                loop.append(nu1)
        except:
            if valid_spins_ == []:
                return state
    state_output = deepcopy(state)
    # Remove non-loop spins and Flip the loop.
    loop_start = False
    for nu in array(loop):
#         print(nu)
        if nu == nu1:
#             print("nu == nu1")
            loop_start = True
        if loop_start == True:
#             print("flip " + str(nu))
            loc = nutoi(nu)
#             print(loc)
            state_output[loc[0],loc[1],loc[2],loc[3]] *= -1
    return state_output

def spin_ice_mmc(state, Nn = 1000, T = zeros(1000), Nx=8,Ny=8,Nz=8, loops = True, verbosity="full"): # Add T array functionality (Epoch, Series)
    Energy = 0.0;
    states = zeros((Nn+1,Nx,Ny,Nz,Nmu),dtype=int); states[0, :,:,:,:] = deepcopy(state);
    energies = zeros((Nn+1));                      energies[0] = deepcopy(Energy);
    energy_sampling = []
    n = 0
    time0 = time.time() #Initial Time for Progress Bar

    # Initialize last_spin outside of bounds
    last_spin = itonu(Nx, Ny, Nz, Nmu)

    while n < Nn:
        # Build a trial state by flipping a random spin
        state_trial = deepcopy(states[n, :,:,:,:])
        x1, y1, z1, mu1 = randint(0,Nx), randint(0,Ny), randint(0,Nz), randint(0, Nmu)
        
        # If spin was just flipped, flip a loop
        if last_spin == itonu(x1, y1, z1, mu1) and loops == True:
            state_trial[x1, y1, z1, mu1] *= -1 #Flip the Spin Back
            state_loop = loopstep(array([x1, y1, z1, mu1]), deepcopy(state_trial[:,:,:,:]))
            if not is_same(state_loop, states[n, :,:,:,:]):
                state_trial = deepcopy(state_loop)
#                 print("Loop Flip")
        # Else flip the spin
        else:
            state_trial[x1, y1, z1, mu1] *= -1 #Flip the Trial Spin
            last_spin = itonu(x1, y1, z1, mu1)
        
        E_trial = energies[n] + calc_energy_dif(x1, y1, z1, mu1, states[n,:,:,:,:])

        if T[n] == 0.0 and E_trial <= Energy:
            Energy = deepcopy(E_trial)
            energies[n+1] = deepcopy(E_trial)
            states[n+1, :,:,:,:] = deepcopy(state_trial)
            n += 1

        elif (E_trial-Energy <= 0) or (np.exp(-(E_trial-Energy)/(kB*T[n])) > rand()):
#             print(str(Energy) +" to "+ str(E_trial))
#             print(str(2.71828**-((E_trial-Energy)/(kB*T[n]))) +" vs "+ str(rand()))
            Energy = deepcopy(E_trial);
            energies[n+1] = deepcopy(E_trial)
            states[n+1, :,:,:,:] = deepcopy(state_trial)
            if verbosity == "full":
                progress(n, Nn, time0)   #Progress bar
            n += 1
        energy_sampling.append(deepcopy(Energy))
    if verbosity == "full" or verbosity == "partial":
        print("System Energy: " + str(Energy))
    return states, array(energy_sampling)

def spin_ice_kmc(initial_state, N_steps = 1000, Temp = np.ones(1000), Nx=8,Ny=8,Nz=8, flip_mode="two-rate", use_gpu=True, verbose = False,
                 # D = 1.3224, J1 = 3.41, J2 = -0.14, J3a = -0.00466, J3b = 0.4
              # D = 1.41; J1 = 3.72; J2 = 0.0; J3a = 0.0; J3b = 0.0               # + # H_s Hertog and Gingras 2000
              # D = 1.3224; J1 = 3.41; J2 = -0.14; J3a = 0.25; J3b = 0.025       # + # H_g
              # D = 1.3224; J1 = 3.41; J2 = -0.14; J3a = 0.03; J3b = 0.031       # + # H_g+ Bramwell 2018
              # D = 1.3224; J1 = 3.41; J2 = 0.006; J3a = -0.01416; J3b = 0.10216 # + # H_OP Samarakoon et al. 2020
              # D = 1.3224, J1 = 3.41, J2 = 0.0, J3a = -0.00466, J3b = 0.0439       # + # H_OP Hallen 2022    <-----------
              D = 1.3224, J1 = 3.41, J2 = 0.0, J3a = 0.0,     J3b = 0.4, R = -20,          # + # H_J' Hallen 2022
              # D = 1.3224; J1 = 3.41; J2 = 0.0; J3a = -0.00466; J3b = 0.4           # + # H_para
              # D = 0.0; J1 = -5.7; J2 = 0.0; J3a = 0.0; J3b = 0.0                 # + # H_NN Hallen 2022
              # D = 1.3224; J1 = 0.0; J2 = 0.0; J3a = 0.0; J3b = 0.0               # + # Test
                hfields = 0):
    check_eta(Nx, Ny, Nz,
              D = D, J1 = J1, J2 = J2, J3a = J3a,     J3b = J3b, R = R,
              verbose = verbose)
    time0 = time.time()
    energies = np.zeros(N_steps+1)
    times = np.zeros(N_steps+1)
    states = np.zeros((N_steps+1, Nx, Ny, Nz, Nmu), dtype=int8)
    states[0,:,:,:,:] = initial_state
    transverse_field = calc_effective_field(initial_state, Nx, Ny, Nz)
    labels = calc_labels(Nx, Ny, Nz)

    if use_gpu == True:
        state_gpu = cp.asarray(initial_state, dtype=int8)
        global eta_gpu
        hfields_gpu = cp.asarray(hfields, dtype=float32)
        hfield_gpu = cp.zeros(3, dtype=float32) # initialize for delta_hfield implementation
    for n in arange(0, N_steps):
        Tn = Temp[n] # this helped
        if type(hfields) != int:
            delta_hfield_gpu = hfields_gpu[n,:] - hfield_gpu
            hfield_gpu = hfields_gpu[n,:]
        if not np.all(delta_hfield_gpu == 0):  # if non-zero change in external field
            add_zeeman_gpu(delta_hfield_gpu, Nx,Ny,Nz) # add zeeman term to eta_gpu global
        
        if flip_mode == "two-rate":
            if use_gpu == True:
                rates_n = flatten_state(rates_mask(cp.asnumpy(calc_rates_gpu(state_gpu, Tn, Nx, Ny, Nz)), states[n, :,:,:,:], h = transverse_field), Nx, Ny, Nz)
            else:
                rates_n = flatten_state(rates_mask(calc_rates(states[n, :,:,:,:], Tn, Nx, Ny, Nz), states[n, :,:,:,:],
                                                   h = transverse_field), Nx, Ny, Nz)
        elif flip_mode == "one-rate":
            if use_gpu == True:
                rates_n = flatten_state(cp.asnumpy(calc_rates_gpu(state_gpu, Tn, Nx, Ny, Nz)), Nx, Ny, Nz)
            else:
                rates_n = flatten_state(calc_rates(states[n, :,:,:,:], Tn, Nx, Ny, Nz), Nx, Ny, Nz)
        
        _ = np.argsort(rates_n) # this helped
        rates_n = np.cumsum(rates_n[_])
        labels_ = labels[_]
        Qk = rates_n[-1]
        
        loc_flip = nutoi( labels_[ np.where((rates_n >= random.rand()*Qk) )[0][0] ], Nx, Ny, Nz)
        x, y, z, mu = loc_flip
        if use_gpu == True:
            E_dif = calc_energy_dif_gpu(x, y, z, mu, state_gpu)
        else:
            E_dif = calc_energy_dif(x, y, z, mu, states[n,:,:,:,:], Nx, Ny, Nz)
        transverse_field = update_effective_field(transverse_field, states[n,:,:,:,:], x,y,z,mu)
        
        dt = ( (Nmu*Nx*Ny*Nz)*1.1e-4*np.exp(6*0.26/(Tn-0.26)) # parameterize time units
               *(Qk**-1)*log(1/random.rand()) )
        times[n+1] = times[n] + dt
        
        states[n+1,:,:,:,:] = states[n,:,:,:,:]
        states[n+1,x,y,z,mu] *= -1
        if use_gpu == True:
            state_gpu[x,y,z,mu] *= -1
        energies[n+1] = energies[n] + E_dif
        
        # if not np.all(hfield_gpu == 0):  # if non-zero external field
        #     add_zeeman_gpu(-hfield_gpu, Nx,Ny,Nz) # subtract zeeman term to eta_gpu global
            
    return states, energies, times