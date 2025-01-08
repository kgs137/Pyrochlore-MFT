from numpy import *
import numpy as np

# Bravais Lattice Vectors
a = 1.0
l0 = array([0.0, 0.0, 0.0]);     l1 = array([0, 0.5*a, 0.5*a]);
l2 = array([0.5*a, 0.0, 0.5*a]); l3 = array([0.5*a, 0.5*a, 0.0]);
l = array([l0, l1, l2, l3])

def calc_labels(start = 0, stop = 100):
    Nt = stop - start
    labels = Nt*np.ones(Nt)

    label = 0
    for d in arange(start, stop):
        if labels[d-start] == Nt:
            labels[d-start] = label
            
            for l in d+np.where(energies[d:stop] == energies[d])[0]:
                if np.array_equiv(states[d-start, :,:,:,:], states[l-start, :,:,:,:]):
                    labels[l-start] = label

            label += 1
    return labels

def find_monopoles(state, Nx, Ny, Nz):
    lattice = zeros((Nx,Ny,Nz,3), dtype=float)
    for x in arange(0, Nx):
        for y in arange(0, Ny):
            for z in arange(0, Nz):
                lattice[x, y, z, :]   = x*l[1] + y*l[2] + z*l[3]
    
    charges = []
    locations = []
    index = []
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                # Explicit Tetrahedra
                net = sum(state[i,j,k,:].reshape(-1))
                if net != 0:
                    charges.append(net)
                    locations.append(lattice[i,j,k,:] + array([1/8, 1/8, 1/8]))
                    index.append(np.array([i, j, k]))
                
                # Implicit Tetrahedra
                net = -1*(state[i,j,k,0] + state[i-1,j,k,1] + state[i,j-1,k,2] + state[i,j,k-1,3])
                if net != 0:
                    charges.append(net)
                    locations.append((lattice[i,j,k,:] + lattice[i-1,j,k,:] + lattice[i,j-1,k,:] + lattice[i,j,k-1,:])/4.0 + array([1/8, 1/8, 1/8]))
                    index.append(-1.0*np.array([i, j, k]))
                
    return charges, np.array(locations), index

def extract_paths(states, Nx, Ny, Nz):
    # list initialization
    poles = []
    paths = []
    active = [False] # skip initial monopoles
    index = []; ind = 0;
    state_index = []; location_index = [];
    c_pair = []; a_pair = [];
    first_active = 0
    Nd = states.shape[0]; Nd_inc = int(round( Nd * 0.01 ))

    # increment through states
    for n in range(0, Nd):
        if n % Nd_inc == 0:
            print(n/Nd, first_active, len(active))
        charges, locations, ijks = find_monopoles(states[n, :,:,:,:], Nx, Ny, Nz)
        if n > 5*Nd_inc:
            try:
                first_active = where(array(active)==True)[0][0] #argsame(True, array(active[first_active:])) # Recently modified
            except:
                first_active = where(array(active)==False)[0][-1]
        else:
            first_active = 0
        active_check = np.zeros(len(active)+len(charges)) #.................................

        # iterate over monopoles
        for q, x, ijk in zip(charges, locations, ijks):
            added_to_existing = False
            
            # Search existing paths (standing)
            for q0, x0, i0 in zip(poles[first_active:], paths[first_active:], index[first_active:]):
                x0 = x0[-1]
                if q0 == q and np.array_equiv(x0,x) and active[i0] == True:
                    paths[i0].append(x)
                    state_index[i0].append(n)
                    location_index[i0].append(ijk)
                    active[i0] = True
                    added_to_existing = True
                    break
            
            # Search existing paths (adjacents)
            if added_to_existing == False:
                for q0, x0, i0 in zip(poles[first_active:], paths[first_active:], index[first_active:]):
                    x0 = x0[-1]
                    if q0 == q and dot(x0-x,x0-x)<=0.5 and active[i0] == True:
                        paths[i0].append(x)
                        state_index[i0].append(n)
                        location_index[i0].append(ijk)
                        active[i0] = True
                        added_to_existing = True
                        break
                    else:
                        active_check[i0] += 1 #............................................
    
            # Check for new paths
            if added_to_existing == False:
                poles.append(q)
                paths.append([x])
                state_index.append([n])
                location_index.append([ijk])
                if ind != 0:
                    active.append(True)
                index.append(ind)
    #             path_labels.append(labels[n])
                # Create pair
                for q0, x0, i0 in zip(poles[first_active:], paths[first_active:], index[first_active:]):
                    x0 = x0[-1]
                    if q0 == -q and dot(x0-x,x0-x)<=0.5 and active[i0] == True:
                        c_pair.append([i0, ind])
                        break
                ind += 1
                
        for i in range(len(active)):
            if active_check[i] >= len(charges): #.........................................
                # Annihilate pair
    #             for q0, x0, i0 in zip(poles[first_active:], paths[first_active:], index[first_active:]):
    #                 x0 = x0[-1]
    #                 if q0 == -q and dot(x0-x,x0-x)<=0.5 and active[i0] == True:
    #                     a_pair.append([i0, i])
    #                     break
                active[i] = False
    return paths