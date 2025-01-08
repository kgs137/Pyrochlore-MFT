global Nmu, im, identity, sigmax, sigmay, sigmaz, local_axis
import numpy as np
Nmu = 4
im = complex(0,1)
identity = np.array([[1.0, 0.0],
                     [0.0, 1.0]], dtype=np.complex64)
sigmax = np.array([[0.0, 1.0],
                   [1.0, 0.0]], dtype=np.complex64)
sigmay = np.array([[0.0, -1*im],
                   [1*im, 0.0]], dtype=np.complex64)
sigmaz = np.array([[1.0, 0.0],
                   [0.0,-1.0]], dtype=np.complex64)
u0 = (3**-0.5)*np.array([-1.0, -1.0, -1.0]); u1 = (3**-0.5)*np.array([-1.0, 1.0, 1.0]);
u2 = (3**-0.5)*np.array([1.0, -1.0, 1.0]);   u3 = (3**-0.5)*np.array([1.0, 1.0, -1.0]);
local_axis = np.array([u0, u1, u2, u3])
