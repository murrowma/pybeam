# ---------------------------------------------------------------------- #
# Import Python packages and C libraries.
# ---------------------------------------------------------------------- #

import cython
cimport cython

from libc.math cimport pow, exp, log, sin, cos
from libc.stdlib cimport srand, rand, RAND_MAX

# ---------------------------------------------------------------------- #
# Model parameter functions.
# ---------------------------------------------------------------------- #

# total number of input parameters
DEF N_phi = 9

# function for the non-decision time
cdef double non_decision_time(double phi[N_phi]):
    return phi[0]

# function for the relative start point
cdef double relative_start(double phi[N_phi]):
    return phi[1]

# function for the drift rate
cdef double drift(double phi[N_phi], double x, double t):
    return phi[2]

# function for the diffusion rate
cdef double diffusion(double phi[N_phi], double x, double t):
    return phi[3]

# function for the upper decision threshold
cdef double upper_decision_threshold(double phi[N_phi], double t):
    return phi[4]*exp(-t/phi[5])

# function for the lower decision threshold
cdef double lower_decision_threshold(double phi[N_phi], double t):
    return -upper_decision_threshold(phi, t)

# function for the contamination strength
cdef double contamination_strength(double phi[N_phi]):
    return phi[6]

# function for the contamination probability
cdef double contamination_probability(double phi[N_phi], double t):
    if (t < phi[7]) or (t > phi[8]):
        return 0.0
    else:
        return 1.0/(phi[8] - phi[7])
        
# ---------------------------------------------------------------------- #
# Function used to modify time step, generally not necessary to modify.
# ---------------------------------------------------------------------- #

# function to modify the time step
cdef double modify_dt(double phi[N_phi], double t):
    return 1.0
