# ---------------------------------------------------------------------- #
# Import Python packages and C libraries.
# ---------------------------------------------------------------------- #

import cython
cimport cython

from libc.math cimport pow, exp, log, sin, cos

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
    cdef double mu
    cdef double l
    cdef double k
    cdef double t0
    cdef double t1

    mu = phi[2]
    l = phi[3]
    k = phi[4]
    t0 = phi[7]
    t1 = phi[8]

    if (t < t0):
        return mu*(1.0 + k*t) + ( k/(1.0 + k*t) - l )*x
    elif (t >= t0) and (t < t1):
        return -mu*(1.0 + k*t) + ( k/(1.0 + k*t) - l )*x
    else:
        return mu*(1.0 + k*t) + ( k/(1.0 + k*t) - l )*x

# function for the diffusion rate
cdef double diffusion(double phi[N_phi], double x, double t):
    cdef double sigma
    cdef double k

    sigma = phi[5]
    k = phi[4]

    return phi[5]*(1.0 + phi[4]*t)

# function for the upper decision threshold
cdef double upper_decision_threshold(double phi[N_phi], double t):
    return phi[6]

# function for the lower decision threshold
cdef double lower_decision_threshold(double phi[N_phi], double t):
    return -upper_decision_threshold(phi, t)

# function for the contamination strength
cdef double contamination_strength(double phi[N_phi]):
    return 0.0

# function for the contamination probability
cdef double contamination_probability(double phi[N_phi], double t):
    return 0.0
        
# ---------------------------------------------------------------------- #
# Function to modify time step with unusual likelihood functions.
# ---------------------------------------------------------------------- #

# function to modify the time step
cdef double modify_dt(double phi[N_phi], double t):
    return 1.0
