# ---------------------------------------------------------------------- #
# Import Python packages and C libraries.
# ---------------------------------------------------------------------- #

#!python
#cython: language_level=3

import cython
cimport cython
import numpy as np

from libc.math cimport ceil, floor, log, sqrt, sin, cos, pi
from libc.stdlib cimport srand, rand, RAND_MAX

include "model.pyx"

# ---------------------------------------------------------------------- #
# solveFP
# ---------------------------------------------------------------------- #

DEF N_deps_max = 501; # max number of spatial mesh points
DEF N_dt = 2500; # maximum number of time steps

DEF N_dt_ini = 25; # number of small initial time steps to dampen Crank-Nicolson oscillations
DEF N_dt_scaleup = 100;
DEF dt_ini_scale = 0.01; # sets initial time step size (dt_ini = dt_ini_scale * dt)

DEF ds_ratio_cutoff = 0.02; # accept time step size if boundary collapse at or beneath this ratio
DEF dt_mod_scale = 0.5; # if boundary collapses too much during step, decrease by this factor

DEF p_fpt_cutoff = 1.0e-5; # stop running solver if fpt drops below this threshold
DEF int_prob_cutoff = 0.25; # make sure integrated probability at least equals this before cutting off (makes sure we don't stop initially when probability equals zero)

DEF threshold_cutoff = 1.0e-4; # stop threshold if it reaches threshold_cutoff

DEF delt = 1.0e-5; # time step size for numerical boundary derivatives

@cython.cdivision(True)
cpdef solveFP(phi, N_deps, dt_scale, fpt_bu, fpt_bl, output_loglike):

    # ----------------------------------- #
    # Variable/array definitions in order of appearance.
    # ----------------------------------- #
    
    # -------------------- #
    cdef double phic[N_phi] # Cython array containing model parameters
    cdef double dt_scalec # Cython variable for input dt_scale 
    # -------------------- #
    cdef double p_n[N_deps_max] # probability at old time step n
    cdef double p_ncn[N_deps_max] # modified probability at old time step n (for CN scheme)
    cdef double p_np1[N_deps_max] # probability at current time step n+1
    cdef double p_fpt[3][N_dt] # first passage time probability
    # -------------------- #
    cdef double s_ini; # initial boundary separation
    cdef double deps; # spatial step size
    cdef double eps[N_deps_max]; # spatial coordinates
    cdef double dt_base, dt_ini, dt;
    # -------------------- #
    cdef double ww;
    cdef double w_dist[2];
    # -------------------- #
    cdef double t_nd;
    cdef double tt;
    cdef double gg;
    cdef double int_prob;
    cdef double p_fpt_max;
    cdef double bu_n, dbudt_n;
    cdef double bl_n, dbldt_n;
    cdef double s_n, dsdt_n;
    cdef double t_in;
    cdef double bu_np1, dbudt_np1;
    cdef double bl_np1, dbldt_np1;
    cdef double s_np1, dsdt_np1;
    # -------------------- #
    cdef double bu_nm1, bl_nm1;
    cdef double s_nm1;
    cdef double ds_ratio, ds_ratio_np1, ds_ratio_nm1;
    # -------------------- #
    cdef double v_n[N_deps_max];
    cdef double v_np1[N_deps_max];
    cdef double sigma_n[N_deps_max];
    cdef double sigma_np1[N_deps_max];
    cdef double sigma2_n[N_deps_max];
    cdef double sigma2_np1[N_deps_max];
    cdef double x_n;
    cdef double x_np1;
    cdef double alpha_n, alpha_np1;
    cdef double beta_n, beta_np1;
    cdef double gamma_n, gamma_np1;
    cdef double AA[N_deps_max];
    cdef double BB[N_deps_max];
    cdef double CC[N_deps_max];
    cdef double DD[N_deps_max];
    cdef double EE[N_deps_max];
    cdef double FF[N_deps_max];
    cdef double WW;
    cdef double g_prob;
    # -------------------- #
    cdef double loglike;
    cdef double fpt_buc;
    cdef double fpt_blc;    
    # -------------------- #
    cdef int ii, jj, kk
    cdef int N_depsc
    cdef int w_ind[2]
    cdef int n_cut
    # -------------------- #
    cdef double max_rtbu
    cdef double max_rtbl
    cdef double max_rt
    
    # ----------------------------------- #
    # Move Python input variables/arrays into Cython variables/arrays for speed.
    # ----------------------------------- #

    max_rtbu = 0.0
    max_rtbl = 0.0
    if (len(fpt_bu) > 0):
        max_rtbu = np.max(fpt_bu)
    if (len(fpt_bl) > 0):
        max_rtbl = np.max(fpt_bl)
    
    if (max_rtbu > max_rtbl):
        max_rt = max_rtbu
    else:
        max_rt = max_rtbl

    if (max_rt == 0.0):
        max_rt = 100.0

    for ii in range(N_phi):
        phic[ii] = phi[ii]
    N_depsc = N_deps
    dt_scalec = dt_scale
    
    # ----------------------------------- #
    # Initialize arrays.
    # ----------------------------------- #
    
    for ii in range(N_deps_max):
        p_n[ii] = 0.0;
        p_ncn[ii] = 0.0;
        p_np1[ii] = 0.0;
        
    for ii in range(N_dt):
        p_fpt[0][ii] = 0.0;
        p_fpt[1][ii] = 0.0;
        p_fpt[2][ii] = 0.0;
        
    # ---------------------------------------- #
    # Set spatial and time step sizes.
    # ---------------------------------------- #
    
    # determine initial boundary separation
    s_ini = upper_decision_threshold(phic, 0.0) - lower_decision_threshold(phic, 0.0)
    
    # calculate spatial step size
    deps = 1.0/(N_depsc - 1.0)
    
    # create array of spatial locations
    eps[0] = 0.0;
    for ii in range(1, N_depsc):
        eps[ii] = eps[ii-1] + deps
        
    # calculate time step size
    dt_base = approx_dt(phic, dt_scalec)
    dt_ini = dt_ini_scale*dt_base
    dt = dt_ini
    
    # ---------------------------------------- #
    # Initialize probability array at start point.
    # ---------------------------------------- #
    
    # determine relative start point
    ww = relative_start(phic)
    
    # determine which array element the walker starts at
    w_ind[0] = int( ceil(ww/deps) )
    w_ind[1] = int( floor(ww/deps) )
    
    # determine distance between relative start point and nearest array elements
    w_dist[0] = abs( ww - deps*w_ind[0] )
    w_dist[1] = abs( ww - deps*w_ind[1] )
    
    # numerically approximate delta function, distribute between nearest array elements
    if w_ind[0] == w_ind[1]:
        p_n[ w_ind[0] ] = 1.0/deps
    else:
        p_n[ w_ind[0] ] = (1.0 - w_dist[0]/deps)/deps
        p_n[ w_ind[1] ] = (1.0 - w_dist[1]/deps)/deps
        
    # ---------------------------------------- #
    # Remaining program initilization.
    # ---------------------------------------- #
    
    t_nd = non_decision_time(phic);
    tt = 0.0;
    gg = contamination_strength(phic);
    int_prob = 0.0;
    p_fpt_max = 0.0;
    n_cut = N_dt;
            
    bu_n = 0.0;
    dbudt_n = 0.0;    
    bl_n = 0.0;
    dbldt_n = 0.0;    
    s_n = 0.0;
    dsdt_n = 0.0
    
    t_in = dt_ini/10.0;    
    bu_np1 = upper_decision_threshold(phic, t_in);
    dbudt_np1 = (upper_decision_threshold(phic, t_in + delt) - bu_np1)/delt;    
    
    bl_np1 = lower_decision_threshold(phic, t_in);
    dbldt_np1 = (lower_decision_threshold(phic, t_in + delt) - bl_np1)/delt;    
    
    s_np1 = bu_np1 - bl_np1;
    dsdt_np1 = dbudt_np1 - dbldt_np1;
    
    # ---------------------------------------- #
    # Solve PDE using Crank-Nicolson.
    # ---------------------------------------- #
    
    for ii in range(N_dt):
        
        # ---------------------------------------- #
        # Set old boundary location.
        # ---------------------------------------- #
        
        bu_n = bu_np1;
        dbudt_n = dbudt_np1;
        bl_n = bl_np1;
        dbldt_n = dbldt_np1;
            
        s_n = bu_n - bl_n;
        dsdt_n = dbudt_n - dbldt_n;
        
        # ---------------------------------------- #
        # Modify time step for moving boundaries.
        # ---------------------------------------- #
    
        if (ii < N_dt_ini):
            dt = dt_ini
        elif (ii >= N_dt_ini) & (ii < N_dt_ini + N_dt_scaleup):
            dt += (dt_base - dt_ini)/N_dt_scaleup
        else:
            dt = dt_base;
    
        jj = 0
        kk = 0;
        while (jj < 1) & (kk < 10):
            t_in = tt + dt;
            bu_np1 = upper_decision_threshold(phic, t_in);                
            bl_np1 = lower_decision_threshold(phic, t_in);

            t_in = tt - dt;
            if (t_in <= 0.0):
                bu_nm1 = upper_decision_threshold(phic, dt_ini/10.0);
                bl_nm1 = lower_decision_threshold(phic, dt_ini/10.0);
            else:
                bu_nm1 = upper_decision_threshold(phic, t_in);
                bl_nm1 = lower_decision_threshold(phic, t_in);

            s_np1 = bu_np1 - bl_np1;
            s_nm1 = bu_nm1 - bl_nm1;

            ds_ratio = 0.0;
            ds_ratio_np1 = abs(s_n - s_np1)/s_n;
            ds_ratio_nm1 = abs(s_n - s_nm1)/s_n;

            if (ds_ratio_np1 < ds_ratio_nm1):
                ds_ratio = ds_ratio_nm1;
            else:
                ds_ratio = ds_ratio_np1;

            if (ds_ratio > ds_ratio_cutoff):
                dt = dt*dt_mod_scale;
            else:
                jj = 1;

            kk += 1;
            
        # ---------------------------------------- #
        # Calculate new boundary location.
        # ---------------------------------------- #
    
        dt = modify_dt(phic, tt)*dt
        tt += dt;
        
        bu_np1 = upper_decision_threshold(phic, tt);
        dbudt_np1 = (upper_decision_threshold(phic, tt+delt) - bu_np1)/delt;

        if (bu_np1 <= threshold_cutoff):
            bu_np1 = threshold_cutoff
            dbudt_np1 = 0.0
        
        bl_np1 = lower_decision_threshold(phic, tt);
        dbldt_np1 = (lower_decision_threshold(phic, tt+delt) - bl_np1)/delt;

        if (bl_np1 >= -threshold_cutoff):
            bl_np1 = -threshold_cutoff
            dbldt_np1 = 0.0
        
        s_np1 = bu_np1 - bl_np1;
        dsdt_np1 = dbudt_np1 - dbldt_np1;
        
        # ---------------------------------------- #
        # Invert matrix using tridiagonal matrix algorithm.
        # ---------------------------------------- #
        
        for jj in range(2):
            
            p_np1[jj] = 0.0;

            if (ii == 0):
                x_n = s_n*eps[jj] + bl_n
                v_n[jj] = drift(phic, x_n, tt-dt);
                sigma_n[jj] = diffusion(phic, x_n, tt-dt);
                sigma2_n[jj] = sigma_n[jj]*sigma_n[jj];
            else:
                x_n = x_np1;
                v_n[jj] = v_np1[jj];
                sigma_n[jj] = sigma_np1[jj];
                sigma2_n[jj] = sigma2_np1[jj];

            x_np1 = s_np1*eps[jj] + bl_np1
            v_np1[jj] = drift(phic, x_np1, tt);        
            sigma_np1[jj] = diffusion(phic, x_np1, tt);
            sigma2_np1[jj] = sigma_np1[jj]*sigma_np1[jj];
        
        for jj in range(1, N_depsc-1):

            p_np1[jj+1] = 0.0;

            if (ii == 0):
                x_n = s_n*eps[jj+1] + bl_n
                v_n[jj+1] = drift(phic, x_n, tt-dt);
                sigma_n[jj+1] = diffusion(phic, x_n, tt-dt);
                sigma2_n[jj+1] = sigma_n[jj+1]*sigma_n[jj+1];
            else:
                x_n = x_np1;
                v_n[jj+1] = v_np1[jj+1];
                sigma_n[jj+1] = sigma_np1[jj+1]
                sigma2_n[jj+1] = sigma2_np1[jj+1]

            x_np1 = s_np1*eps[jj+1] + bl_np1
            v_np1[jj+1] = drift(phic, x_np1, tt);        
            sigma_np1[jj+1] = diffusion(phic, x_np1, tt);
            sigma2_np1[jj+1] = sigma_np1[jj+1]*sigma_np1[jj+1];
            
            alpha_n = dt/(4.0*s_n*deps);
            alpha_np1 = dt/(4.0*s_np1*deps);
            
            beta_n = eps[jj]*dsdt_n + dbldt_n;
            beta_np1 = eps[jj]*dsdt_np1 + dbldt_np1;
            
            gamma_n = dt/(4.0*s_n*s_n*deps*deps);
            gamma_np1 = dt/(4.0*s_np1*s_np1*deps*deps);
            
            AA[jj] = alpha_np1*beta_np1 - alpha_np1*v_np1[jj-1] - gamma_np1*sigma2_np1[jj-1];
            BB[jj] = 1.0 + 2.0*gamma_np1*sigma2_np1[jj];
            CC[jj] = -alpha_np1*beta_np1 + alpha_np1*v_np1[jj+1] - gamma_np1*sigma2_np1[jj+1];
            
            DD[jj] = -alpha_n*beta_n + alpha_n*v_n[jj-1] + gamma_n*sigma2_n[jj-1];
            EE[jj] = 1.0 - 2.0*gamma_n*sigma2_n[jj];
            FF[jj] = alpha_n*beta_n - alpha_n*v_n[jj+1] + gamma_n*sigma2_n[jj+1];
            
            p_ncn[jj] = DD[jj]*p_n[jj-1] + EE[jj]*p_n[jj] + FF[jj]*p_n[jj+1];
            
        for jj in range(2, N_depsc - 1):
            WW = AA[jj]/BB[jj-1];
            BB[jj] = BB[jj] - WW*CC[jj-1];
            p_ncn[jj] = p_ncn[jj] - WW*p_ncn[jj-1];
            
        p_np1[N_depsc-2] = p_ncn[N_depsc-2]/BB[N_depsc-2];
        p_n[N_depsc-2] = p_np1[N_depsc-2];
        
        for jj in range(N_depsc - 3, 0, -1):
            p_np1[jj] = (p_ncn[jj] - CC[jj]*p_np1[jj+1])/BB[jj];
            p_n[jj] = p_np1[jj];
            
        # ---------------------------------------- #
        # Calculate first passage times.
        # ---------------------------------------- #

        p_fpt[0][ii] = tt + t_nd;           
        p_fpt[1][ii] = abs( 0.5 * sigma2_np1[N_depsc-1] * (4.0*p_np1[N_depsc-2] - p_np1[N_depsc-3])/(2.0*deps*s_np1*s_ini) );
        p_fpt[2][ii] = abs( 0.5 * sigma2_np1[N_depsc-1] * (4.0*p_np1[1] - p_np1[2])/(2.0*deps*s_np1*s_ini) );

        if (p_fpt[1][ii] > p_fpt_max):
            p_fpt_max = p_fpt[1][ii];
        if (p_fpt[2][ii] > p_fpt_max):
            p_fpt_max = p_fpt[2][ii];

        int_prob += p_fpt[1][ii]*dt + p_fpt[2][ii]*dt;
        g_prob = gg*contamination_probability(phic, p_fpt[0][ii])/2.0;

        if (p_fpt[0][ii] > max_rt):
            p_fpt[1][ii] = (1.0 - gg)*p_fpt[1][ii] + g_prob;
            p_fpt[2][ii] = (1.0 - gg)*p_fpt[2][ii] + g_prob;
            n_cut = ii;
            break;
        elif (int_prob > 0.25) & (p_fpt[1][ii] <= p_fpt_cutoff) and (p_fpt[2][ii] <= p_fpt_cutoff):
            p_fpt[1][ii] = (1.0 - gg)*p_fpt[1][ii] + g_prob;
            p_fpt[2][ii] = (1.0 - gg)*p_fpt[2][ii] + g_prob;
            n_cut = ii;
            break;
        else:
            p_fpt[1][ii] = (1.0 - gg)*p_fpt[1][ii] + g_prob;
            p_fpt[2][ii] = (1.0 - gg)*p_fpt[2][ii] + g_prob;

        if (p_fpt[1][ii] <= p_fpt_cutoff):
            p_fpt[1][ii] = p_fpt_cutoff
        if (p_fpt[2][ii] <= p_fpt_cutoff):
            p_fpt[2][ii] = p_fpt_cutoff
                
    # ---------------------------------------- #
    # Calculate loglikelihood of input data.
    # ---------------------------------------- #
    
    loglike = 0.0;
    
    if (output_loglike == True):
    
        for ii in range(len(fpt_bu)):

            fpt_buc = fpt_bu[ii];
            jj = n_cut;
            kk = 0;

            while (kk <= n_cut):
                if (fpt_buc < p_fpt[0][kk]):
                    jj = kk - 1;
                    kk = n_cut;
                kk += 1;

            if (fpt_buc < p_fpt[0][0]):
                g_prob = gg*contamination_probability(phic, fpt_buc)/2.0;
                if (g_prob > 0.0):
                    loglike += log(g_prob);
                else:
                    loglike += log(p_fpt[1][0]);
            elif (jj < n_cut):
                loglike += log(p_fpt[1][jj] + (fpt_buc - p_fpt[0][jj])*(p_fpt[1][jj+1] - p_fpt[1][jj])/(p_fpt[0][jj+1] - p_fpt[0][jj]));
            else:
                loglike += log(p_fpt[1][n_cut]);

        for ii in range(len(fpt_bl)):

            fpt_blc = fpt_bl[ii];
            jj = n_cut;
            kk = 0;

            while (kk <= n_cut):
                if (fpt_blc < p_fpt[0][kk]):
                    jj = kk - 1;
                    kk = n_cut;
                kk += 1;

            if (fpt_blc < p_fpt[0][0]):
                g_prob = gg*contamination_probability(phic, fpt_blc)/2.0;
                if (g_prob > 0.0):
                    loglike += log(g_prob);
                else:
                    loglike += log(p_fpt[2][0]);
            elif (jj < n_cut):
                loglike += log(p_fpt[2][jj] + (fpt_blc - p_fpt[0][jj])*(p_fpt[2][jj+1] - p_fpt[2][jj])/(p_fpt[0][jj+1] - p_fpt[0][jj]));
            else:
                loglike += log(p_fpt[2][n_cut]);

    if (output_loglike == False):
        return p_fpt;
    elif (output_loglike == True):
        return loglike;
    
# ---------------------------------------------------------------------- #
# approx_dt
# ---------------------------------------------------------------------- #
        
DEF N_w = 10; # number of accumulators simulated
DEF max_run_time = 100.0; # max run time of accumulators (seconds)
DEF sim_dt = 0.025; # multiplied by sigma to determine simulation step size (dt_xim = sim_dt_scale * sigma)
    
@cython.cdivision(True)
cdef double approx_dt(double phi[N_phi], double dt_scale):
    
    # ----------------------------------- #
    # Variable/array definitions in order of appearance.
    # ----------------------------------- #
    
    cdef double bu; # upper boundary location
    cdef double bl; # lower boundary location
    cdef double ww; # relative start point
    cdef double zz; # start point
    cdef double fpt_avg; # average first passage time
    cdef double tt; # current time
    cdef double x_n; # location of walker at previous time step
    cdef double vv;
    cdef double sigma;
    cdef double dt;
    
    cdef int ii, jj, kk, ll; # indices
    cdef int rr;
    
    # set seed for random number generator
    srand(90210);
    
    # calculate initial boundary location and start point
    bu = upper_decision_threshold(phi, 0.0);
    bl = lower_decision_threshold(phi, 0.0);
    ww = relative_start(phi);
    zz = bl + ww*(bu - bl);
    
    fpt_avg = 0.0;
    
    for ii in range(N_w):
        
        tt = 0.0
        x_n = zz;
        jj = 0;
        
        while (jj < 1) & (tt <= max_run_time):
            
            tt += sim_dt;
            vv = drift(phi, x_n, tt);
            sigma = diffusion(phi, x_n, tt);
            bu = upper_decision_threshold(phi, tt);
            bl = lower_decision_threshold(phi, tt);
            
            rr = -1 +  2*( rand() % 2 );
            
            x_np1 = x_n + sim_dt*vv + (sim_dt**0.5)*sigma*rr;
            x_n = x_np1;
            
            if (x_np1 > bu) or (x_np1 < bl):                
                fpt_avg += tt;
                jj = 1;
                
    fpt_avg = fpt_avg/N_w;
    dt = dt_scale * fpt_avg;
    
    if (dt >= 0.05):
        dt = 0.05;
    
    return dt;

# ---------------------------------------------------------------------- #
# Simulate model.
# ---------------------------------------------------------------------- #

@cython.cdivision(True)
cpdef simulate_model_c(N_sims, dt, phi, seed):
    
    # ----------------------------------- #
    # Variable/array definitions in order of appearance.
    # ----------------------------------- #
    
    cdef double dt_c; # C variable for the input time step
    cdef double phi_c[N_phi]; # C array for ddm model parameters
    cdef double bu; # upper boundary location
    cdef double bl; # lower boundary location
    cdef double ww; # relative start point
    cdef double zz; # start point
    cdef double gg;
    cdef double tt; # current time
    cdef double x_n; # location of walker at previous time step
    cdef double guess;
    cdef double t_guess
    cdef double side;
    cdef double vv; # drift rate
    cdef double sigma;
    cdef double u1
    cdef double u2
    cdef double uu
    
    cdef int ii, jj; # indices
    cdef int N_sims_c; # number of simulations
    cdef int seedc
    
    # set seed for random number generator
    seedc = seed
    srand(seedc);
    np.random.seed(seed)
    
    # load inputs into C arrays for speed
    N_sims_c = N_sims;
    dt_c = dt;
    for ii in range(N_phi):
        phi_c[ii] = phi[ii];
    
    # calculate initial boundary location and start point
    bu = upper_decision_threshold(phi_c, 0.0);
    bl = lower_decision_threshold(phi_c, 0.0);
    ww = relative_start(phi_c);
    zz = bl + ww*(bu - bl);

    # contamination strength
    gg = contamination_strength(phi_c)
    
    # allocate space for first passage time values and initialize random numbers
    fpt = np.zeros(N_sims)
    guess = 0.0
    side = 0.0
    u1 = 0.0
    u2 = 0.0
    uu = 0.0

    # approximate the max value of the guessing model distribution
    
    for ii in range(N_sims_c):
        
        tt = 0.0
        x_n = zz;
        jj = 0;         

        guess = 0.0
        side = 0.0

        guess = rand() / (RAND_MAX*1.0)

        if (guess <= gg):

            t_guess = np.random.uniform(phi[14], phi[15])
            side = rand() / (RAND_MAX*1.0)
            if (side <= 0.5):
                fpt[ii] = t_guess
            else:
                fpt[ii] = -t_guess

        else:
            
            while (jj < 1) & (tt <= max_run_time):
                
                tt += dt_c;
                vv = drift(phi_c, x_n, tt);
                sigma = diffusion(phi_c, x_n, tt);
                bu = upper_decision_threshold(phi_c, tt);
                bl = lower_decision_threshold(phi_c, tt);

                u1 = 0.0
                u2 = 0.0
                uu = 0.0

                u1 = rand() / (RAND_MAX*1.0)
                u2 = rand() / (RAND_MAX*1.0)
                uu = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
    
                x_np1 = x_n + dt_c*vv + (dt_c**0.5)*sigma*uu;
                x_n = x_np1;
                
                if (x_np1 >= bu):                
                    fpt[ii] = phi_c[0] + tt;
                    jj = 1;
                elif (x_np1 <= bl):                
                    fpt[ii] = -(phi_c[0] + tt);
                    jj = 1;
                
    return fpt;
    
# ---------------------------------------------------------------------- #
# Model check.
# ---------------------------------------------------------------------- #

def model_check(dict_length):
    if (dict_length != N_phi):
        return(False)

# ---------------------------------------------------------------------- #
# Model functions check.
# ---------------------------------------------------------------------- #

def model_functions_check(function, phi, x, t):

    cdef double phic[N_phi]
    for ii in range(N_phi):
        phic[ii] = phi[ii]

    functions_list = ['non_decision_time',
                      'relative_start',
                      'drift',
                      'diffusion',
                      'upper_decision_threshold',
                      'lower_decision_threshold',
                      'contamination_strength',
                      'contamination_probability',
                      'modify_dt']

    if (function == functions_list[0]):
        return non_decision_time(phic)
    elif (function == functions_list[1]):
        return relative_start(phic)
    elif (function == functions_list[2]):
        return drift(phic, x, t)
    elif (function == functions_list[3]):
        return diffusion(phic, x, t)
    elif (function == functions_list[4]):
        return upper_decision_threshold(phic, t)
    elif (function == functions_list[5]):
        return lower_decision_threshold(phic, t)
    elif (function == functions_list[6]):
        return contamination_strength(phic)
    elif (function == functions_list[7]):
        return contamination_probability(phic, t)
    elif (function == functions_list[8]):
        return modify_dt(phic, t)
