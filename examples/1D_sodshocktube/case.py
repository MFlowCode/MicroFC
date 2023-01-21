#!/usr/bin/env python3

import math
import json

# Numerical setup
Nx = 399
dx = 1./(1.*(Nx+1))

Tend, Nt = 0.1, 1000
mydt     = Tend/(1.*Nt)

# Configuring case dictionary
print(json.dumps({
    # Logistics ================================================================
    'case_dir'                     : '\'.\'',
    'run_time_info'                : 'T',
    # ==========================================================================

    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : 0.E+00,
    'x_domain%end'                 : 1.E+00,
    'm'                            : Nx,
    'n'                            : 0,
    'dt'                           : mydt,
    't_step_start'                 : 0,
    't_step_stop'                  : int(Nt),
    't_step_save'                  : int(math.ceil(Nt/10.)),
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'num_fluids'                   : 1,
    'time_stepper'                 : 3,
    'weno_vars'                    : 2,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -3,
    'bc_x%end'                     : -3,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'T',
    # ==========================================================================
                                                            
    # Patch 1 L ================================================================
    'patch_icpp(1)%geometry'       : 1,
    'patch_icpp(1)%x_centroid'     : 0.25,
    'patch_icpp(1)%length_x'       : 0.5,
    'patch_icpp(1)%vel(1)'         : 0.0,
    'patch_icpp(1)%pres'           : 1.0,
    'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
    'patch_icpp(1)%alpha(1)'       : 1.,
    # ==========================================================================

    # Patch 2 R ================================================================
    'patch_icpp(2)%geometry'       : 1,
    'patch_icpp(2)%x_centroid'     : 0.75,
    'patch_icpp(2)%length_x'       : 0.5,
    'patch_icpp(2)%vel(1)'         : 0.0,
    'patch_icpp(2)%pres'           : 0.1,
    'patch_icpp(2)%alpha_rho(1)'   : 0.125E+00,
    'patch_icpp(2)%alpha(1)'       : 1.,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
    'fluid_pp(1)%pi_inf'           : 0.0,
    # ==========================================================================
}))

# ==============================================================================
