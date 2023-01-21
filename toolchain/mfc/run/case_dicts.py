from ..      import common
from ..state import ARG


COMMON = [
    "parallel_io",
    "case_dir", 
    "m", "adv_alphan", "num_fluids", "model_eqns",
    "weno_order", "n", "precision"
]


PRE_PROCESS = COMMON + [
    'old_grid', 'old_ic', 't_step_old',
    'fluid_rho', 'num_patches'
]

for cmp in ["x", "y"]:
    for prepend in ["domain%beg", "domain%end", "a", "b"]:
        PRE_PROCESS.append(f"{cmp}_{prepend}")

    for append in ["stretch", "a", "loops"]:
        PRE_PROCESS.append(f"{append}_{cmp}")

    PRE_PROCESS.append(f"bc_{cmp}%beg")
    PRE_PROCESS.append(f"bc_{cmp}%end")

for f_id in range(1, 10+1):
    PRE_PROCESS.append(f'fluid_rho({f_id})')

    for attribute in ["gamma", "pi_inf"]:
        PRE_PROCESS.append(f"fluid_pp({f_id})%{attribute}")

for p_id in range(1, 10+1):
    for attribute in ["geometry", "radius", "radii", "epsilon", "beta",
                      "normal", "smoothen", "smooth_patch_id", "alpha_rho",
                      "smooth_coeff", "rho", "vel", "pres", "alpha", "gamma",
                      "pi_inf"]:
        PRE_PROCESS.append(f"patch_icpp({p_id})%{attribute}")

    for cmp_id, cmp in enumerate(["x", "y"]):
        cmp_id += 1
        PRE_PROCESS.append(f'patch_icpp({p_id})%{cmp}_centroid')
        PRE_PROCESS.append(f'patch_icpp({p_id})%length_{cmp}')

        for append in ["radii", "normal", "vel"]:
            PRE_PROCESS.append(f'patch_icpp({p_id})%{append}({cmp_id})')

    for arho_id in range(1, 10+1):
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha({arho_id})')
        PRE_PROCESS.append(f'patch_icpp({p_id})%alpha_rho({arho_id})')

    if p_id >= 2:
        PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch')

        for alter_id in range(1, p_id):
            PRE_PROCESS.append(f'patch_icpp({p_id})%alter_patch({alter_id})')


SIMULATION = COMMON + [
    'run_time_info', 't_step_old', 'dt', 't_step_start',
    't_step_stop', 't_step_save', 'time_stepper', 'weno_eps',
    'weno_avg', 'weno_Re_flux', 'fd_order', 'num_probes', 'probe_wrt', 'cu_mpi'
]

for cmp in ["x", "y"]:
    SIMULATION.append(f'bc_{cmp}%beg')
    SIMULATION.append(f'bc_{cmp}%end')

for wrt_id in range(1,10+1):
    for cmp in ["x", "y"]:
        SIMULATION.append(f'probe_wrt({wrt_id})%{cmp}')

for probe_id in range(1,3+1):
    for cmp in ["x", "y"]:
        SIMULATION.append(f'probe({probe_id})%{cmp}')

for f_id in range(1,10+1):
    for attribute in ["gamma", "pi_inf"]:
        SIMULATION.append(f"fluid_pp({f_id})%{attribute}")

    for re_id in [1, 2]:
        SIMULATION.append(f"fluid_pp({f_id})%Re({re_id})")

POST_PROCESS = COMMON + [
    't_step_start', 't_step_stop', 't_step_save', 
    'format', 'coarsen_silo', 'fourier_decomp',
    'fourier_modes%beg', 'fourier_modes%end', 'alpha_rho_wrt', 'rho_wrt',
    'mom_wrt', 'vel_wrt', 'flux_lim', 'flux_wrt', 'E_wrt', 'pres_wrt',
    'alpha_wrt', 'kappa_wrt', 'gamma_wrt', 'heat_ratio_wrt', 'pi_inf_wrt',
    'pres_inf_wrt', 'cons_vars_wrt', 'prim_vars_wrt', 'c_wrt', 'omega_wrt',
    'schlieren_wrt', 'schlieren_alpha', 'fd_order'
]

for cmp_id in range(1,2+1):
    cmp = ["x", "y"][cmp_id-1]

    POST_PROCESS.append(f'bc_{cmp}%beg')
    POST_PROCESS.append(f'bc_{cmp}%end')

    for attribute in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
        POST_PROCESS.append(f'{attribute}({cmp_id})')

for fl_id in range(1,10+1):
    for append in ["schlieren_alpha", "alpha_rho_wrt", "alpha_wrt", "kappa_wrt"]:
        POST_PROCESS.append(f'{append}({fl_id})')

    for attribute in ["gamma", "pi_inf"]:
        POST_PROCESS.append(f"fluid_pp({fl_id})%{attribute}")


CASE_OPTIMIZATION = [ "weno_order" ]


def get_input_dict_keys(target_name: str) -> list:
    result = None
    if target_name == "pre_process":  result = PRE_PROCESS.copy()
    if target_name == "simulation":   result = SIMULATION.copy()
    if target_name == "post_process": result = POST_PROCESS.copy()

    if result == None:
        raise common.MFCException(f"[INPUT DICTS] Target {target_name} doesn't have an input dict.")

    if not ARG("case_optimization") or target_name != "simulation":
        return result
    
    return [ x for x in result if x not in CASE_OPTIMIZATION ]
