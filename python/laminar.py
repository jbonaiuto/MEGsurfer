import numpy as np

from invert import invert_ebb


def model_comparison(subj_id, out_dir, nas, lpa, rpa, mri_fname, mesh_fnames, sim_fname, patch_size=5,
                     n_temp_modes=4, woi=[-np.inf, np.inf]):
    f_vals=[]
    for mesh_fname in mesh_fnames:
        [_, f_val] = invert_ebb(subj_id, out_dir, nas, lpa, rpa, mri_fname, mesh_fname, sim_fname, 1,
                                patch_size=patch_size, n_temp_modes=n_temp_modes, woi=woi)
        f_vals.append(f_val)
    return f_vals