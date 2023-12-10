import os
import matlab.engine
import numpy as np

from surf import smoothmesh_multilayer_mm


def invert_ebb(subj_id, out_dir, nas, lpa, rpa, mri_fname, mesh_fname, data_file, n_layers, patch_size=5,
               n_temp_modes=4, woi=[-np.inf, np.inf]):
    subject_out_dir = os.path.join(out_dir, subj_id, 'inv')

    print(f'Smoothing {mesh_fname}')
    _ = smoothmesh_multilayer_mm(mesh_fname, patch_size, n_layers, n_jobs=-1)

    parasite = matlab.engine.start_matlab()
    mesh_base=os.path.split(os.path.splitext(mesh_fname)[0])[-1]
    data_fname=os.path.split(data_file)[-1]
    coreg_fname = os.path.join(subject_out_dir, f'{mesh_base}.{data_fname}')

    F=parasite.invert_ebb(data_file, coreg_fname, mri_fname, mesh_fname, matlab.double(nas), matlab.double(lpa),
                          matlab.double(rpa), float(patch_size), float(n_temp_modes), matlab.double(woi), nargout=1)
    parasite.close()

    return [coreg_fname, F]
