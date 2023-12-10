import matlab.engine
import numpy as np

def run_sinusoidal_simulation(data_file, mri_fname, mesh_fname, nas, lpa, rpa, sim_vertex, freq, dipole_moment,
                              sim_patch_size, SNR, invwoi=[-np.inf, np.inf]):
    parasite = matlab.engine.start_matlab()
    sim_fname=parasite.simulate_sinusoidal_patch(data_file, mri_fname, mesh_fname, nas, lpa, rpa, sim_vertex, invwoi,
                                                 float(freq), float(dipole_moment), float(sim_patch_size), float(SNR),
                                                 nargout=1)
    parasite.close()

    return sim_fname