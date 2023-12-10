import matlab.engine
import numpy as np

def run_current_density_simulation(data_file, prefix, sim_vertices, sim_signals, dipole_moments, sim_patch_sizes, SNR,
                                   sim_woi=[-np.inf, np.inf]):
    if np.isscalar(sim_vertices):
        sim_vertices=[sim_vertices]
    if np.isscalar(dipole_moments):
        dipole_moments=[dipole_moments]
    if np.isscalar(sim_patch_sizes):
        sim_patch_sizes=[sim_patch_sizes]
    parasite = matlab.engine.start_matlab()
    sim_fname=parasite.simulate(
        data_file,
        prefix,
        matlab.double(sim_vertices),
        matlab.double(sim_woi),
        matlab.double(sim_signals.tolist()),
        matlab.double([]),
        matlab.double(dipole_moments),
        matlab.double(sim_patch_sizes),
        float(SNR),
        nargout=1
    )
    parasite.close()

    return sim_fname


def run_dipole_simulation(data_file, prefix, sim_vertices, sim_signals, dipole_orientations, dipole_moments, sim_patch_sizes,
                          SNR, sim_woi=[-np.inf, np.inf]):
    if np.isscalar(sim_vertices):
        sim_vertices=[sim_vertices]
    if np.isscalar(dipole_moments):
        dipole_moments=[dipole_moments]
    if np.isscalar(sim_patch_sizes):
        sim_patch_sizes=[sim_patch_sizes]
    parasite = matlab.engine.start_matlab()
    sim_fname=parasite.simulate(
        data_file,
        prefix,
        matlab.double(sim_vertices),
        matlab.double(sim_woi),
        matlab.double(sim_signals.tolist()),
        matlab.double(dipole_orientations),
        matlab.double(dipole_moments),
        matlab.double(sim_patch_sizes),
        float(SNR),
        nargout=1
    )
    parasite.close()

    return sim_fname