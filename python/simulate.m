function sim_file=simulate(data_file, prefix, sim_vertices, sim_woi,...
    sim_signals, dipole_orientations, dipole_moments, sim_patch_sizes, SNR)

addpath('/home/bonaiuto/DANC_spm12/spm12');

% Start SPM
spm('defaults','eeg');
spm_jobman('initcfg');

Dmesh=spm_eeg_load(data_file);

%% get location to simulate dipole on this mesh
simpos=Dmesh.inv{1}.mesh.tess_mni.vert(sim_vertices,:); 

% Make signal have unit variance
sim_signals=sim_signals./repmat(std(sim_signals'),size(sim_signals,2),1)';

[Dnew,~]=spm_eeg_simulate({Dmesh}, 'prefix', prefix,...
    'patchmni', simpos, 'simsignal', sim_signals,...
    'ormni', dipole_orientations, 'woi', sim_woi, 'SNRdB', SNR,...
    'dipfwhm', sim_patch_sizes, 'nAmdipmom', dipole_moments);
[a1, ~, ~]=fileparts(data_file);
sim_file=fullfile(a1, Dnew.fname);