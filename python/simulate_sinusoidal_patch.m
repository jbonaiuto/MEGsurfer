function sim_file=simulate_sinusoidal_patch(data_file, mri_fname,...
    mesh_fname, nas, lpa, rpa, sim_vertex, invwoi, freq, dipole_moment,...
    sim_patch_size, SNR)

addpath('/home/bonaiuto/Dropbox/Toolboxes/DANC_spm12/spm12');

% Start SPM
spm('defaults','eeg');
spm_jobman('initcfg');

invwoi = cell2mat(invwoi);
nas=cell2mat(nas);
lpa=cell2mat(lpa);
rpa=cell2mat(rpa);

%% coregister to correct mesh
matlabbatch=[];
matlabbatch{1}.spm.meeg.source.headmodel.D = {data_file};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.mri = {[mri_fname ',1']};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.cortex = {mesh_fname};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.iskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.oskull = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.custom.scalp = {''};
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'nas';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = nas;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'lpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = lpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'rpa';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = rpa;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
spm_jobman('run', matlabbatch);

Dmesh=spm_eeg_load(data_file);

%% get location to simulate dipole on this mesh
simpos=Dmesh.inv{1}.mesh.tess_mni.vert(sim_vertex,:); 
prefix=sprintf('sim_mesh.%d.',sim_vertex);

% % Simulate source 
% matlabbatch=[];
% matlabbatch{1}.spm.meeg.source.simulate.D = {data_file};
% matlabbatch{1}.spm.meeg.source.simulate.val = 1;
% matlabbatch{1}.spm.meeg.source.simulate.prefix = prefix;
% matlabbatch{1}.spm.meeg.source.simulate.whatconditions.all = 1;
% matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.woi = invwoi;
% matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.isSin.foi = mean(invfoi);
% matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.dipmom = [dipole_moment sim_patch_size];
% matlabbatch{1}.spm.meeg.source.simulate.isinversion.setsources.locs = simpos;
% matlabbatch{1}.spm.meeg.source.simulate.isSNR.setSNR = SNR;               
% [a,b]=spm_jobman('run', matlabbatch);
% 
% sim_file=a{1}.D{1};

simsignal=sin((Dmesh.time- Dmesh.time(1))*freq*2*pi);
simsignal=simsignal./repmat(std(simsignal'),size(simsignal,2),1)';

[Dnew,meshsourceind]=spm_eeg_simulate({Dmesh}, prefix, simpos,...
    simsignal, [], invwoi, [], SNR, [], [], sim_patch_size,...
    dipole_moment, []);
sim_file=Dnew.fname;