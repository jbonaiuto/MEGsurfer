function postprocess_freesurfer_surfaces(subj_id, varargin)

% Parse inputs
defaults = struct('subjects_dir','/usr/local/freesurfer/subjects', ...
    'combine_hemispheres', true, 'downsample', true, 'combine_layers', true,...
    'inflate', true, 'extract_subcortical_surfs',true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

hemispheres={'lh','rh'};
surfaces={'pial','white','sphere'};
    
% Read RAS offset from freesurfer volume
[status, out]=system(['mri_info --cras ' fullfile(params.subjects_dir,...
    subj_id, 'mri', 'orig.mgz')]);
cols=strsplit(out,' ');
ras_offset=[str2num(cols{1}) str2num(cols{2}) str2num(cols{3})]

% Convert freesurfer surface files to gifti
for s_idx=1:length(surfaces)
    for h_idx=1:length(hemispheres)    
        orig_name=fullfile(params.subjects_dir, subj_id,...
            'surf', sprintf('%s.%s', hemispheres{h_idx}, surfaces{s_idx}));
        new_name=fullfile(params.subjects_dir, subj_id,...
            'surf', sprintf('%s.%s.gii', hemispheres{h_idx}, surfaces{s_idx}));
        system(sprintf('mris_convert %s %s', orig_name, new_name));
   
        % Read in each hemisphere's gifti file and adjust for RAS offset
        g=gifti(new_name);
        % Set transformation matrix to identiy
        g.mat=eye(4);
        g=set_mat(g,'NIFTI_XFORM_UNKNOWN','NIFTI_XFORM_TALAIRACH');
        % Apply RAS offset
        g.vertices=g.vertices+repmat(ras_offset,size(g.vertices,1),1);
        save(g, new_name);
    end
    
    % combine hemispheres
    if params.combine_hemispheres
        lh=fullfile(params.subjects_dir, subj_id, 'surf', sprintf('lh.%s.gii', surfaces{s_idx}));
        rh=fullfile(params.subjects_dir, subj_id, 'surf', sprintf('rh.%s.gii', surfaces{s_idx}));
        combined=fullfile(params.subjects_dir, subj_id, 'surf', sprintf('%s.gii', surfaces{s_idx}));
        combine_surfaces(lh, rh, combined);
    end
end

bem_surfaces={'inner_skull.surf','outer_skull.surf','outer_skin.surf'};

for s_idx=1:length(bem_surfaces)
    orig_name=fullfile(params.subjects_dir, subj_id,...
        'bem', bem_surfaces{s_idx});
    if exist(orig_name,'file')==2
        new_name=fullfile(params.subjects_dir, subj_id,...
            'bem', sprintf('%s.gii', bem_surfaces{s_idx}));
        system(sprintf('mris_convert %s %s', orig_name, new_name));

        % Read in each hemisphere's gifti file and adjust for RAS offset
        g=gifti(new_name);
        % Set transformation matrix to identiy
        g.mat=eye(4);
        g=set_mat(g,'NIFTI_XFORM_UNKNOWN','NIFTI_XFORM_TALAIRACH');
        % Apply RAS offset
        g.vertices=g.vertices+repmat(ras_offset,size(g.vertices,1),1);
        save(g, new_name);    
    end
end

% downsample
if params.downsample
    pial_surf=fullfile(params.subjects_dir, subj_id, 'surf', 'pial.gii');
    ds_pial_surf=fullfile(params.subjects_dir, subj_id, 'surf', 'pial.ds.gii');
    white_surf=fullfile(params.subjects_dir, subj_id, 'surf', 'white.gii');
    ds_white_surf=fullfile(params.subjects_dir, subj_id, 'surf', 'white.ds.gii');
    decimate_two_surfaces(white_surf, pial_surf, ds_white_surf, ds_pial_surf, 0.1);
end

if params.extract_subcortical_surfs
    subsurfs = struct();
    subsurfs(1).code=17;
    subsurfs(1).name='lh.hipp';
    subsurfs(2).code=53;
    subsurfs(2).name='rh.hipp';
    subsurfs(3).code=10;
    subsurfs(3).name='lh.thal';
    subsurfs(4).code=49;
    subsurfs(4).name='rh.thal';
    subsurfs(5).code=16;
    subsurfs(5).name='brainstem';
    
    % Convert freesurfer aseg file to nifti
    system(['mri_convert ' fullfile(params.subjects_dir, subj_id, 'mri', 'aseg.hires.mgz') ' ' fullfile(params.subjects_dir, subj_id, 'mri', 'aseg.hires.nii')]);

    for surf_idx=1:length(subsurfs)
        out_mesh_file=sprintf('%s.gii', subsurfs(surf_idx).name);
        extract_subcortical_surface(subj_id, subsurfs(surf_idx).code,...
            out_mesh_file, 'subjects_dir', params.subjects_dir);

        % downsample
        ds_mesh_file=sprintf('%s.ds.gii', subsurfs(surf_idx).name);
        decimate_surface(fullfile(params.subjects_dir, subj_id, 'surf', out_mesh_file),...
            fullfile(params.subjects_dir, subj_id, 'surf', ds_mesh_file),0.1);
    end
end

% combine layers
if params.combine_layers
    white=fullfile(params.subjects_dir, subj_id, 'surf', 'white.ds.gii');
    pial=fullfile(params.subjects_dir, subj_id, 'surf', 'pial.ds.gii');
    combined=fullfile(params.subjects_dir, subj_id, 'surf', 'white.ds-pial.ds.gii');
    combine_surfaces(white, pial, combined);
end

% inflate
if params.inflate
    surfaces={'pial','white'};
    oldenv=getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH','');
    for s_idx=1:length(surfaces)
        surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.ds.gii',surfaces{s_idx}));
        inflated_surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.ds.inflated.gii',surfaces{s_idx}));
        very_inflated_surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.ds.veryinflated.gii',surfaces{s_idx}));
        system(['wb_command -surface-generate-inflated ' surface_fname ' ' inflated_surface_fname ' ' very_inflated_surface_fname]);
        
        surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.gii',surfaces{s_idx}));
        inflated_surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.inflated.gii',surfaces{s_idx}));
        very_inflated_surface_fname=fullfile(params.subjects_dir, subj_id, 'surf',...
            sprintf('%s.veryinflated.gii',surfaces{s_idx}));
        system(['wb_command -surface-generate-inflated ' surface_fname ' ' inflated_surface_fname ' ' very_inflated_surface_fname]);
    end
    setenv('LD_LIBRARY_PATH',oldenv);
end

