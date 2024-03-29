function create_multilayer_surface(subj_id, mri_path, surf_path)

% Parse inputs
defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

hemispheres={'lh','rh'};
layers=[0.1:0.1:0.9];
    
  
for l=1:length(layers)
    layer=layers(l);
    for h=1:length(hemispheres)
        hemi=hemispheres{h};
        wm_file=fullfile(surf_path, sprintf('%s.white',hemi));
        out_file=fullfile(surf_path, sprintf('%s.%.1f',hemi,layer));
        [status, out]=system(['mris_expand -thickness ' wm_file ' ' sprintf('%.1f',layer) ' ' out_file]);        
    end
end

% Read RAS offset from freesurfer volume
[status, out]=system(['mri_info --cras ' fullfile(mri_path, 'orig.mgz')]);
cols=strsplit(out,' ');
ras_offset=[str2num(cols{1}) str2num(cols{2}) str2num(cols{3})]

% Convert freesurfer surface files to gifti
for l=1:length(layers)
    layer=layers(l);
    for h_idx=1:length(hemispheres)    
        hemi=hemispheres{h_idx};
        orig_name=fullfile(surf_path, sprintf('%s.%.1f', hemi, layer));
        new_name=fullfile(surf_path, sprintf('%s.%.1f.gii', hemi, layer));
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
    lh=fullfile(surf_path sprintf('lh.%.1f.gii', layer));
    rh=fullfile(surf_path, sprintf('rh.%.1f.gii', layer));
    combined=fullfile(surf_path, sprintf('%.1f.gii', layer));
    combine_surfaces({lh, rh}, combined);
end

% downsample
in_surfs={fullfile(surf_path, 'white.gii')};
out_surfs={fullfile(surf_path, 'white.ds.gii')};
for l=1:length(layers)
    layer=layers(l);
    in_surfs{end+1}=fullfile(surf_path, sprintf('%.1f.gii', layer));
    out_surfs{end+1}=fullfile(surf_path, sprintf('%.1f.ds.gii', layer));
end
in_surfs{end+1}=fullfile(surf_path, 'pial.gii');
out_surfs{end+1}=fullfile(surf_path, 'pial.ds.gii');
decimate_multiple_surfaces(in_surfs, out_surfs, 0.1);

combined=fullfile(surf_path, 'multilayer.ds.gii');
% reverse order so surface order matches electrode order in laminar recordings
out_surfs(end:-1:1) = out_surfs(:);
combine_surfaces(out_surfs, combined);

