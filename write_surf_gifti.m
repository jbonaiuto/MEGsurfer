function write_surf_gifti(filename, vertices, faces, varargin)

% Parse inputs
defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end


s=gifti();
% Set transformation matrix to identiy
s.mat=eye(4);
s=set_mat(s,'NIFTI_XFORM_UNKNOWN','NIFTI_XFORM_TALAIRACH');
s.vertices=vertices;
s.faces=faces;
save(s, filename);

