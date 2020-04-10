function combine_surfaces(surface_one_fname, surface_two_fname, out_fname, varargin)
% function combine_surfaces(surface_one_fname, surface_two_fname, out_fname, varargin)
% Combines gifti surfaces into single gifti
% INPUT: 
%   surface_one_fname: surface filename
%   surface_two_fname: surface fielname
%   out_fname: new surface filename
% ---------------------------
% v1.0 James Bonaiuto (james.bonaiuto@isc.cnrs.fr)
% 

% Parse inputs
defaults = struct();  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Read in each surface's gifti file
mesh_one=gifti(surface_one_fname);
mesh_two=gifti(surface_two_fname);

% mesh two face number will start with max mesh one face number+1
face_offset=max(mesh_one.faces(:));

% Combine vertices and faces
combined_vertices=[mesh_one.vertices; mesh_two.vertices];
combined_faces=[mesh_one.faces; mesh_two.faces+face_offset];
combined_normals=[];
if isfield(mesh_one,'normals') && isfield(mesh_two,'normals')
    combined_normals=[mesh_one.normals; mesh_two.normals];
end

% Create and save combined gifti
write_surf_gifti(out_fname, combined_vertices, combined_faces, 'normals', combined_normals);
