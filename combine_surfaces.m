function combine_surfaces(surfaces, out_fname, varargin)
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
face_offset=0;
mesh_one=gifti(surfaces{1});
% Combine vertices and faces
combined_vertices=[mesh_one.vertices];
combined_faces=[mesh_one.faces+face_offset];
combined_normals=[];
if isfield(mesh_one,'normals')
    combined_normals=[mesh_one.normals];
end
% mesh two face number will start with max mesh one face number+1
face_offset=face_offset+max(mesh_one.faces(:));

for s=2:length(surfaces)
    mesh_two=gifti(surfaces{s});

    % Combine vertices and faces
    combined_vertices=[combined_vertices; mesh_two.vertices];
    combined_faces=[combined_faces; mesh_two.faces+face_offset];
    if length(combined_normals)>0 && isfield(mesh_two,'normals')
        combined_normals=[combined_normals; mesh_two.normals];
    end
    face_offset=face_offset+max(mesh_two.faces(:));
end

% Create and save combined gifti
write_surf_gifti(out_fname, combined_vertices, combined_faces, 'normals', combined_normals);
