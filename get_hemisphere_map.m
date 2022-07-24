function hemisphere_map=get_hemisphere_map(surface_name, orig_surface_name, varargin)
% 1=left, 2=right

[meshpath mesh_file ext]=fileparts(surface_name);

surface=gifti(surface_name);
n_vertices=size(surface.vertices,1);
hemisphere_map=ones(1,n_vertices);
orig_surface=gifti(orig_surface_name);

[meshpath orig_mesh_file ext]=fileparts(orig_surface_name);
lh_orig_surface=gifti(fullfile(meshpath,sprintf('lh.%s.gii',orig_mesh_file)));
n_lh_vertices=size(lh_orig_surface.vertices,1);

% Map from downsampled pial to original pial surface
surf_ds_surf_map=knnsearch(orig_surface.vertices,surface.vertices);
hemisphere_map(find(surf_ds_surf_map>n_lh_vertices))=2;

