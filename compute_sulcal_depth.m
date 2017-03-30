function [depth,HS]=compute_sulcal_depth(surface_file)
% COMPUTE_SULCAL_DEPTH  Compute sulcal depth for each vertex in a surface
%
% Use as
%   [depth,HS]=compute_sulcal_depth(S)
% where the first argument is the surface gifti object. Returns
% a vector of depth values with an element for each vertex, and the hull
% surface used to compute depth.
% Requires cat toolbox for spm (http://www.neuro.uni-jena.de/cat/)

spm('defaults','eeg');
S=gifti(surface_file);

[path file ext]=fileparts(surface_file);
hull_file=fullfile(path, sprintf('%s_hull.gii',file));
if exist(hull_file,'file')~=2
    % Compute hull surface
    hull=cat_surf_fun('hull',S);
    HS=gifti();
    HS.vertices=hull.vertices;
    HS.faces=hull.faces;
    save(HS, hull_file);
else
    HS=gifti(hull_file);
end
% Compute mapping from surface to hull
mapping=dsearchn(HS.vertices,S.vertices);
% Compute distance for each vertex
depth=sqrt(sum((S.vertices-HS.vertices(mapping,:)).^2,2));
