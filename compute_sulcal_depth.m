function [depth,HS]=compute_sulcal_depth(S)
% COMPUTE_SULCAL_DEPTH  Compute sulcal depth for each vertex in a surface
%
% Use as
%   [depth,HS]=compute_sulcal_depth(S)
% where the first argument is the surface gifti object. Returns
% a vector of depth values with an element for each vertex, and the hull
% surface used to compute depth.
% Requires cat toolbox for spm (http://www.neuro.uni-jena.de/cat/)

% Compute hull surface
HS=cat_surf_fun('hull',S);
% Compute mapping from surface to hull
mapping=dsearchn(HS.vertices,S.vertices);
% Compute distance for each vertex
depth=sqrt(sum((S.vertices-HS.vertices(mapping,:)).^2,2));