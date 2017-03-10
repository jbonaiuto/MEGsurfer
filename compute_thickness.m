function [pial_thickness, wm_thickness]=compute_thickness(ds_pial_file, ...
    ds_white_file, orig_pial_file, orig_white_file, varargin)
% COMPUTE_THICKNESS  Compute cortical thickness for each vertex in a pial
% surface
%
% Use as
%   [pial_thickness wm_thickness]=compute_thickness(ds_pial_file, ...
%       ds_white_file, orig_pial_file, orig_white_file, varargin)
% where the first argument is the downsampled pial surface filename, the
% second is the downsampled white matter surface filename, the third is the
% original pial surface filename and the fourth is the original white
% matter surface filename. Returns a vector of thickness values with an 
% element for each downsampled pial vertex, and a vector of thickness
% values with an element for each downsampled white matter vertex
%
%   compute_thickness(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * smooth - true (default) or false - whether or not to smooth the
%    thickness values

% Parse inputs
defaults = struct('smooth', true);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Compute mapping between surfaces
pial_white_map=map_pial_to_white(ds_white_file, ds_pial_file, ...
        'mapType', 'link', 'origPial', orig_pial_file, ...
        'origWhite', orig_white_file);
white_pial_map=map_white_to_pial(ds_white_file, ds_pial_file, ...
    'mapType', 'link', 'origPial', orig_pial_file, ...
    'origWhite', orig_white_file);
    
% Load suface
pial=gifti(ds_pial_file);
wm=gifti(ds_white_file);

% Compute thickness from each surface
pial_thickness=sqrt(sum((pial.vertices-wm.vertices(pial_white_map,:)).^2,2));
wm_thickness=sqrt(sum((pial.vertices(white_pial_map,:)-wm.vertices).^2,2));

% Smooth thickness
if params.smooth
    pial_thickness=spm_mesh_smooth(pial, pial_thickness, 5);
    wm_thickness=spm_mesh_smooth(wm, wm_thickness, 5);
end