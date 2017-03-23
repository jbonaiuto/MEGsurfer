function [pial_mask,wm_mask,mask]=get_pial_wm_mask(pial_metric, wm_metric,...
    threshold, pial_white_map, varargin)
% GET_PIAL_WM_MASK  Create a mask of pial, white matter and combined
% surface based on threshold
%
% Use as
%   [pial_mask,wm_mask,mask]=get_pial_wm_mask(pial_metric, wm_metric, threshold, pial_white_map)
% where the first argument is the pial surface data and the second is the
% white matter surface data, the third is the threshold, and the fifth is a
% mapping from the pial to white matter surface. It returns a mask for the 
% pial surface (where the data on that surface exceeds the threshold, one
% for the white matter surface(where the data on that surface exceeds the
% threshold), and one fo the combine surface (where vertices on either
% surface exceed the threshold).
%   get_pial_wm_mask(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * thresh_type - 'lower' (default) or 'upper' - whether threshold is an
%    upper or lower limit

% Parse inputs
defaults = struct('thresh_type', 'lower');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Combined mask
mask=[];
% Pial mask
pial_mask=[];
% White matter mask
wm_mask=[];

% If threshold is not []
if length(threshold)>0
    % Mapped white matter data
    mapped_wm_metric=wm_metric(pial_white_map);
    
    % Find vertices greater than threshold
    if strcmp(params.thresh_type,'lower')
        % Create pial and white masks and mapped white mask
        pial_mask=find(pial_metric>threshold);
        wm_mask=find(wm_metric>threshold);
        mapped_wm_mask=find(mapped_wm_metric>threshold);
    % Find vertices less than threshold
    else
        % Create pial and white maskss and mapped white mask
        pial_mask=find(pial_metric<threshold);
        wm_mask=find(wm_metric<threshold);        
        mapped_wm_mask=find(mapped_wm_metric<threshold);        
    end
    % Combined mask is the union of the pial and mapped white mask
    mask=union(pial_mask,mapped_wm_mask);
end