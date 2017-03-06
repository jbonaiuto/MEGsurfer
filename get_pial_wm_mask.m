function [pial_mask,wm_mask,mask]=get_pial_wm_mask(pial_metric, wm_metric,...
    threshold, pial_white_map)
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

% Combined mask
mask=[];
% Pial mask
pial_mask=[];
% White matter mask
wm_mask=[];

% If threshold is not 0
if abs(threshold)>0
    % Find vertices greater than threshold
    if threshold>0
        % Create pial and white masks
        pial_mask=find(pial_metric>=threshold);
        wm_mask=find(wm_metric>=threshold);
    % Find vertices less than threshold
    else
        % Create pial and white masks
        pial_mask=find(pial_metric<threshold);
        wm_mask=find(wm_metric<threshold);        
    end
    % Combined mask is the union of the two
    mask=union(pial_mask,wm_mask(pial_white_map));
end