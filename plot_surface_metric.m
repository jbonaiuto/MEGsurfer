function [ax, metric_data]=plot_surface_metric( surface, metric_data, varargin)
% PLOT_SURFACE_METRIC  Plots data overlayed on a surface mesh
%
% Use as
%   ax=plot_surface_metric(g, metric_data)
% where the first argument is the surface gifti object, the second is a 
% vector of data to overlay (with one element per surface vertex. Returns
% handle to the axis plotted on and metric data after masking and
% thresholding
%
%   plot_surface_metric(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * output_file - '' (default) or string - name of the file to write the figure
%        to. Only writes if specified
%    * output_format -'png' (default) or string - file format to write to
%    * ax - 0 (default) or handle - handle to the axis to plot on. Creates a new
%        figure if equal to 0
%    * limits - [] (default) or two element vector - Limits for overlay
%    color scale
%    * surface_alpha - 0.1 (default) or double - surface alpha value. Passed as FaceAlpha
%        to patch
%    * mask - [] (default) or array - list of vertex indices to include in
%        overlay. Other vertices are plotted as gray
%    * title - '' (default) or string - plot title
%    * threshold - 0 (default) or float - only plot vertices with values
%        greater than (if positive) or less than (if negative) threshold
%    * clip_vals - true (default) or false - whether or not to clip values
%        in color scale by 2 and 98 percentile of overlay data
%    * clabel - '' (default) or string - label for color map
%    * custom_cm - true (default) or false - whether or not to use a custom color map
%    * specular_strength - 0.0 (default) - light specular strength
%    * ambient_strength - 0.75 (default) - light ambient strength
%    * diffuse_strength - 1.0 (default) - light diffuse strength
%    * face_lighting - phong (default) - face lighting

% Parse inputs
defaults = struct('output_file', '', 'output_format', 'png', 'ax', 0,...
    'limits', [], 'mask', [], 'title', '', 'threshold', 0, ...
    'clip_vals', true, 'clabel' ,'', 'custom_cm', true, ...
    'specular_strength', 0.0, 'ambient_strength', 0.75,...
    'diffuse_strength', 1.0, 'face_lighting', 'phong');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Create figure and axis if not specified
if params.ax==0
    fig=figure('Renderer','OpenGL', 'Color',[1 1 1]);
    params.ax = axes('Parent', fig);
end
set(params.ax,'Visible','off');
set(params.ax,'Projection','perspective');
set(params.ax,'PlotBoxAspectRatio',[1.268 1 1.129]);
set(params.ax,'DataAspectRatio',[1 1 1]);


% Plot surface
hp = patch('vertices', surface.vertices, 'faces', surface.faces,...
    'EdgeColor', 'none', 'Parent', params.ax, 'FaceColor','interp',...
    'linestyle','none','AmbientStrength',params.ambient_strength,...,
    'DiffuseStrength', params.diffuse_strength,...
    'SpecularStrength',params.specular_strength);
if length(params.face_lighting)
    set(hp, 'FaceLighting', params.face_lighting);
end
 
% Create light
light('Parent',params.ax,'Style','local','Position',[749 868.1 1263]);
 
% Create light
light('Parent',params.ax,'Style','local','Position',[-1611 -265.3 288.5]);
 
% Create light
light('Visible','off','Parent',params.ax,'Style','local','Position',[776.2 899.5 1309]);

% Create color maps
% Control points for positive color map
pos_control_pts=zeros(3,5);
pos_control_pts(1,:)=[340 225 90 10 10];
pos_control_pts(2,:)=[340 340 340 270 10];
pos_control_pts(3,:)=[340 340 340 340 340];

% Control points for negative color map
neg_control_pts=zeros(3,5);
neg_control_pts(1,:)=[340 264 264 180 340];
neg_control_pts(2,:)=[10 92 92 264 340];
neg_control_pts(3,:)=[10 264 264 92 264];

% If metric value limits are not specified
if length(params.limits)==0

    % If clip values by percentile
    if params.clip_vals
        % Split metric data into positive and negative values
        pos_vals=find(metric_data>0);
        neg_vals=find(metric_data<0);

        % Compute percentiles
        if length(neg_vals) && length(pos_vals)
            pos_percentiles=tiedrank(metric_data(pos_vals))/length(pos_vals);
            neg_percentiles=tiedrank(metric_data(neg_vals))/length(neg_vals);

            % Min and max are 2% and 98% of neg and pos values
            min_clipped_val=min(metric_data(neg_vals(find(neg_percentiles>=.02))));
            max_clipped_val=max(metric_data(pos_vals(find(pos_percentiles<=.98))));
            params.limits=[min_clipped_val max_clipped_val];            
        else
            percentiles=tiedrank(metric_data)/length(metric_data);

            % Min and max are 2% and 98% of neg and pos values
            min_clipped_val=min(metric_data(find(percentiles>=.02)));
            max_clipped_val=max(metric_data(find(percentiles<=.98)));
            params.limits=[min_clipped_val max_clipped_val];
        end
    else
        params.limits=[min(metric_data) max(metric_data)];
    end

end

% Number of colors overall
num_colors=255;
% Compute ratio of positive to negative values
if params.limits(1)<0
    if params.limits(2)<0
        ratio=1;
    else
        % Compute ratio of positive to negative values
        ratio=abs(params.limits(1))/(abs(params.limits(1))+params.limits(2));
    end
else
    ratio=0;
end
% Compute number of postive and negative colots
neglen=round(num_colors*ratio);
poslen=num_colors-neglen;

% Generate positive and negative color maps
C_pos=generate_color_map(poslen, pos_control_pts);
C_neg=generate_color_map(neglen, neg_control_pts);

% Set color map
if params.custom_cm
    cm=colormap([C_neg;C_pos]/255);
else
    cm=colormap();
end

% Show colorbar
set(params.ax, 'Clim', params.limits);
cb=colorbar();
if length(params.clabel)
    ylabel(cb,params.clabel)
end
freezeColors(params.ax);
drawnow();
tick_percentages=(get(cb,'Ticks')-params.limits(1))/(params.limits(2)-params.limits(1));
tick_labels=get(cb,'TickLabels');

% Generate display version of metric data
display_metric_data=metric_data;
% Clip range at the limits
display_metric_data(find(display_metric_data<params.limits(1)))=params.limits(1);
display_metric_data(find(display_metric_data>params.limits(2)))=params.limits(2);

% Calculate color number that each intensity maps to
numcol = size(colormap,1);
cfrac = (display_metric_data - params.limits(1)) ./ (params.limits(2) - params.limits(1));
display_metric_data = 1 + floor(numcol * (cfrac * (1-2*eps)));
% Set mapping to direct - intensity values are index of colors in map
set(hp,'CDataMapping', 'direct');


% Apply mask
if length(params.mask)>0 || abs(params.threshold)>0
    % Add extra color (gray to the end of the map)
    colormap([colormap; .5 .5 .5]);
    % Get new size of color map
    numcol = size(colormap,1);

    % Apply mask
    if length(params.mask)>0
        % Apply mask to metric data to return
        metric_data=metric_data(params.mask);
        % Index of all faces
        all_faces=[1:length(display_metric_data)];
        % Get faces that aren't masked
        unmasked_faces=setdiff(all_faces,params.mask);
        % For plotting set unmasked data to extra color - gray
        display_metric_data(unmasked_faces)=numcol+2;
    end    

    % Apply threshold
    if abs(params.threshold)>0
        if params.threshold<0
            display_metric_data(find(metric_data>params.threshold))=numcol;
            metric_data=metric_data(find(metric_data<=params.threshold));
        else
            display_metric_data(find(metric_data<params.threshold))=numcol;
            metric_data=metric_data(find(metric_data>=params.threshold));
        end
    end
end

% Set face color data
set(hp,'FaceVertexCData', display_metric_data);
new_limits=get(cb,'Limits');
new_ticks=(new_limits(2)-new_limits(1))*tick_percentages+new_limits(1);
set(cb,'Ticks',new_ticks);
set(cb,'TickLabels',tick_labels);
axes(params.ax);
cameramenu;
freezeColors(params.ax);

freezeColors(params.ax);

ax=params.ax;

% Set title
if length(params.title)
    annotation(fig,'textbox',[0.15 0.85 0.2 .15],'String',{params.title},'FitBoxToText','on');
end

% Save plot to file
if length(params.output_file)>0
    if strcmp(params.output_format,'eps')
        figure2eps(fig, params.output_file, 10, '-opengl');
    else
        saveas(fig, params.output_file, params.output_format);
    end
end
