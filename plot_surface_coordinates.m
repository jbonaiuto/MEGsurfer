function ax=plot_surface_coordinates(g, coords, coord_color, varargin)
% PLOT_SURFACE_COORDINATES  Plots a surface mesh along with coordinates
%
% Use as
%   ax=plot_surface_coordinates(g, coords, coord_color)
% where the first argument is the surface gifti object, the second is a Nx3 array of
% coordinates, and the third is the color of the coordinates (passed as FaceColor
% to surface). Returns handle to the axis plotted on
%
%   plot_surface_coordinates(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * output_file - '' (default) or string - name of the file to write the figure
%        to. Only writes if specified
%    * output_format -'png' (default) or string - file format to write to
%    * ax - 0 (default) or handle - handle to the axis to plot on. Creates a new
%        figure if equal to 0
%    * title - '' (default) or string - Plot title
%    * surface_alpha - 0.1 (default) or double - surface alpha value. Passed as FaceAlpha
%        to patch
%    * coord_radius - 2 (default) or double - radius of plotted coordinates

% Parse inputs
defaults = struct('output_file', '', 'output_format', 'png', 'ax', 0, 'title', '', 'surface_alpha', 0.1, 'coord_radius', 2);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Plot surface
ax=plot_surface(g, 'ax', params.ax, 'title', params.title, 'surface_alpha', params.surface_alpha);
fig=get(ax,'Parent');

% Plot coordinates
[x,y,z]=sphere();
for s=1:size(coords,1)
    hp=surface(x.*params.coord_radius+double(coords(s,1)),y.*params.coord_radius+double(coords(s,2)),z.*params.coord_radius+double(coords(s,3)),...
       'FaceColor',coord_color,'EdgeColor','none','linestyle','none','FaceLighting','phong');
end

% Save plot to file
if length(params.output_file)>0
    if params.output_format=='eps'
        figure2eps(fig, params.output_file, 10, '-opengl');
    else
        saveas(fig, params.output_file, params.output_format);
    end
end

