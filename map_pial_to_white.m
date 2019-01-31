function pial_white_map=map_pial_to_white(white_name, pial_name, varargin)
% MAP_PIAL_TO_WHITE  Create a mapping from pial to white matter surface vertices
%
% Use as
%   pial_white_map=map_pial_to_white(white, pial)
% where the first argument is the white matter surface gifti object and the second is the pial surface gifti
% object. It returns a vector with length equal to the number of vertices in the pial surface. Each vector value
% is the index of the corresponding vertex of the white matter surface.
%
%   map_pial_to_white(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * mapType - 'nearest' (default) or 'normal' - how to map to the
%        nearest vertex in the white matter surface:
%        - 'nearest' means map to the closest vertex in any direction
%        - 'normal' means map to the closest vertex in the direction of
%             the normal vector
%        - 'link' means using the correspondence between non-downsampled surfaces
%    * origPial - '' (default) filename for original pial surface, only used when mapType=link
%    * origWhite - '' (default) filename for original white surface, only used when mapType=link
%    * recompute - false (default) whether or not to recompute the mapping even if it already exists

% Parse inputs
defaults = struct('mapType','link','origPial','','origWhite','','recompute',false);  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

[meshpath white_file ext]=fileparts(white_name);
[meshpath pial_file ext]=fileparts(pial_name);
map_file=fullfile(meshpath, sprintf('map_pial_white_%s_%s-%s.mat', params.mapType, pial_file, white_file));
if exist(map_file,'file')~=2 || params.recompute
    disp('Recomputing mapping');
    white=gifti(white_name);
    pial=gifti(pial_name);
    n_vertices=size(pial.vertices,1);
    pial_white_map=[1:n_vertices];
    switch params.mapType
        case 'nearest' % Maps to the nearest vertex on the white surface
            % Get index of nearest white vertex for each pial vertex
            pial_white_map=knnsearch(white.vertices,pial.vertices); 
        case 'normal' % Maps to the nearest intersection in the normal vector direction on the other surface
            % Compute normal vectors for each pial surface vertex - pointing outward
            % TODO: how to make sure that they're pointing outward?
            normal=spm_mesh_normals(struct('vertices',pial.vertices,'faces',pial.faces),true);
        
            % Get x, y, z for each white matter surface vertex
            vert1 = white.vertices(white.faces(:,1),:);
            vert2 = white.vertices(white.faces(:,2),:);
            vert3 = white.vertices(white.faces(:,3),:);
        
            pial_white_map=knnsearch(white.vertices,pial.vertices);

            % Iterate over each pial vertex
            for i=1:n_vertices

                % Get intersection of ray originating from pial vertex and pointing in normal direction with the pial
                % surface
                orig=pial.vertices(i,:);

                [intersected,t,u,v,xcoor] = TriangleRayIntersection(orig, normal(i,:), vert1, vert2, vert3, 'border', 'inclusive' ,'lineType', 'line');
                % If there are any intersections
                intersected_idx=find(intersected>0);
                if length(intersected_idx)
                    % Get closest intersection point on white mesh
                    nearest_intersection_idx=knnsearch(xcoor(intersected_idx,:),orig);
                    % Get closest white vertex to intersection
                    closest_white=knnsearch(white.vertices,xcoor(intersected_idx(nearest_intersection_idx),:));
                    pial_white_map(i)=closest_white;
                end            
            end
        case 'link' % Maps to the corresponding vertex according to freesurfer (for downsampled surfaces where correspondence is destroyed)
            orig_pial=gifti(params.origPial);
            orig_white=gifti(params.origWhite);
            % Map from downsampled pial to original pial surface
            pial_ds_pial_map=knnsearch(orig_pial.vertices,pial.vertices);
            % Map from original white surface to downsampled white surface
            ds_white_white_map=knnsearch(white.vertices,orig_white.vertices);
            % Apply original white -> downsampled white mapping to downsampled pial -> original pial mapping
            pial_white_map=ds_white_white_map(pial_ds_pial_map);
    end
    save(map_file,'pial_white_map');
else
    a=load(map_file);
    pial_white_map=a.pial_white_map;
end
