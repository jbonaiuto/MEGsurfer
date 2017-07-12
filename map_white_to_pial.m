function white_pial_map=map_white_to_pial(white_name, pial_name, varargin)
% MAP_WHITE_TO_PIAL  Create a mapping from white matter to pial surface vertices
%
% Use as
%   white_pial_map=map_white_to_pial(white, pial)
% where the first argument is the pial surface gifti object and the second is the white matter surface gifti
% object. It returns a vector with length equal to the number of vertices in the white matter surface. Each vector value
% is the index of the corresponding vertex of the pial surface.
%
%   map_white_to_pial(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * mapType - 'nearest' (default) or 'normal' - how to map to the
%        nearest vertex in the pial surface:
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
map_file=fullfile(meshpath, sprintf('map_white_pial_%s_%s-%s.mat', params.mapType, pial_file, white_file));
if exist(map_file,'file')~=2 || params.recompute
    disp('Recomputing mapping');
    white=gifti(white_name);
    pial=gifti(pial_name);
    n_vertices=size(pial.vertices,1);
    white_pial_map=[1:n_vertices];
    switch params.mapType
        case 'nearest' % Maps to the nearest vertex on the pial surface
            % Get index of nearest pial vertex for each white vertex
            white_pial_map=dsearchn(pial.vertices,white.vertices);
        case 'normal' % Maps to the nearest intersection in the normal vector direction on the other surface
            % Compute normal vectors for each white surface vertex - pointing outward
            % TODO: how to make sure that they're pointing outward?
            normal=spm_mesh_normals(struct('vertices',white.vertices,'faces',white.faces),true);
        
            % Get x, y, z for each pial surface vertex
            vert1 = pial.vertices(pial.faces(:,1),:);
            vert2 = pial.vertices(pial.faces(:,2),:);
            vert3 = pial.vertices(pial.faces(:,3),:);
        
            white_pial_map=dsearchn(pial.vertices,white.vertices);

            % Iterate over each white vertex
            for i=1:n_vertices

                % Get intersection of ray originating from white vertex and pointing in normal direction with the white matter
                % surface
                orig=white.vertices(i,:);

                [intersected,t,u,v,xcoor] = TriangleRayIntersection(orig, normal(i,:), vert1, vert2, vert3, 'border', 'inclusive' ,'lineType', 'line');
                % If there are any intersections
                intersected_idx=find(intersected>0);
                if length(intersected_idx)
                    % Get closest intersection point on pial mesh
                    nearest_intersection_idx=dsearchn(xcoor(intersected_idx,:),orig);
                    % Get closest pial vertex to intersection
                    closest_pial=dsearchn(pial.vertices,xcoor(intersected_idx(nearest_intersection_idx),:));
                    white_pial_map(i)=closest_pial;
                end            
            end
        case 'link' % Maps to the corresponding vertex according to freesurfer (for downsampled surfaces where correspondence is destroyed)
            orig_pial=gifti(params.origPial);
            orig_white=gifti(params.origWhite);
            % Map from downsampled white to original white surface
            white_ds_white_map=dsearchn(orig_white.vertices,white.vertices);
            % Map from original pial surface to downsampled pial surface
            ds_pial_pial_map=dsearchn(pial.vertices,orig_pial.vertices);
            % Apply original pial -> downsampled pial mapping to downsampled white -> original white mapping
            white_pial_map=ds_pial_pial_map(white_ds_white_map);
    end
    save(map_file,'white_pial_map');
else
    a=load(map_file);
    white_pial_map=a.white_pial_map;
end
