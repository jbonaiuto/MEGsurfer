function pial_white_map=map_pial_to_white(white, pial, varargin)
% MAP_PIAL_TO_WHITE  Create a mapping from pial to white matter surface vertices based on vertex normal vectors
%
% Use as
%   map_pial_to_white(white, pial)
% where the first argument is the white matter surface gifti object and the second is the pial surface gifti
% object
%
%   map_pial_to_white(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * mapType - 'nearest' (default) or 'normal' - how to map to the
%        nearest vertex in the white matter surface:
%        - 'nearest' means map to the closest vertex in any direction
%        - 'normal' means map to the closest vertex in the direction of
%             the normal vector

% Parse inputs
defaults = struct('mapType','nearest');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

switch params.mapType
    case 'nearest'
        pial_white_map=dsearchn(white.vertices,pial.vertices);
    case 'normal'
        % Compute normal vectors for each pial surface vertex - pointing inward
        normal_vecs=spm_mesh_normals(struct('faces',pial.faces,'vertices',pial.vertices),true);

        % Get x, y, z for each white matter surface vertex
        vert1 = white.vertices(white.faces(:,1),:);
        vert2 = white.vertices(white.faces(:,2),:);
        vert3 = white.vertices(white.faces(:,3),:);

        n_vertices=size(pial.vertices,1);
        pial_white_map=zeros(n_vertices,1);

        % Iterate over each pial vertex
        for i=1:n_vertices

            % Get intersection of ray originating from pial vertex and pointing in normal direction with the white matter
            % surface
            orig=pial.vertices(i,:);
            [intersected,t,u,v,xcoor] = TriangleRayIntersection(orig, normal_vecs(i,:), vert1, vert2, vert3, 'border', 'inclusive');

            % If there are any intersections
            intersected_idx=find(intersected>0);
            if length(intersected_idx)
                % Get coordinates and distances of intersections
                intersections=xcoor(intersected_idx,:);
                dists=t(intersected_idx);
                % Get closest interection
                [min_dist,closest_intersection]=min(dists);
                % Get white matter surface vertex closest to closest intersection
                closest_white=dsearchn(white.vertices,intersections(closest_intersection,:));        
            end
            if length(intersected_idx)==0 || sqrt(sum((closest_white-orig).^2))>50
                % If no intersections, get white matter vertex closest to pial vertex
                closest_white=dsearchn(white.vertices,orig);
            end
            pial_white_map(i)=closest_white;    
        
        end
end
