function pial_white_map=map_pial_to_white(white, pial)
% MAP_PIAL_TO_WHITE  Create a mapping from pial to white matter surface vertices based on vertex normal vectors
%
% Use as
%   map_pial_to_white(white, pial)
% where the first argument is the white matter surface gifti object and the second is the pial surface gifti
% object

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
    else
        % If no intersections, get white matter vertex closest to pial vertex
        closest_white=dsearchn(white.vertices,orig);
    end
    pial_white_map(i)=closest_white;    
        
end
