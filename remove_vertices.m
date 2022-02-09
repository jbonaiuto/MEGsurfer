function surf=remove_vertices(surf, vertices_to_remove)

vertices_to_keep=setdiff([1:size(surf.vertices,1)], vertices_to_remove);

% Create new array of vertices
new_vertices=surf.vertices(vertices_to_keep,:);

% Find which faces to keep - ones that point to kept vertices
face_x=find(ismember(surf.faces(:,1),vertices_to_keep));
face_y=find(ismember(surf.faces(:,2),vertices_to_keep));
face_z=find(ismember(surf.faces(:,3),vertices_to_keep));
faces_to_keep=intersect(intersect(face_x,face_y),face_z);

% Re-index faces
x_faces=reshape(surf.faces(faces_to_keep,:),1,length(faces_to_keep)*3);
[y,idxs]=ismember(x_faces,vertices_to_keep);
new_faces=reshape(idxs,length(faces_to_keep),3);

% Create new surface
surf.vertices=new_vertices;
surf.faces=int32(new_faces);
if isfield(surf,'normals')
    surf.normals=surf.normals(vertices_to_keep,:);
end