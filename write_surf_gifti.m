function write_surf_gifti(filename, vertices, faces)

s=gifti();
% Set transformation matrix to identiy
s.mat=eye(4);
s=set_mat(s,'NIFTI_XFORM_UNKNOWN','NIFTI_XFORM_TALAIRACH');
s.vertices=vertices;
s.faces=faces;
save(s, filename);

