function decimate_two_surfaces(in_surf1, in_surf2, out_surf1, out_surf2, ratio)

surf1=gifti(in_surf1);
surf2=gifti(in_surf2);

FV=[];
FV.faces=surf1.faces;
FV.vertices=surf1.vertices;
FV_reduced=reducepatch(FV, ratio);

orig_vert_idx=dsearchn(surf1.vertices,FV_reduced.vertices);

write_surf_gifti(out_surf1, FV_reduced.vertices, FV_reduced.faces);
write_surf_gifti(out_surf2, surf2.vertices(orig_vert_idx,:), FV_reduced.faces);

end

