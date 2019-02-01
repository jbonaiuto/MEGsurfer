function decimate_two_surfaces(in_surf1, in_surf2, out_surf1, out_surf2, ratio)

surf1=gifti(in_surf1);
surf2=gifti(in_surf2);

FV_reduced = reducepatch(patch('Faces', surf1.faces, 'Vertices', surf1.vertices), ratio);

[orig_vert_idx, decim_orig_dist] = knnsearch(surf1.vertices, FV_reduced.vertices);
disp(['',num2str(mean(~decim_orig_dist)*100),'% of the vertices in the decimated first surface belong to the initial first surface vertices.'])

write_surf_gifti(out_surf1, FV_reduced.vertices, FV_reduced.faces);
write_surf_gifti(out_surf2, surf2.vertices(orig_vert_idx,:), FV_reduced.faces);

end

