function decimate_surface(in_surf1, out_surf1, ratio)

surf1=gifti(in_surf1);

FV_reduced = reducepatch(patch('Faces', surf1.faces, 'Vertices', surf1.vertices), ratio);

write_surf_gifti(out_surf1, FV_reduced.vertices, FV_reduced.faces);

end

