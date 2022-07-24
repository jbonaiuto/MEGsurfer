function decimate_multiple_surfaces(in_surfs, out_surfs, ratio)

primary_surf=gifti(in_surfs{1});
FV_reduced = reducepatch(patch('Faces', primary_surf.faces, 'Vertices', primary_surf.vertices), ratio);
[orig_vert_idx, decim_orig_dist] = knnsearch(primary_surf.vertices, FV_reduced.vertices);
disp(['',num2str(mean(~decim_orig_dist)*100),'% of the vertices in the decimated first surface belong to the initial first surface vertices.'])

write_surf_gifti(out_surfs{1}, FV_reduced.vertices, FV_reduced.faces);

for s=2:length(in_surfs)
    surf2=gifti(in_surfs{s});
    write_surf_gifti(out_surfs{s}, surf2.vertices(orig_vert_idx,:), FV_reduced.faces);
end

