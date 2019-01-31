function convert_obj_to_gifti(in_name, out_name)

orig=readObj(in_name);
write_surf_gifti(out_name,orig.v,orig.f.v);

