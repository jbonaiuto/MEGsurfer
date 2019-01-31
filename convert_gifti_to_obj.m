function convert_gifti_to_obj(in_name, out_name)

orig=gifti(in_name);
v=orig.vertices;
f=orig.faces; 
vertface2obj(v,f,out_name);
