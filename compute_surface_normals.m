function norm_vecs=compute_surface_normals(subjects_dir, subj_id, surface, method)

surf_dir=fullfile(subjects_dir, subj_id, 'surf');

if strcmp(surface,'white-pial')
    ds_surface=gifti(fullfile(surf_dir, 'white.ds-pial.ds.gii'));
else
    ds_surface=gifti(fullfile(surf_dir, sprintf('%s.ds.gii', surface)));
end
norm_vecs = spm_mesh_normals(struct('faces',ds_surface.faces,...
    'vertices',ds_surface.vertices),true);

switch method
    case 'cps'
        if strcmp(surface,'white-pial')
            ds_white=gifti(fullfile(surf_dir, 'white.ds.gii'));
            orig_white=gifti(fullfile(surf_dir, 'white.gii'));
            ds_white_idx=knnsearch(ds_white.vertices,orig_white.vertices);
            orig_white_norm=spm_mesh_normals(struct('faces',orig_white.faces,...
                'vertices',orig_white.vertices),true);
            white_norm_vecs=zeros(size(ds_white.vertices,1),3);
            for v_idx=1:size(white_norm_vecs,1)
                white_norm_vecs(v_idx,:)=mean(orig_white_norm(ds_white_idx==v_idx,:),1);
            end
            
            ds_pial=gifti(fullfile(surf_dir, 'pial.ds.gii'));
            orig_pial=gifti(fullfile(surf_dir, 'pial.gii'));
            ds_pial_idx=knnsearch(ds_pial.vertices,orig_pial.vertices);
            orig_pial_norm=spm_mesh_normals(struct('faces',orig_pial.faces,...
                'vertices',orig_pial.vertices),true);
            pial_norm_vecs=zeros(size(ds_pial.vertices,1),3);
            for v_idx=1:size(pial_norm_vecs,1)
                pial_norm_vecs(v_idx,:)=mean(orig_pial_norm(ds_pial_idx==v_idx,:),1);
            end
            norm_vecs=[white_norm_vecs; pial_norm_vecs];
        else
            orig_surface=gifti(fullfile(surf_dir, sprintf('%s.gii', surface)));
            ds_idx=knnsearch(ds_surface.vertices,orig_surface.vertices);
            orig_norm=spm_mesh_normals(struct('faces',orig_surface.faces,...
                'vertices',orig_surface.vertices),true);            
            for v_idx=1:size(norm_vecs,1)
                norm_vecs(v_idx,:)=mean(orig_norm(ds_idx==v_idx,:),1);
            end
        end
        
    %% Orig surface normals
    case 'orig_surf_norm'
        if strcmp(surface,'white-pial')
            ds_white=gifti(fullfile(surf_dir, 'white.ds.gii'));
            orig_white=gifti(fullfile(surf_dir, 'white.gii'));
            orig_white_idx=knnsearch(orig_white.vertices,ds_white.vertices);
            white_norm_vecs=spm_mesh_normals(struct('faces',orig_white.faces,...
                'vertices',orig_white.vertices),true);
            white_norm_vecs=white_norm_vecs(orig_white_idx,:);
            
            ds_pial=gifti(fullfile(surf_dir, 'pial.ds.gii'));
            orig_pial=gifti(fullfile(surf_dir, 'pial.gii'));
            orig_pial_idx=knnsearch(orig_pial.vertices,ds_pial.vertices);
            pial_norm_vecs=spm_mesh_normals(struct('faces',orig_pial.faces,...
                'vertices',orig_pial.vertices),true);
            pial_norm_vecs=pial_norm_vecs(orig_pial_idx,:);
            norm_vecs=[white_norm_vecs; pial_norm_vecs];
        else
            orig_surface=gifti(fullfile(surf_dir, sprintf('%s.gii', surface)));
            orig_idx=knnsearch(orig_surface.vertices,ds_surface.vertices);
            norm_vecs=spm_mesh_normals(struct('faces',orig_surface.faces,...
                'vertices',orig_surface.vertices),true);
            norm_vecs=norm_vecs(orig_idx,:);
        end
    case 'link_vector'
        switch surface
            case 'pial'
                ds_white=gifti(fullfile(surf_dir, 'white.ds.gii'));
                norm_vecs=ds_white.vertices-ds_surface.vertices;
            case 'white'
                ds_pial=gifti(fullfile(surf_dir, 'pial.ds.gii'));
                % Multiply by -1 so vectors point inward (toward inside of
                % brain)
                norm_vecs=(ds_pial.vertices-ds_surface.vertices).*-1;
            case 'white-pial'
                ds_pial=gifti(fullfile(surf_dir, 'pial.ds.gii'));
                ds_white=gifti(fullfile(surf_dir, 'white.ds.gii'));
                % Multiply white norms by -1 so vectors point inward
                % (toward inside of brain)
                norm_vecs=[(ds_pial.vertices-ds_white.vertices).*-1; ds_white.vertices-ds_pial.vertices];
        end        

    case 'variational'
        if strcmp(surface, 'white-pial')
            ds_white=gifti(fullfile(surf_dir, 'white.ds.gii'));
            orig_white=gifti(fullfile(surf_dir, 'white.gii'));
            orig_white_idx=knnsearch(orig_white.vertices,ds_white.vertices);
            lh_white_normals=dlmread(fullfile(surf_dir, 'lh.white.variational_normals.txt'));
            rh_white_normals=dlmread(fullfile(surf_dir, 'rh.white.variational_normals.txt'));
            white_norm_vecs=[lh_white_normals; rh_white_normals].*-1;        
            white_norm_vecs=white_norm_vecs(orig_white_idx,:);
            
            ds_pial=gifti(fullfile(surf_dir, 'pial.ds.gii'));
            orig_pial=gifti(fullfile(surf_dir, 'pial.gii'));
            orig_pial_idx=knnsearch(orig_pial.vertices,ds_pial.vertices);
            lh_pial_normals=dlmread(fullfile(surf_dir, 'lh.pial.variational_normals.txt'));
            rh_pial_normals=dlmread(fullfile(surf_dir, 'rh.pial.variational_normals.txt'));
            pial_normals=[lh_pial_normals; rh_pial_normals];        
            pial_normals=pial_normals(orig_pial_idx,:);
            
            norm_vecs=[white_norm_vecs; pial_normals];
        else 
            orig_surface=gifti(fullfile(surf_dir, sprintf('%s.gii',surface)));
            orig_idx=knnsearch(orig_surface.vertices,ds_surface.vertices);
            lh_norm_vecs=dlmread(fullfile(surf_dir, sprintf('lh.%s.variational_normals.txt', surface)));
            rh_norm_vecs=dlmread(fullfile(surf_dir, sprintf('rh.%s.variational_normals.txt', surface)));
            norm_vecs=[lh_norm_vecs; rh_norm_vecs];        
            % Multiply by -1 so vectors point inward (toward inside of
            % brain)
            if strcmp(surface,'white')
                norm_vecs=-1.*norm_vecs;
            end
            norm_vecs=norm_vecs(orig_idx,:);
        end         
        
        normN = sqrt(sum(norm_vecs.^2,2));
        normN(normN < eps)=1;
        norm_vecs = bsxfun(@rdivide,norm_vecs,normN);
        norm_vecs=double(norm_vecs);        
                
        % Replace where 0 with face normal
        norm_vecs2 = spm_mesh_normals(struct('faces',ds_surface.faces,...
            'vertices',ds_surface.vertices),true);
        z=find(sum(norm_vecs,2)==0);
        norm_vecs(z,:)=norm_vecs2(z,:);
end