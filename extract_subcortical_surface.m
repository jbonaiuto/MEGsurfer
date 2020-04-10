function extract_subcortical_surface(subj_id, structure_id, out_mesh_file, varargin)
% function extract_subcoortical_surface(subj_id, structure_id, out_mesh_file, varargin)
% out_mesh_file='hippl.gii';
% Read in aseg (automatic segmentation volume) file and extract a structure. Dec 2013
%% based on register_freesurf_mesh.m:
%% This code is to register a freesurfer defined mesh into the native space of the spm nifti file
% INPUT: 
%   subj_id: subject initials
%   structure_id: ID of structure in aseg (r hipp=53, l hipp=17, l thal=10)
%   out_mesh_file: name of mesh file to output ('hippl.gii')
%   subjects_dir: freesurfer's SUBJECT_DIR
% ---------------------------
% GRB Oct 2012
% JJB (j.bonaiuto@ucl.ac.uk) Jul 2015 - no need for coregistration - read and apply ras offset
% 

% Parse inputs
defaults = struct('subjects_dir',fullfile('/usr/local/freesurfer/subjects'));  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Read in aseg volume and find voxels labeled with the given structure_id
Vaseg=spm_vol(fullfile(params.subjects_dir, subj_id, 'mri', 'aseg.nii'));
[Yaseg,xyz]=spm_read_vols(Vaseg);
useind=find(Yaseg==structure_id);
coords=[xyz(1,useind)',xyz(2,useind)',xyz(3,useind)'];
indices=(pinv(Vaseg.mat)*[xyz(:,useind);ones(1,size(xyz(:,useind),2))])';

% Border of two voxels around object
border=2;
% No resampling
vstep=1;
% Limits of grid to create mesh in
Nx=[min(indices(:,1))-border:vstep:max(indices(:,1))+border];
Ny=[min(indices(:,2))-border:vstep:max(indices(:,2))+border];
Nz=[min(indices(:,3))-border:vstep:max(indices(:,3))+border];

% Create mesh grid of computed size
[X Y Z] = meshgrid(Nx, Ny, Nz);
% Create volume for structure
V=zeros(size(X));
for x1=Nx,
    for y1=Ny,
        for z1=Nz,
            if Yaseg(round(x1),round(y1),round(z1))==structure_id, %% note x and y switched,
                V(round(y1-min(Ny)+1),round(x1-min(Nx)+1),round(z1-min(Nz)+1))=1;
            end;
        end;
    end;
end;

% Construct surface in meshgrid from volume
FV = isosurface(X,Y,Z,V); 
figure;
trisurf(FV.faces,FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));
hold on;
plot3(indices(:,1),indices(:,2),indices(:,3),'wx');

%% now have the mesh defined in indices, want to convert these back to aseg mri coordinates
FVc=[];
FVc.faces=FV.faces;
FVc.vertices=(Vaseg.mat*[FV.vertices';ones(1,size(FV.vertices,1))])';
FVc.vertices=FVc.vertices(:,1:3);%+repmat(ras_offset,size(FVc.vertices,1),1);;
figure;
trisurf(FVc.faces,FVc.vertices(:,1),FVc.vertices(:,2),FVc.vertices(:,3));
hold on;
plot3(coords(:,1),coords(:,2),coords(:,3),'wx');

% Convert to gifti and save
gFV=gifti(FVc);
save(gFV,fullfile(params.subjects_dir, subj_id, 'surf', out_mesh_file));
trisurf(gFV.faces,gFV.vertices(:,1),gFV.vertices(:,2),gFV.vertices(:,3));

%figure;
%spmmesh =export(gifti(fullfile(subjects_dir, subj_id, 'surf', out_mesh_file)));
%trisurf(spmmesh.faces,spmmesh.vertices(:,1),spmmesh.vertices(:,2),spmmesh.vertices(:,3));
%hold on;
%spmmesh2 =export(gifti(fullfile(subjects_dir, subj_id, 'surf', 'white.gii')));
%trisurf(spmmesh2.faces,spmmesh2.vertices(:,1),spmmesh2.vertices(:,2),spmmesh2.vertices(:,3));

