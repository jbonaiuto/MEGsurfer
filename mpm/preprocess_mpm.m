function preprocess_mpm(data_dir, varargin)

% Parse inputs
defaults = struct('untar',true,'convert_dicom',true,'create_mpm',true,...
    'robust_bias_correction',true,'align_to_headcast_t1',true,...
    'headcast_t1','','mt_prefix','mtw_mfc_3dflash','pd_prefix','pdw_mfc_3dflash',...
    't1_prefix','t1w_mfc_3dflash','b1_map_prefix','mfc_seste_b1map',...
    'gre_map_prefix','gre_field_mapping');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)'
    if ~isfield(params, f{1})
        params.(f{1}) = defaults.(f{1});
    end
end

%%% Finding tar files and untarring them automatically
if params.untar
    h = waitbar(0,'Untarring ...');
    fileList = dir(fullfile(data_dir, '*.tar'));
    if ~isempty(fileList) % if there are tar files in the determined directory
        for i = 1 : length(fileList)
            untar(fullfile(data_dir, fileList(i).name),...
                fullfile(data_dir, fileList(i).name(1:end-4)))
            waitbar(i/length(fileList))
        end
        close(h)
    end
end

spm('defaults', 'EEG');

mpm_dir=fullfile(data_dir, 'mpm'); % The resulted files will be stored in 'mpm' folder

%%% Import DICOM
if params.convert_dicom
    spm_jobman('initcfg'); 
    clear jobs
    matlabbatch={};
    batch_idx=1;
    matlabbatch{batch_idx}.spm.util.import.dicom.data={};
    [~,dirs] = spm_select('List', data_dir);
    for d=1:size(dirs,1)
        if dirs(d,1)=='M'
            scan_dir=deblank(dirs(d,:));
            [scan_files, ~] = spm_select('List', fullfile(data_dir,scan_dir));
            for f=1:size(scan_files,1)
                filename=deblank(scan_files(f,:));
                matlabbatch{batch_idx}.spm.util.import.dicom.data{end+1}=fullfile(data_dir,scan_dir,filename);
            end
        end
    end
    matlabbatch{batch_idx}.spm.util.import.dicom.data=matlabbatch{batch_idx}.spm.util.import.dicom.data';
    matlabbatch{batch_idx}.spm.util.import.dicom.root = 'series';
    matlabbatch{batch_idx}.spm.util.import.dicom.outdir = {data_dir};
    matlabbatch{batch_idx}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{batch_idx}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{batch_idx}.spm.util.import.dicom.convopts.icedims = 0;
    spm_jobman('run',matlabbatch);
end

%%% Create MPMs
if params.create_mpm
    spm_jobman('initcfg'); 
    clear jobs
    matlabbatch={};
    batch_idx=1;

    if exist(mpm_dir,'dir')~=7
        mkdir(mpm_dir);
    end

    matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.output.outdir = {mpm_dir};
    if ~params.robust_bias_correction
        matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';
    else
        Maps = {params.mt_prefix, 'raw_sens_MT' ; 
                params.pd_prefix, 'raw_sens_PD';
                params.t1_prefix, 'raw_sens_T1'};
        
        for i = 1 : size(Maps , 1)
            matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.(Maps{i , 2})={};
            % Figure out scan number of map
            fileList=dir(fullfile(data_dir, [Maps{i,1} '*']));
            if ~isempty(fileList)
                name_parts=strsplit(fileList(1).name,'_');
                scan_num=str2num(name_parts(end));
                fileList=dir(fullfile(data_dir, sprintf('%s*%d','mfc_smaps_v1a_Array', scan_num-2)));
                if ~isempty(fileList)
                    for ii=1:length(fileList)
                        seq_dir=fullfile(data_dir, fileList(ii).name);
                        [files,~]=spm_select('List', seq_dir);
                        for f=1:size(files,1)
                            filename=deblank(files(f,:));
                            matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.(Maps{i , 2}){end+1}=fullfile(seq_dir,filename);
                        end
                    end
                end
                fileList=dir(fullfile(data_dir, sprintf('%s*%d','mfc_smaps_v1a_QBC', scan_num-1)));
                if ~isempty(fileList)
                    for ii=1:length(fileList)
                        seq_dir=fullfile(data_dir, fileList(ii).name);
                        [files,~]=spm_select('List', seq_dir);
                        for f=1:size(files,1)
                            filename=deblank(files(f,:));
                            matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.(Maps{i , 2}){end+1}=fullfile(seq_dir,filename);
                        end
                    end
                end
                matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.(Maps{i , 2})=matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_per_contrast.(Maps{i , 2})';
            end
        end
    end
    
    % Set B1 and B0 input
    Maps = {params.b1_map_prefix, 'b1input' ; 
            params.gre_map_prefix, 'b0input'};
        
    for i = 1 : size(Maps , 1)
        matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.(Maps{i , 2})= {}; % Initializing
        fileList=dir(fullfile(data_dir, [Maps{i , 1} , '*']));
        if ~isempty(fileList)
            for ii=1:length(fileList)
                seq_dir=fullfile(data_dir, fileList(ii).name);
                [files,~]=spm_select('List', seq_dir);
                for f=1:size(files,1)
                    filename=deblank(files(f,:));
                    matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.(Maps{i , 2}){end+1}=fullfile(seq_dir,filename);
                end
            end
            matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.(Maps{i , 2})=matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.(Maps{i , 2})';
        end
    end
    
    matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1parameters.b1metadata = 'yes';
    
    Maps = {params.mt_prefix, 'MT'; 
            params.pd_prefix, 'PD';
            params.t1_prefix, 'T1'};
        
    % Set map input
    for i = 1 : size(Maps , 1)
        matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.raw_mpm.(Maps{i , 2})= {}; % Initializing
        fileList=dir(fullfile(data_dir, [Maps{i , 1} , '*']));
        if ~isempty(fileList)
            for ii=1:length(fileList)
                seq_dir=fullfile(data_dir, fileList(ii).name);
                [files,~]=spm_select('List', seq_dir);
                for f=1:size(files,1)
                    filename=deblank(files(f,:));
                    matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.raw_mpm.(Maps{i , 2}){end+1}=fullfile(seq_dir,filename);
                end
            end
            matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.raw_mpm.(Maps{i , 2})=matlabbatch{batch_idx}.spm.tools.hmri.create_mpm.subj.raw_mpm.(Maps{i , 2})';
        end
    end        
    
    spm_jobman('run',matlabbatch);
end

if params.align_to_headcast_t1
    spm_jobman('initcfg'); 
    clear jobs
    matlabbatch={};
    batch_idx=1;
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.ref = {params.headcast_t1};
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.source = {};
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.other = {};
    [files,~] = spm_select('List', fullfile(mpm_dir,'Results'));
    for f=1:size(files,1)
        filename=deblank(files(f,:));
        if contains(filename,'_R1.nii')
            matlabbatch{batch_idx}.spm.spatial.coreg.estimate.source = {fullfile(mpm_dir,'Results',[filename ',1'])};
        elseif contains(filename,'_MT.nii') || contains(filename,'_PD.nii') || contains(filename,'_R1.nii') || contains(filename,'_R2s_OLS.nii')
            matlabbatch{batch_idx}.spm.spatial.coreg.estimate.other{end+1}=fullfile(mpm_dir,'Results',[filename ',1']);
        end
    end
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.other=matlabbatch{batch_idx}.spm.spatial.coreg.estimate.other';
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{batch_idx}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    spm_jobman('run',matlabbatch);
end
