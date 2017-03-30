function plot_surface_statistic(statistic, subj_info, varargin)

% Parse inputs
defaults = struct('surf_dir', 'd:\pred_coding\surf', 'clip_vals', true, 'limits', [], ...
    'plot_dir', 'C:\Users\jbonai\Dropbox\meg\layer_sim');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

orig_white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','white.hires.deformed.surf.gii');
white_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed.surf.gii');

orig_pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','pial.hires.deformed.surf.gii');
pial_mesh=fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed.surf.gii');

allmeshes=strvcat(white_mesh, pial_mesh);
Nmesh=size(allmeshes,1);

wm=gifti(white_mesh);
pial=gifti(pial_mesh);
wm_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_white.hires.deformed_inflated.surf.gii'));
pial_inflated=gifti(fullfile(params.surf_dir,[subj_info.subj_id subj_info.birth_date '-synth'],'surf','ds_pial.hires.deformed_inflated.surf.gii'));

pial_white_map=map_pial_to_white(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);
white_pial_map=map_white_to_pial(white_mesh, pial_mesh, ...
        'mapType', 'link', 'origPial', orig_pial_mesh, ...
        'origWhite', orig_white_mesh);

switch statistic
    case 'thickness'
        [pial_statistic, wm_statistic]=compute_thickness(pial_mesh, ...
            white_mesh, orig_pial_mesh, orig_white_mesh);
    case 'curvature'
        pial_statistic=compute_curvature(pial_mesh, 'curvature_type', 'mean');
        wm_statistic=compute_curvature(white_mesh, 'curvature_type', 'mean');
    case 'depth'
        [pial_statistic,HS]=compute_sulcal_depth(pial_mesh);
        mapping=dsearchn(HS.vertices,wm.vertices);
        wm_statistic=sqrt(sum((wm.vertices-HS.vertices(mapping,:)).^2,2));
    case 'lead_field_norm'
        D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_pial.mat');
        pial_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';

        D=spm_eeg_load('C:\Users\jbonai\Dropbox\meg\layer_sim\rgb_1_white.mat');
        wm_statistic=sqrt(sum(D.inv{1}.inverse.L.^2,1))';
end

switch statistic
    case 'thickness'
        fig=figure();
        bins=linspace(min(pial_statistic),max(pial_statistic),100);
        [n,b] = histc(pial_statistic,bins);
        bar(bins,n./sum(n)*100.0,'histc');
        xlabel('Thickness');
        ylabel('% Vertices');
    otherwise
        fig=figure();
        plot(pial_statistic,wm_statistic(pial_white_map),'.');
        hold on
        min_val=min([pial_statistic(:) wm_statistic(:)]);
        max_val=max([pial_statistic(:) wm_statistic(:)]);
        plot([min_val max_val],[min_val max_val],'r-');
        xlabel(sprintf('Pial %s',statistic));
        ylabel(sprintf('White matter %s',statistic));
end
%saveas(fig, fullfile(out_dir, sprintf('pial_white_%s.png', statistic)), 'png');
%figure2eps(fig, fullfile(out_dir, sprintf('pial_white_%s.eps', statistic)), 10, '-opengl');
        
[ax, metric_data]=plot_surface_metric(pial, pial_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('pial_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(pial_inflated, pial_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('pial_inflated_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(wm, wm_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('wm_%s.png', statistic)), 'png');

[ax, metric_data]=plot_surface_metric(wm_inflated, wm_statistic, 'clip_vals', params.clip_vals, 'limits', params.limits);
set(ax,'CameraViewAngle',5.338);
set(ax,'CameraTarget',[20.523 26.884 0.768]);
set(ax,'CameraUpVector',[-0.052 0.926 0.375]);
set(ax,'CameraPosition',[57.457 -626.649 1618.708]);
fig=get(ax,'Parent');
%saveas(fig, fullfile(out_dir, sprintf('wm_inflated_%s.png', statistic)), 'png');
