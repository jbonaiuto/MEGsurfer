function curvature=compute_curvature(surface_file, varargin)
% COMPUTE_CURVATURE  Compute curvature for each vertex in a surface
%
% Use as
%   curvature=compute_curvature(surface_file)
% where the first argument is the surface filename. Returns
% a vector of curvature values with an element for each vertex
%
%   compute_curvature(...,'param','value','param','value'...) allows
%    additional param/value pairs to be used. Allowed parameters:
%    * curvature_type - 'k1' (default) or string - type of curvature to compute
%       - k1
%       - k2
%       - mean
%       - gaussian
%       - ici
%       - nici
%       - gln
%       - aici
%       - mci
%       - nmci
%       - mln
%       - amci
%       - fi
%       - ci

% Parse inputs
defaults = struct('curvature_type', 'k1');  %define default values
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

addpath('C:\Users\jbonai\Documents\MATLAB\curvature');
[path file ext]=fileparts(surface_file);
k1_file=fullfile(path, sprintf('%s_k1',file));
k2_file=fullfile(path, sprintf('%s_k2',file));

% Check if k1 and k2 files already exist - if not, recompute
if exist(sprintf('%s.shape.gii',k1_file),'file')~=2 || exist(sprintf('%s.shape.gii',k2_file),'file')~=2
    surface=gifti(surface_file);
    %% calculate curvatures
    getderivatives=0;
    [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( surface ,getderivatives);

    k1=PrincipalCurvatures(1,:)';
    k2=PrincipalCurvatures(2,:)';
    write_metric_gifti(k1_file, k1);
    write_metric_gifti(k2_file, k2);
% Otherwise load from files
else
    k1_surf=gifti(sprintf('%s.shape.gii',k1_file));
    k1=k1_surf.cdata(:);
    k2_surf=gifti(sprintf('%s.shape.gii',k2_file));
    k2=k2_surf.cdata(:);
end

K=k1.*k2;
H=mean([k1 k2],2);
% Depending on curvature type - compute curvature metric
switch params.curvature_type
    case 'k1'
        curvature=k1;
    case 'k2'
        curvature=k2;
    case 'mean'
        curvature=H;
    case 'gaussian'
        curvature=K;
    case 'ici'        
        curvature=max(K,0);
    case 'nici'
        curvature=min(K,0);
    case 'gln'
        curvature=K.*K;
    case 'aici'
        curvature=abs(K);
    case 'mci'
        curvature=max(H,0);
    case 'nmci'
        curvature=min(H,0);
    case 'mln'
        curvature=H.*H;
    case 'amci'
        curvature=abs(H);
    case 'fi'
        curvature=abs(k1).*(abs(k1)-abs(k2));
    case 'ci'
        curvature=sqrt((k1.*k1+k2.*k2)/2.0);
end