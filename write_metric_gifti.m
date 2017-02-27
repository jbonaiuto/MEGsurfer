function write_metric_gifti(filename, data)
% WRITE_METRIC_GIFTI  Write data to a metric gifti file
%
% Use as
%   write_metric_gifti(filename, data)
% where the first argument is the name of the file to write to and the second
% argument is the data to write (NxN, where N=the number of vertices)

[filepath,filename,ext]=fileparts(filename);
c=file_array(fullfile(filepath,sprintf('%s.dat',filename)), size(data),'FLOAT32-LE',0,1,0);
c(:,:)=data;
g = gifti;
g.cdata = c;
save(g, fullfile(filepath,sprintf('%s.gii',filename)), 'ExternalFileBinary');
