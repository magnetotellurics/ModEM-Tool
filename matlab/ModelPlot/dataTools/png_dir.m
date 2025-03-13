function png_dir(ext,dirname)
% png_dir(ext,dirname)
% ext = 'edi','xml','zrr','zmm','zss'
% if no dirname, use current directory
if nargin < 2
    dirname = '.';
end
fname = findfiles(ext,dirname);
for i = 1:length(fname)
    emtf_png(fname{i});
end
end
