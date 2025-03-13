function emtf_png(filename)
[~,~,myext] = fileparts(filename);
if contains(myext,'xml')
    format = 'XML';
elseif contains(myext,'edi')
    format = 'EDI';
elseif contains(myext,'z')
    format = 'Z';
else
    error(['Unknown EMTF file format: ',filename]);
end
tf = mttf.read(filename,format,'Full_Impedance');
tf.apresplt;
fname=filename;
oname=strrep(fname,myext,'.png');
fprintf('%s \n',oname);
com=['print -dpng ''' oname ''''];
eval(com);close all;
return
