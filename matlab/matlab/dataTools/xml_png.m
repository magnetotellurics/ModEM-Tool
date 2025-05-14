function xml_png(filename)
tf = mttf.read(filename,'xml','Full_Impedance');
tf.apresplt;
fname=filename;
oname=strrep(fname,'.xml','.png');
fprintf('%s \n',oname);
com=['print -dpng ' oname];
eval(com);close all;
return
