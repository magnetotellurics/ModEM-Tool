function xml_png_dir(dirname);
[a,b]=unix(['ls -1 ' dirname '/*.xml']);
ii=findstr(b,'.xml');nf=length(ii);
i1=1;
for k=1:nf
 i2=ii(k)+3;
 fnames{k}=b(i1:i2);
 i1=i2+2;
end
for k=1:nf
 fname=char(fnames{k});
 try
  xml_png(fname);
 catch
 fprintf('ERROR: INCOMPLETE XML %s\n',fname);
 end
end
return
end
