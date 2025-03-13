function [result]  = fourOReight(cfile);
%   returns length (4 or 8) in bytes of reacord header length
%     for sequential binary fortran file cfile; if neither case
%     appears correct (e.g., not a sequential fortran binary file)
%     returns -1

fid = fopen(cfile);
lrec = fread(fid,1,'integer*4');
status = fseek(fid,lrec,'cof');
if status == -1
    result = -1;
    fclose(fid);
    return
end
nrec = fread(fid,1,'integer*4');
if nrec==lrec
    result = 4;
    fclose(fid);
    return
else
    %   rewind file, try 8 bytes
    fseek(fid,0,-1);
    lrec = fread(fid,1,'integer*8');
    status = fseek(fid,lrec,'cof');
    if status == -1
        result = -1;
        fclose(fid);
        return
    end
    nrec = fread(fid,1,'integer*8');
    if nrec==lrec
        result = 8;
    else
        result = -1;
    end
end
fclose(fid);