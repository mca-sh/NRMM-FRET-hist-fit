function fle = addnbtofilename(fle)

[src,fname,fext] = fileparts(fle);

outfname = [fname,fext];

number = 1;
while exist([src,filesep,outfname],'file')
    number = number + 1;
    outfname = [fname,'(' num2str(number) ')',fext];
end

fle = [src,filesep,outfname];