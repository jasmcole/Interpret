function I = ReadRAW16bit(filename)

bands = 1;
precision = 'uint16';
offset = 0;
interleave = 'bip'; % Also try bil or bip
byteorder = 'ieee-be';%Little endian
fileinfo = dir(filename);

if(fileinfo.bytes == 4456448)
    width = 2048;
    height = 1088;
    byteorder = 'ieee-le';%Little endian
elseif(fileinfo.bytes > 1228800)
    width = 1280;
    height = 960;
elseif(fileinfo.bytes < 1228800)
    width = 640;
    height = 480;
end

I = multibandread(filename, [height width bands], precision, offset, interleave, byteorder);
I = flipud(I);

end