function updateInterpret(savepath)

optionsText = weboptions('username','jasmcole', 'password', 'f56003177f78f893519e7eabc0a77957d3c86bc4');
optionsBinary = weboptions('username','jasmcole', 'password', 'f56003177f78f893519e7eabc0a77957d3c86bc4', 'ContentType', 'binary');

% Get list of source files
srcFiles = webread(['https://api.github.com/repos/jasmcole/Interpret/contents/InterpretSource'], optionsText);

hWait = waitbar(0, 'Downloading update...');

for n = 1:length(srcFiles)

    fname = srcFiles(n).name;

    srcFile = webread(['https://api.github.com/repos/jasmcole/Interpret/contents/InterpretSource/' fname], optionsText);
    dlURL = srcFile.download_url;
    fileContents = webread(dlURL, optionsBinary);
    waitbar(n/length(srcFiles), hWait, ['Downloaded ' fname])
    
    fid = fopen([savepath fname],'w');
    fwrite(fid, fileContents);
    fclose(fid);
end

delete(hWait)

end