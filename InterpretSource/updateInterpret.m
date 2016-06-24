function updateInterpret(savepath)

optionsText = weboptions('username','jasmcole', 'password', 'dd20afd5570e01e418a0e59ea906b0e8e71f2652');
optionsBinary = weboptions('username','jasmcole', 'password', 'dd20afd5570e01e418a0e59ea906b0e8e71f2652', 'ContentType', 'binary');

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