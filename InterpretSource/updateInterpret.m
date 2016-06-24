function updateInterpret(savepath)

token = load([srcPrefix 'update.mat']);
optionsText = weboptions('username','jasmcole', 'password', token.token);
optionsBinary = weboptions('username','jasmcole', 'password', token.token, 'ContentType', 'binary');

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