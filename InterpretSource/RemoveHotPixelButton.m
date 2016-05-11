function phase2 = RemoveHotPixelButton(handles)
phase = handles.phase;

% Median filtered phase
pmed = medfilt2(phase, [3 3]);

% Get Laplacians of original and median filtered phases
pl = del2(phase);
pmedl = del2(pmed);

% Compute relative 'spikiness' of phase image
% Add regularisation parameter 1 to avoid division by zeros
spikiness = abs((1 + pl)./(1 + pmedl));

% Threshold spikiness to replace pixels
level = 1;
phase2 = phase;
phase2(spikiness > level) = pmed(spikiness > level);

end