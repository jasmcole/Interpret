function phase = RemoveHotPixelButton(handles)
phase = handles.phase;
phase = RemoveHotPixels(phase, 5, 1);
end