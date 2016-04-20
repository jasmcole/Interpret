function success = saveGUIState(handles)
    guistate = struct();
    try
        guistate.experiment = get(handles.ExpPop,'Value');
        guistate.calibration = get(handles.CalibPop,'Value');
        guistate.phasemethod = get(handles.PhasePop,'Value');
        guistate.densitymethod = get(handles.DensityPop,'Value');
        guistate.medium = get(handles.MediumPop,'Value');
        guistate.month = str2num(get(handles.MonthBox, 'String'));
        guistate.day   = str2num(get(handles.DayBox, 'String'));
        guistate.run   = str2num(get(handles.RunBox, 'String'));
        guistate.shot  = str2num(get(handles.ShotBox, 'String'));
        guistate.year  = str2num(get(handles.yearBox, 'String'));
        guistate.wavelength = str2num(get(handles.WavelengthBox, 'String'));
        save([handles.introot filesep 'GUIstate.mat'], 'guistate')
        success = true;
    catch
        success = false;
    end
    
    if(~success)
        choice = questdlg('Couldn''t save GUI state. Close anyway?', 'Close dialog', 'Yes','No','No');
        
        switch choice
            case 'Yes'
                success = true;
            case 'No'
                success = false;
        end
    end
end