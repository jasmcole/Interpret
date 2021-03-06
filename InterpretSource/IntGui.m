function varargout = IntGui(varargin)
% INTGUI MATLAB code for IntGui.fig
%      INTGUI, by itself, creates a new INTGUI or raises the existing
%      singleton*.
%
%      H = INTGUI returns the handle to a new INTGUI or the handle to
%      the existing singleton*.
%
%      INTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTGUI.M with the given input arguments.
%
%      INTGUI('Property','Value',...) creates a new INTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IntGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IntGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IntGui

% Last Modified by GUIDE v2.5 23-Sep-2017 13:03:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @IntGui_OpeningFcn, ...
    'gui_OutputFcn',  @IntGui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before IntGui is made visible.
function IntGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IntGui (see VARARGIN)

% Choose default command line output for IntGui
handles.output = hObject;

mroot = fileparts(which('Interpret.m'));
introot = [mroot filesep 'Interpret'];
handles.introot = introot;
srcPrefix = [handles.introot filesep 'InterpretSource' filesep];

% Update handles structure
guidata(hObject, handles);

try
    dategithub = GetGithubRepoLastCommitTime('Interpret', srcPrefix);

    if(exist([introot filesep 'dateinstalled.mat']) == 2)
       dateinstalled = load([introot filesep 'dateinstalled.mat']);
       dateinstalled = dateinstalled.dateinstalled;
    else
       dateinstalled = datestr(now, 'yyyy-mm-dd');
       save([introot filesep 'dateinstalled.mat'], 'dateinstalled') 
    end

    if(datetime(dateinstalled) < datetime(dategithub))
        button = questdlg('There is a newer version of Interpret. Download now?', 'Interpret Update', 'Yes', 'Not now', 'Skip this version', 'Yes');
        if(strcmp(button, 'Yes'))
            button = questdlg(['Installing files to ' srcPrefix '. Is this correct?']);
            if(strcmp(button, 'Yes'))
                updateInterpret(srcPrefix);
                dateinstalled = datestr(now, 'yyyy-mm-dd');
                save([introot filesep 'dateinstalled.mat'], 'dateinstalled')
                UpdateStatus('Update successful. Please restart Interpret.', handles)
            end
        elseif(strcmp(button, 'Skip this version'))
            dateinstalled = datestr(now, 'yyyy-mm-dd');
            save([introot filesep 'dateinstalled.mat'], 'dateinstalled') 
        else      
            UpdateStatus(['Source last updated ' dategithub char(10) 'Installation date ' dateinstalled], handles)
            
        end
    else
        UpdateStatus(['Source last updated ' dategithub char(10) 'Installation date ' dateinstalled], handles)
    end

catch
    UpdateStatus('There was an authentication error with the Interpret Github repository.', handles)
end
    
try
    load([handles.introot filesep 'GUIstate.mat'])
    set(handles.yearBox,'String', guistate.year)
    set(handles.MonthBox,'String', guistate.month)
    set(handles.DayBox,'String', guistate.day)
    set(handles.RunBox,'String', guistate.run)
    set(handles.ShotBox,'String', guistate.shot)
    set(handles.WavelengthBox,'String', guistate.wavelength)
    set(handles.ExpPop, 'Value', guistate.experiment)
    set(handles.CalibPop, 'Value', guistate.calibration)
    set(handles.PhasePop, 'Value', guistate.phasemethod)
    set(handles.DensityPop, 'Value', guistate.densitymethod)
    set(handles.MediumPop, 'Value', guistate.medium)
catch
    UpdateStatus('There was an error loading the GUI state. If this error persists try deleting GUIstate.m .', handles)
end

try
    updateCalibPop(handles.CalibPop, handles);
catch
    UpdateStatus('There was an error loading the Calibration Database.', handles)
end



function varargout = IntGui_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
success = saveGUIState(handles);
if(success)
    delete(handles.figure1)
end


% --- Executes on selection change in ExpPop.
function ExpPop_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ExpPop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

data = read_mixed_csv('ExperimentDatabase.csv', ',');
[nrecords, nparams] = size(data);
exps = cell(nrecords-1,1);
for n = 2:nrecords
    exps{n-1} = data{n,1};
end
set(hObject, 'String', exps);

function MonthBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function MonthBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DayBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function DayBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RunBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function RunBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ShotBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ShotBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalibBut.
function CalibBut_Callback(hObject, eventdata, handles)
set(handles.StatusBox,'String','Reading calibration'); drawnow
calibration = feval(@(x) x{1}{x{2}},get(handles.CalibPop,{'String','Value'}));
calibdata = CalibrationDatabase(calibration);
set(handles.StatusBox,'String','Reading calibration'); drawnow
if(isfield(calibdata, 'Name'))
    set(handles.CalibCheck, 'Value',1,'String','Found calibration', 'ForegroundColor', 'green')
    set(handles.CalibrationTable, 'RowName', fieldnames(calibdata))
    set(handles.CalibrationTable, 'ColumnName', 'Calibration Data')
    
    members = fieldnames(calibdata);
    params = cell(length(members),1);
    for n = 1:length(members)
        params{n} = calibdata.(members{n});
    end
    set(handles.CalibrationTable, 'Data', params)
    
    handles.calibdata = calibdata;
    guidata(hObject,handles)
    
    axes(handles.InterferogramAxes);
    I = handles.originaldataimage;
    I = imrotate(I, calibdata.rotation);
    
    if (calibdata.x == 0)
        calibdata.x = 1;
    end
    if (calibdata.y == 0)
        calibdata.y = 1;
    end
    if (calibdata.w == 0)
        calibdata.w = length(I(1,:))-1;
    end
    if (calibdata.h == 0)
        calibdata.h = length(I(:,1))-1;
    end
    
    try
        I = I(calibdata.y:calibdata.y+calibdata.h, calibdata.x:calibdata.x+calibdata.w);
    catch
        calibdata.x = 1;
        calibdata.y = 1;
        calibdata.w = length(I(1,:))-1;
        calibdata.h = length(I(:,1))-1;
        uiwait(msgbox('Failed to apply ROI. Reset ROI by clicking one of x, y, w, h.', 'Warning'));
    end

    imagesc(I)
    axis image xy
    caxis([0 2*mean(mean(I))])
    handles.dataimage = I;
    guidata(hObject,handles)
    
    axes(handles.ReferenceAxes);
    UpdateStatus('Loading reference', handles)
        
    try
        if (strcmp(calibdata.reference(end-2:end), 'raw'))
            I = ReadRAW16bit(calibdata.reference);
        else
            I = imread(calibdata.reference);
            if(length(size(I)) > 2)
                I = rgb2gray(I);
            end
        end
        I = I - min(min(I));
        UpdateStatus('Reference loaded', handles)
    catch
        I = ones(size(handles.originaldataimage));
        UpdateStatus('Reference not found', handles)
    end
    
    I = imrotate(I, calibdata.rotation);
    I = I(calibdata.y:calibdata.y+calibdata.h, calibdata.x:calibdata.x+calibdata.w);
    
    handles.refimage = I;
    guidata(hObject,handles)
    
    imagesc(I); axis image xy
    caxis([0 2*mean(mean(I))])
    drawnow
    
end
if(~isfield(calibdata, 'Name'))
    set(handles.CalibCheck, 'Value',0,'String','No calibration', 'ForegroundColor', 'red')
end

% --- Executes on button press in FileBut.
function FileBut_Callback(hObject, eventdata, handles)

experiment = feval(@(x) x{1}{x{2}},get(handles.ExpPop,{'String','Value'}));
month = str2num(get(handles.MonthBox, 'String'));
day   = str2num(get(handles.DayBox, 'String'));
run   = str2num(get(handles.RunBox, 'String'));
shot  = str2num(get(handles.ShotBox, 'String'));
year  = str2num(get(handles.yearBox, 'String'));

data = read_mixed_csv('ExperimentDatabase.csv', ',');
[nrecords nparams] = size(data);

for n = 2:nrecords
    if(strcmp(experiment, data{n,1}))
        fileexpression = data{n,2};
    end
end

filename = ParseExperimentPath(fileexpression, year, month, day, run, shot);

handles.datafile = filename;
guidata(hObject,handles)
[filenamewrap,newpos] = textwrap(handles.FileTxt,{filename}, 50);
set(handles.FileTxt,'String',filenamewrap, 'Position', newpos)
if (exist(filename, 'file') == 2)
    set(handles.FileCheck, 'Value',1,'String','Found file', 'ForegroundColor', 'green')
    set(handles.CalibCheck, 'Value',0,'String','Reload calibration?', 'ForegroundColor', 'blue')
end
if (exist(filename, 'file') ~= 2)
    set(handles.FileCheck,'Value',0,'String','No file - Drobo mounted?', 'ForegroundColor', 'red')
    return
end

axes(handles.InterferogramAxes);
UpdateStatus('Reading data file', handles)

if (strcmp(filename(end-2:end), 'raw'))
    I = ReadRAW16bit(filename);
else
    I = imread(filename);
    if(length(size(I)) > 2)
        I = rgb2gray(I);
    end
end
I(isnan(I)) = 0;
I = I - min(min(I));

handles.dataimage = I;
handles.originaldataimage = I;
guidata(hObject,handles)

colormap(gray)
caxis([0 2*mean(mean(I))])
UpdateStatus('Data file loaded', handles)
if(isfield(handles, 'calibdata'))
    CalibBut_Callback(gcbo, eventdata, handles)
else
    imagesc(I); axis image xy
end

% --- Executes on button press in FileCheck.
function FileCheck_Callback(hObject, eventdata, handles)

% --- Executes on button press in CalibCheck.
function CalibCheck_Callback(hObject, eventdata, handles)

function FileBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function FileBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WavelengthBox_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function WavelengthBox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in MediumPop.
function MediumPop_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function MediumPop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PhasePop.
function PhasePop_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function PhasePop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DensityPop.
function DensityPop_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function DensityPop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PhaseBut.
function PhaseBut_Callback(hObject, eventdata, handles)

RetrievalType = feval(@(x) x{1}{x{2}},get(handles.PhasePop,{'String','Value'}));
axes(handles.PhaseAxes)

switch RetrievalType
    case 'FFT'
        UpdateStatus('Beginning FFT phase retrieval', handles)
        phase = PhaseRetrieve(handles);
    case 'CWT'
        UpdateStatus('Beginning CWT phase retrieval', handles)
        [phase, mask] = CWTPhaseRetrieve(handles);
        UpdateStatus('Finished CWT phase retrieval', handles)
    case 'Hilbert'
        UpdateStatus('Beginning Hilbert phase retrieval', handles)
        phase = HilbertPhaseRetrieve(handles);
        UpdateStatus('Finished Hilbert phase retrieval', handles)
    case 'CWT2D'
        set(handles.StatusBox, 'String', 'Beginning CWT2D phase retrieval'); drawnow
        [phase, mask] = CWT2DPhaseRetrieve(handles);
        set(handles.StatusBox, 'String', 'Finshed CWT2D phase retrieval'); drawnow
end

if (mean(phase(:)) < 0)
    phase = -phase;
end

phase(isnan(phase)) = 0;

if(~exist('mask', 'var'))
    mask = ones(size(phase));
end

handles.phase = phase.*mask;
handles.mask = mask;

guidata(hObject, handles)


% --- Executes on button press in DensityBut.
function DensityBut_Callback(hObject, eventdata, handles)
% hObject    handle to DensityBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calibdata = handles.calibdata;
lambda = str2num(get(handles.WavelengthBox, 'String'));
phase = handles.phase;
medium = feval(@(x) x{1}{x{2}},get(handles.MediumPop,{'String','Value'}));
rhoplotflag = 1;
axes(handles.DensityAxes)
RetrievalType = feval(@(x) x{1}{x{2}},get(handles.DensityPop,{'String','Value'}));
switch RetrievalType
    case 'Abel Inversion'
        [rho,xaxis,yaxis,rhomean] = AbelInversion(phase, calibdata, medium, lambda, rhoplotflag, handles);
    case 'Trapezoidal Fit'
        rho = FitToPlateau(phase, calibdata, medium, lambda, rhoplotflag);
    case 'Moire Abel Inversion'
        [rho,xaxis,yaxis,rhomean] = AbelInversionMoire(phase, calibdata, medium, lambda, rhoplotflag, handles);
end

assignin('base','rho',rho)
assignin('base','phase',phase)
assignin('base','xaxis',xaxis)
assignin('base','yaxis',yaxis)

handles.rho = rho;
handles.rhomean = rhomean;
handles.xaxis = xaxis;
handles.yaxis = yaxis;
guidata(hObject, handles)

% --- Executes on button press in OpenBtn.
function OpenBtn_Callback(hObject, eventdata, handles)
% hObject    handle to OpenBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.*', 'Open image file');

if(filename ~= 0)
    filename = [pathname filename];

    handles.datafile = filename;
    guidata(hObject,handles)
    [filenamewrap,newpos] = textwrap(handles.FileTxt,{filename}, 50);
    set(handles.FileTxt,'String',filenamewrap, 'Position', newpos)
    if (exist(filename, 'file') == 2)
        set(handles.FileCheck, 'Value',1,'String','Found file', 'ForegroundColor', 'green')
        set(handles.CalibCheck, 'Value',0,'String','Reload calibration?', 'ForegroundColor', 'blue')
    end
    if (exist(filename, 'file') ~= 2)
        set(handles.FileCheck,'Value',0,'String','No file', 'ForegroundColor', 'red')
    end

    axes(handles.InterferogramAxes);
    set(handles.StatusBox,'String','Reading data file'); drawnow

    if(strcmp(filename(end-2:end), 'raw'))
        I = ReadRAW16bit(filename);
    else
        I = imread(filename);
        if(length(size(I)) > 2)
            I = rgb2gray(I);
        end
    end
    I = I - min(I(:));
    imagesc(I)
    colormap jet
    caxis([0 2*mean(mean(I))])
    handles.dataimage = I;
    handles.originaldataimage = I;
    guidata(hObject,handles)
    set(handles.StatusBox,'String','Data file loaded'); drawnow
end


% --- Executes on button press in SavecalibBut.
function SavecalibBut_Callback(hObject, eventdata, handles)
% hObject    handle to SavecalibBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UpdateCalibrationDatabase(handles);
set(handles.StatusBox,'String','Applying updated calibration'); drawnow
CalibBut_Callback(hObject, eventdata, handles)
set(handles.StatusBox,'String','Calibration applied'); drawnow


% --- Executes on button press in GetavBut.
function GetavBut_Callback(hObject, eventdata, handles)
% hObject    handle to GetavBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xaxis = handles.xaxis;
yaxis = handles.yaxis;
rhomean = handles.rhomean;
set(handles.StatusBox, 'String', 'Click between two points')
axes(handles.DensitydiagAxes)
[x y] = ginput(2);
[val, xind1] = min(abs(xaxis - x(1)));
[val, xind2] = min(abs(xaxis - x(2)));
rhoav = mean(rhomean(xind1:xind2))/1e24;
rhoerr = std(rhomean(xind1:xind2))/1e24;
hold on
line([min(xaxis) max(xaxis)], [rhoav*1e24 rhoav*1e24])
fill([x(1) x(2) x(2) x(1)], [1e24*(rhoav-rhoerr) 1e24*(rhoav-rhoerr) 1e24*(rhoav+rhoerr) 1e24*(rhoav+rhoerr)], 'r')
hold off
set(gca,'children',flipud(get(gca,'children')))

assignin('base','rhoav',rhoav)
assignin('base','rhoerr',rhoerr)
set(handles.StatusBox, 'String', ['Average density/10^18 ' num2str(rhoav) ' with error ' num2str(rhoerr)])


% --- Executes on selection change in CalibPop.
function CalibPop_Callback(hObject, eventdata, handles)
% hObject    handle to CalibPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CalibPop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CalibPop


% --- Executes during object creation, after setting all properties.
function CalibPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalibPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateCalibPop(hObject, handles)
data = read_mixed_csv([handles.introot filesep 'CalibrationDatabase.csv'], ',');
[nparams, nrecords] = size(data);
exps = cell(nrecords-1,1);
for n = 2:nrecords
    exps{n-1} = data{1,n};
end
set(hObject, 'String', exps);


% --- Executes on button press in unwrapTopBtn.
function unwrapTopBtn_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapTopBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = handles.phase';
handles.phase = unwrap(handles.phase);
handles.phase = handles.phase';
handles.phase = flipud(unwrap(flipud(handles.phase)));
handles.phase(handles.mask == 0) = 0;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in unwrapBottomBtn.
function unwrapBottomBtn_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapBottomBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = handles.phase';
handles.phase = unwrap(handles.phase);
handles.phase = handles.phase';
handles.phase = unwrap(handles.phase);
handles.phase(handles.mask == 0) = 0;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)

% --- Executes on button press in unwrapLeftBtn.
function unwrapLeftBtn_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapLeftBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = flipud(unwrap(flipud(handles.phase)));
handles.phase = handles.phase';
handles.phase = unwrap(handles.phase);
handles.phase = handles.phase';
handles.phase(handles.mask == 0) = 0;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)

% --- Executes on button press in unwrapRightBtn.
function unwrapRightBtn_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapRightBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = flipud(unwrap(flipud(handles.phase)));
handles.phase = imrotate(handles.phase, 90);
handles.phase = unwrap(handles.phase);
handles.phase = imrotate(handles.phase, -90);
handles.phase(handles.mask == 0) = 0;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)

% --- Executes on button press in savePhaseBtn.
function savePhaseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to savePhaseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
output = struct;
output.phase = handles.phase;
output.datafile = handles.datafile;
output.reffile = handles.calibdata.reference;
if(isfield(handles, 'rho'))
    output.rho = handles.rho;
end
if(isfield(handles, 'xaxis'))
    output.xaxis = handles.xaxis;
    output.yaxis = handles.yaxis;
end
assignin('base','InterpretData',output)
guidata(hObject, handles)


% --- Executes on button press in savePhaseFileBtn.
function savePhaseFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to savePhaseFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uiputfile('phase.mat', 'Save as mat file');
InterpretData = struct;
InterpretData.phase = handles.phase;
InterpretData.datafile = handles.datafile;
InterpretData.reffile = handles.calibdata.reference;
InterpretData.rho = handles.rho;
InterpretData.xaxis = handles.xaxis;
InterpretData.yaxis = handles.yaxis;
save([pathname filename], 'InterpretData')
set(handles.StatusBox,'String',['Saved data to ' pathname filename]); drawnow


% --- Executes on button press in rLeftBtn.
function rLeftBtn_Callback(hObject, eventdata, handles)
% hObject    handle to rLeftBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = imrotate(handles.phase, -90);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in rRightBtn.
function rRightBtn_Callback(hObject, eventdata, handles)
% hObject    handle to rRightBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = imrotate(handles.phase, 90);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in goldsteinBtn.
function goldsteinBtn_Callback(hObject, eventdata, handles)
% hObject    handle to goldsteinBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StatusBox,'String','Select point of known phase'); drawnow
handles.phase = GoldsteinUnwrap2D(handles.phase, ones(size(handles.phase)), handles);
handles.phase = handles.phase .* handles.mask;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes when entered data in editable cell(s) in CalibrationTable.
function CalibrationTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to CalibrationTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if(get(handles.autoUpdateCheck, 'Value'))
    SavecalibBut_Callback(hObject, eventdata, handles)
else
    set(handles.CalibCheck, 'Value',0,'String','Update calibration!', 'ForegroundColor', 'blue')
end


% --- Executes on button press in autoUpdateCheck.
function autoUpdateCheck_Callback(hObject, eventdata, handles)
% hObject    handle to autoUpdateCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoUpdateCheck


% --- Executes on button press in newCalibBtn.
function newCalibBtn_Callback(hObject, eventdata, handles)
% hObject    handle to newCalibBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = read_mixed_csv([handles.introot filesep 'CalibrationDatabase.csv'], ',');
[nparams nrecords] = size(data);
exps = cell(nrecords-1,1);
for n = 2:nrecords
    exps{n-1} = data{1,n};
end

newcalfig = figure();
newcalfig.MenuBar = 'none';
newcalfig.ToolBar = 'none';
newcalfig.Name = 'Copy calibration';
newcalfig.NumberTitle = 'off';
pos = get(newcalfig, 'Position');
width = 300;
height = 200;
set(newcalfig, 'Position', [pos(1) pos(2) width height])
newcaltxt = uicontrol('Style', 'text', 'parent', newcalfig, 'string', 'Name of new calibration', 'FontSize', 14, 'HorizontalAlignment', 'left');
newcaltxt.Position = [10 height-30 180 20];
newcaledt = uicontrol('Style', 'edit', 'parent', newcalfig, 'FontSize', 14, 'HorizontalAlignment', 'left');
newcaledt.Position = [11 height-55 180 20];
newcaltxt = uicontrol('Style', 'text', 'parent', newcalfig, 'string', 'Copy from', 'FontSize', 14, 'HorizontalAlignment', 'left');
newcaltxt.Position = [10 height-90 180 20];
newcalpop = uicontrol('Style', 'popupmenu', 'parent', newcalfig, 'string', exps, 'FontSize', 14);
newcalpop.Position = [7 height-100 width-20 1];
newcalbtn = uicontrol('Style', 'pushbutton', 'parent', newcalfig, 'string', 'OK', 'FontSize', 14, 'Callback', @newCalibOK_Callback, 'UserData', handles);
newcalbtn.Position = [10 height-160 100 20];

function newCalibOK_Callback(hObject, eventdata)

handles = hObject.UserData;

data = read_mixed_csv([handles.introot filesep 'CalibrationDatabase.csv'], ',');
[nparams nrecords] = size(data);
copyfromcolumn = hObject.Parent.Children(2).Value + 1;
newname = hObject.Parent.Children(4).String;

disp(['Copying from ' data{1,copyfromcolumn} ' to ' newname])

data{1,nrecords+1} = hObject.Parent.Children(4).String;

for n = 2:nparams
    data{n,nrecords+1} = data{n,copyfromcolumn};
end

copyfile([handles.introot filesep 'CalibrationDatabase.csv'], [handles.introot filesep 'CalibrationDatabase_backup.csv'])
cell2csv([handles.introot filesep 'CalibrationDatabase.csv'], data, ',')
updateCalibPop(handles.CalibPop, handles)
set(handles.CalibPop, 'Value', nrecords)
CalibBut_Callback(gcbo, eventdata, handles)
close(hObject.Parent)


% --- Executes on button press in imcontrastBtn.
function imcontrastBtn_Callback(hObject, eventdata, handles)
% hObject    handle to imcontrastBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.InterferogramAxes)
imcontrast(gca)


% --- Executes on button press in subtractBtn.
function subtractBtn_Callback(hObject, eventdata, handles)
% hObject    handle to subtractBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
phases = double(csaps2(handles.phase, 1e-10));

axes(handles.PhasediagAxes)
% plot(sum(handles.phase))
% hold all
% plot(sum(phases))
% hold off
surf(handles.phase, 'EdgeColor', 'none')
hold all
surf(phases, 'EdgeColor', 'none')
hold off

handles.phase = handles.phase - phases;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in invertBtn.
function invertBtn_Callback(hObject, eventdata, handles)
% hObject    handle to invertBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = -handles.phase;
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in roiBtn.
function roiBtn_Callback(hObject, eventdata, handles)
% hObject    handle to roiBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.PhaseAxes)
rect = getrect;
x = round(rect(1));
y = round(rect(2));
w = round(rect(3));
h = round(rect(4));

[Ny, Nx] = size(handles.phase);
if (x < 0); x = 1; end
if (x+w > Nx); w = Nx-x; end
if (y < 0); y = 1; end
if (y+h > Ny); h = Ny-y; end

try
    handles.phase = handles.phase(y:y+h, x:x+w);
    handles.mask = handles.mask(y:y+h, x:x+w);
    imagesc(handles.phase); axis image xy;
    guidata(hObject,handles)
catch
    set(handles.StatusBox,'String','Select inside figure'); drawnow
end



function yearBox_Callback(hObject, eventdata, handles)
% hObject    handle to yearBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yearBox as text
%        str2double(get(hObject,'String')) returns contents of yearBox as a double


% --- Executes during object creation, after setting all properties.
function yearBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yearBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in lineoutBtn.
function lineoutBtn_Callback(hObject, eventdata, handles)
% hObject    handle to lineoutBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.DensityAxes)
imagesc(handles.xaxis, handles.yaxis, handles.rho); axis image xy; colorbar;drawnow
pos = ginput(1);
[val xind] = min(abs(handles.xaxis - pos(1)));
[val yind] = min(abs(handles.yaxis - pos(2)));
vline = handles.rho(:,xind);
hline = handles.rho(yind,:);
axes(handles.DensitydiagAxes)
plot(handles.xaxis,hline,handles.yaxis,vline)
legend('Horizontal', 'Vertical')
axes(handles.DensityAxes)
imagesc(handles.xaxis, handles.yaxis, handles.rho); axis image xy; colorbar
line([pos(1) pos(1)], [min(handles.yaxis) max(handles.yaxis)], 'LineStyle', '--', 'Color', 'g')
line([min(handles.xaxis) max(handles.xaxis)],[pos(2) pos(2)], 'LineStyle', '--', 'Color', 'b')


% --- Executes on button press in droboBtn.
function droboBtn_Callback(hObject, eventdata, handles)
% hObject    handle to droboBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MountNAS()


% --- Executes on button press in umountDrobo.
function umountDrobo_Callback(hObject, eventdata, handles)
% hObject    handle to umountDrobo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.StatusBox,'String','Button undefined'); drawnow


% --- Executes on button press in phaseContrastBtn.
function phaseContrastBtn_Callback(hObject, eventdata, handles)
% hObject    handle to phaseContrastBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.PhaseAxes)
imcontrast


% --- Executes on button press in smoothPhaseBtn.
function smoothPhaseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to smoothPhaseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = csaps2(handles.phase, 0.1);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy
guidata(hObject,handles)


% --- Executes on button press in BatchBtn.
function BatchBtn_Callback(hObject, eventdata, handles)
% hObject    handle to BatchBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
startshot = str2num(get(handles.BatchStartBox, 'string'));
endshot = str2num(get(handles.BatchEndBox, 'string'));

shots = startshot:1:endshot;

InterpretDataBatch = struct;
InterpretDataBatch.shots = shots;

for n = 1:length(shots)
    set(handles.ShotBox, 'string', num2str(shots(n)))
    
    FileBut_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
    
    PhaseBut_Callback(hObject, eventdata, handles);
    handles = guidata(hObject);
    
    commandlist = EditableBatchCommands();
    
    for c = 1:length(commandlist)
        switch(commandlist{c})
            case 'Rotate left'
                rLeftBtn_Callback(hObject,eventdata,handles);
            case 'Rotate right'
                rRightBtn_Callback(hObject,eventdata,handles);
            case 'Unwrap left'
                unwrapLeftBtn_Callback(hObject,eventdata,handles);
            case 'Unwrap right'
                unwrapRightBtn_Callback(hObject,eventdata,handles);
            case 'Unwrap top'
                unwrapTopBtn_Callback(hObject,eventdata,handles);
            case 'Unwrap bottom'
                unwrapBottomBtn_Callback(hObject,eventdata,handles);
            case 'Volkov unwrap'
                volkovBtn_Callback(hObject, eventdata, handles);
            case 'Smooth phase'
                smoothPhaseBtn_Callback(hObject, eventdata, handles);
            case 'Subtract gradient'
                subtractBtn_Callback(hObject, eventdata, handles)
        end
        handles = guidata(hObject);
    end
    
    DensityBut_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);
    
    rhobatch{n} = handles.rho;
    phibatch{n} = handles.phase;
    InterpretDataBatch.rho = rhobatch;
    InterpretDataBatch.phi = phibatch;
    assignin('base', 'InterpretDataBatch', InterpretDataBatch)
end


function BatchStartBox_Callback(hObject, eventdata, handles)
% hObject    handle to BatchStartBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BatchStartBox as text
%        str2double(get(hObject,'String')) returns contents of BatchStartBox as a double


% --- Executes during object creation, after setting all properties.
function BatchStartBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BatchStartBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BatchEndBox_Callback(hObject, eventdata, handles)
% hObject    handle to BatchEndBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BatchEndBox as text
%        str2double(get(hObject,'String')) returns contents of BatchEndBox as a double


% --- Executes during object creation, after setting all properties.
function BatchEndBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BatchEndBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in volkovBtn.
function volkovBtn_Callback(hObject, eventdata, handles)
% hObject    handle to volkovBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = VolkovUnwrap(handles);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes on button press in batchEditBtn.
function batchEditBtn_Callback(hObject, eventdata, handles)
% hObject    handle to batchEditBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

edit([handles.introot filesep 'InterpretSource' filesep 'EditableBatchCommands.m'])


% --- Executes on button press in prevBtn.
function prevBtn_Callback(hObject, eventdata, handles)
% hObject    handle to prevBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shot = str2num(get(handles.ShotBox, 'String'));
if(shot > 1)
    set(handles.ShotBox, 'String', num2str(shot - 1))
    FileBut_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in nextBtn.
function nextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shot = str2num(get(handles.ShotBox, 'String'));
if(shot > 0)
    set(handles.ShotBox, 'String', num2str(shot + 1))
    FileBut_Callback(hObject, eventdata, handles);
end


% --- Executes when selected cell(s) is changed in CalibrationTable.
function CalibrationTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to CalibrationTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if(length(eventdata.Indices) > 0)
    rowname = eventdata.Source.RowName{eventdata.Indices(1)};
    if(strcmp(rowname, 'reference'))
        
        experiment = feval(@(x) x{1}{x{2}},get(handles.ExpPop,{'String','Value'}));
        month = str2num(get(handles.MonthBox, 'String'));
        day   = str2num(get(handles.DayBox, 'String'));
        run   = str2num(get(handles.RunBox, 'String'));
        shot  = str2num(get(handles.ShotBox, 'String'));
        year  = str2num(get(handles.yearBox, 'String'));
        data = read_mixed_csv('ExperimentDatabase.csv', ',');
        [nrecords nparams] = size(data);
        
        for n = 2:nrecords
            if(strcmp(experiment, data{n,1}))
                fileexpression = data{n,2};
            end
        end
        
        filepath = GetFilePath(fileexpression, year, month, day, run, shot);
        
        [filename, pathname, filterindex] = uigetfile([filepath filesep '*.*'], 'Click cancel to clear reference');
        
        if(filename == 0)
            button = questdlg('Clear the reference?');
            if(strcmp(button, 'Yes'))
                handles.CalibrationTable.Data{eventdata.Indices(1),1} = 'None';
                guidata(hObject,handles)
                SavecalibBut_Callback(hObject, eventdata, handles)
            end
        else
            handles.CalibrationTable.Data{eventdata.Indices(1),1} = fullfile(pathname, filename);
            guidata(hObject,handles)
            SavecalibBut_Callback(hObject, eventdata, handles)
        end
    end
    
    if(strcmp(rowname, 'x') || strcmp(rowname, 'y') || strcmp(rowname, 'w') || strcmp(rowname, 'h'))
        axes(handles.InterferogramAxes)
        imagesc(imrotate(handles.originaldataimage, handles.calibdata.rotation))
        handles.calibdata
        if(handles.calibdata.x == 0 || ...
           handles.calibdata.y == 0 || ...
           handles.calibdata.w == 0 || ...
           handles.calibdata.h == 0)
            xl = get(gca, 'XLim');
            yl = get(gca, 'YLim');
            xc = mean(xl);
            yc = mean(yl);
            himrect = imrect(handles.InterferogramAxes, [xc/2 yc/2 xc yc]);
        else
            himrect = imrect(handles.InterferogramAxes, [handles.calibdata.x handles.calibdata.y handles.calibdata.w handles.calibdata.h]);
        end
        set(handles.StatusBox,'String','Resize ROI and double click inside when done'); drawnow
        position = wait(himrect);
        delete(himrect)
        
        rotatedimage = imrotate(handles.originaldataimage, handles.calibdata.rotation);
        [Ny, Nx ] = size(rotatedimage);
        x = round(position(1));
        y = round(position(2));
        w = round(position(3));
        h = round(position(4));
        if(x < 1); x = 1; end
        if(x+w>Nx); w = Nx-x; end
        if(y < 1); y = 1; end
        if(y+h>Ny); h = Ny-y; end
        
        handles.calibdata.x = x;
        handles.calibdata.y = y;
        handles.calibdata.w = w;
        handles.calibdata.h = h;
        
        handles.dataimage = rotatedimage(handles.calibdata.y:handles.calibdata.y+handles.calibdata.h, handles.calibdata.x:handles.calibdata.x+handles.calibdata.w);
        imagesc(handles.dataimage)
        axis image xy
        caxis([0 2*mean(mean(handles.dataimage))])
        
        % Update the table (urgh)
        members = fieldnames(handles.calibdata);
        newcalibdata = cell(length(members),1);
        for n = 1:length(members)
            newcalibdata{n} = handles.calibdata.(members{n});
        end
        set(handles.CalibrationTable, 'Data', newcalibdata)
        guidata(hObject,handles)
        SavecalibBut_Callback(hObject, eventdata, handles)
    end
end


% --- Executes on button press in costantiniBtn.
function costantiniBtn_Callback(hObject, eventdata, handles)
% hObject    handle to costantiniBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.phase = cunwrap(handles);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
success = saveGUIState(handles);
if(success)
    delete(hObject);
end


% --- Executes on button press in removeHotPixelsBtn.
function removeHotPixelsBtn_Callback(hObject, eventdata, handles)
% hObject    handle to removeHotPixelsBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for n = 1:10
    handles.phase = RemoveHotPixelButton(handles);
    axes(handles.PhaseAxes)
    imagesc(handles.phase); axis image xy; drawnow
    guidata(hObject,handles)
    UpdateStatus(['Iteration ' num2str(n) ' of 10.'], handles)
end
UpdateStatus('Done', handles)


% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
% hObject    handle to EditMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% This is horrible. TODO: fix
function ReduceFontSize_Callback(hObject, eventdata, handles)

fig = gcf;

for i = 1:length(fig.Children)
    
    if length(fig.Children(i).Children) > 0
        for j = 1:length(fig.Children(i).Children)
            try
                fig.Children(i).Children(j).FontSize = fig.Children(i).Children(j).FontSize/1.5;
            end
        end
    end
    
    try
        fig.Children(i).FontSize = fig.Children(i).FontSize/1.5;
    end
end


% This is horrible. TODO: fix
function IncreaseFontSize_Callback(hObject, eventdata, handles)

fig = gcf;

for i = 1:length(fig.Children)
    
    if length(fig.Children(i).Children) > 0
        for j = 1:length(fig.Children(i).Children)
            try
                fig.Children(i).Children(j).FontSize = fig.Children(i).Children(j).FontSize*1.5;
            end
        end
    end
    
    try
        fig.Children(i).FontSize = fig.Children(i).FontSize*1.5;
    end
end

function Help_Callback(hObject, eventdata, handles)

function OnlineHelp_Callback(hObject, eventdata, handles)
web('https://github.com/jasmcole/Interpret', '-browser')


function DeveloperThesis_Callback(hObject, eventdata, handles)
web('https://spiral.imperial.ac.uk:8443/handle/10044/1/42222', '-browser')


% --- Executes on button press in nonlocalunwrapBtn.
function nonlocalunwrapBtn_Callback(hObject, eventdata, handles)
handles.phase = NonlocalUnwrap(handles);
axes(handles.PhaseAxes)
imagesc(handles.phase); axis image xy;
guidata(hObject,handles)
