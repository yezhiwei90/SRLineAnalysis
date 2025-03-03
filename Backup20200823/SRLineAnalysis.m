%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024

function varargout = SRLineAnalysis(varargin)
% SRLINEANALYSIS MATLAB code for SRLineAnalysis.fig
%      SRLINEANALYSIS, by itself, creates a new SRLINEANALYSIS or raises the existing
%      singleton*.
%
%      H = SRLINEANALYSIS returns the handle to a new SRLINEANALYSIS or the handle to
%      the existing singleton*.
%
%      SRLINEANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SRLINEANALYSIS.M with the given input arguments.
%
%      SRLINEANALYSIS('Property','Value',...) creates a new SRLINEANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SRLineAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SRLineAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SRLineAnalysis

% Last Modified by GUIDE v2.5 23-Aug-2020 18:16:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SRLineAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @SRLineAnalysis_OutputFcn, ...
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


% --- Executes just before SRLineAnalysis is made visible.
function SRLineAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SRLineAnalysis (see VARARGIN)

% Choose default command line output for SRLineAnalysis
handles.output = hObject;
pathname = 'D:\Figure\LocalizationMethod';
handles.pathname = pathname;
filename = {};
setappdata(handles.SRLineAnalysis, 'Filename_F',filename);
setappdata(handles.SRLineAnalysis, 'Filename_SR',filename);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SRLineAnalysis wait for user response (see UIRESUME)
% uiwait(handles.SRLineAnalysis);


% --- Outputs from this function are returned to the command line.
function varargout = SRLineAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Adress1.
function Adress1_Callback(hObject, eventdata, handles)
% hObject    handle to Adress1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = handles.pathname;
[FileName,PathName] = uigetfile('*.csv','Select the SR line data', pathname);
handles.pathname = PathName;
setappdata(handles.SRLineAnalysis, 'Filename_SR',[PathName FileName]);
guidata(hObject, handles);

% --- Executes on button press in Adress2.
function Adress2_Callback(hObject, eventdata, handles)
% hObject    handle to Adress2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = handles.pathname;
[FileName,PathName] = uigetfile('*.csv','Select the wide field line data', pathname);
handles.pathname = PathName;
setappdata(handles.SRLineAnalysis, 'Filename_F',[PathName FileName]);
guidata(hObject, handles);


function magnification_Callback(hObject, eventdata, handles)
% hObject    handle to magnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magnification as text
%        str2double(get(hObject,'String')) returns contents of magnification as a double

% --- Executes during object creation, after setting all properties.
function magnification_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magnification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mag = str2double(handles.magnification.String);
adrsr = getappdata(handles.SRLineAnalysis,'Filename_SR');
adrf = getappdata(handles.SRLineAnalysis,'Filename_F');
pixsize = str2double(handles.PixSizeEdit.String);
gaussian_chk = handles.GaussianCheck.Value;
srcounts_chk = handles.SRCountsCheck.Value;
gaussian_num = str2double(handles.GaussianNumEdit.String);
plotsize = str2range(handles.PlotSizeEdit.String);
plotstyle = handles.Stylepopupmenu.Value;

if gaussian_chk
    R = LineAnalysis(adrsr,adrf,mag);
    ylim([0 1.1]);
    xlim([0 10000]);
    %R.figure.CurrentAxes.XTick = [0 300 600];
    pathname = handles.pathname;
    [FileName,PathName] = uiputfile('*.tif','Save image as...', [pathname 'SRLinePlot.tif']);
    xlabel('Distance (nm)');
    ylabel('Normalized Intensity');
    mgraphics(R.figure,plotsize);
    R.figure = ModifyPaperSizes(R.figure);
    %mSignificantDigits(R.figure.CurrentAxes);
    export_fig(R.figure,[PathName, FileName],'-png','-m4');
    saveas(R.figure, [PathName, FileName(1:end-4)],'pdf');
    fid=fopen([PathName, FileName(1:end-4),'_AnalysisResults.csv'],'w+');
    num = numel(R.width_WF);
    for i = 1 : num
        fprintf(fid,'Wide field FWHM (Normally its unit was um): %.2d\r\nWF fit Center Position (Normally its unit was um): %.2d\r\n',R.width_WF(i),R.center_WF(i));
        fprintf(fid, 'Super resolution FWHM (Normally its unit was um): %.2d\r\nSR fit Center Position (Normally its unit was um): %.2d\r\n\r\n', R.width_SR(i), R.center_SR(i));
    end
    fclose(fid);
else
    data_sr = imagejread(adrsr); [sr_x,sr_y] = datareader(data_sr);
    if ~srcounts_chk
        sr_y = sr_y/max(sr_y);
    end
    data_wf = imagejread(adrf); [wf_x,wf_y] = datareader(data_wf); wf_y = wf_y/max(wf_y);
    [ fit_type,fit_option ] = GaussianFitInitializer( max(wf_x), gaussian_num );
    [fw,goffw] = fit(wf_x,wf_y,fit_type, fit_option);
    disp(goffw.adjrsquare); width_WF = zeros(gaussian_num,1);
    for i = 1 : gaussian_num
       eval(['width_WF(i) = fw.c' num2str(i) ' * 2 * sqrt(log(2));']);
    end
    xfit=(1:1000)-1;
    wf_xfit = xfit/999 * max(wf_x);
    wf_yfit = fw(wf_xfit);
    wf_xfit = wf_xfit*100 * pixsize;
    [ bin_x, bin_y ] = Binning_WideFieldData(wf_x, wf_y, mag);
    bin_x = bin_x*100  * pixsize; sr_x = sr_x.*1000; % converting to nanometers;
    if ~srcounts_chk
        h = figure; ax = axes; hold on;
        plot(ax, bin_x, bin_y, 'x', 'color', Colorslib('grey')/255,'MarkerSize',2);
        bar(ax, sr_x, sr_y, 'FaceColor', Colorslib('red')./255, ...
            'EdgeColor',Colorslib('black')./255);
        plot(ax, wf_xfit, wf_yfit, '-', 'color', Colorslib('grey')/255);
        xlabel('Distance (nm)');
        ylabel('Normalized Intensity');
        legend('WF','SR','WF_fit');
        h = mgraphics(h);
        ylim([0 1.0]);
        xlim([0 max(sr_x)]);
    else
        h = figure; ax1 = axes; hold on; ax2 = axes; hold on;
        ax2.Color = 'none';ax2.YAxisLocation = 'right';
        plot(ax1, bin_x, bin_y, 'x', 'color', Colorslib('grey')/255,'MarkerSize',2);
        bar(ax2, sr_x, sr_y, 'FaceColor', Colorslib('red')./255, ...
            'EdgeColor',Colorslib('black')./255);
        plot(ax1, wf_xfit, wf_yfit, '-', 'color', Colorslib('grey')/255);
        legend(ax1, 'WF','WF_fit','Location', 'northeast');legend(ax2,'SR','Location','northwest');
        ax1.XLim = [0 max(sr_x)];
        ax2.XLim = [0 max(sr_x)];ax2.XTick = [];ax2.XColor = 'none';
        pos = ax1.Position;
        ax2.Position(1:2) = pos(1:2);
        xlabel(ax1, 'Distance (nm)');
        ylabel(ax1,'Normalized Intensity');
        ylabel(ax2,'Average Counts');
        mSignificantDigits(ax1);
        h = mgraphics(h, plotsize);
    end
pathname = handles.pathname;
[FileName,PathName] = uiputfile('*.png','Save image as...', [pathname 'SRLinePlot']);
export_fig(h,[PathName, FileName],'-png','-m4');
fid=fopen([PathName, FileName,'_AnalysisResults.csv'],'w+');
for i = 1 : gaussian_num
    fprintf(fid,'Wide field FWHM: %.2d\r\n',width_WF(1));
end
fclose(fid);
end


function [x,y] = datareader(data)
descript = data.header;
num = numel(descript);
for i = 1 : num
    descriptitem = descript{i};
    descriptitem = removenumalphabetelemtns(descriptitem);
    switch descriptitem
        case 'X'
            x = data.data(:,i);
        case 'Y'
            y = data.data(:,i);
        otherwise
            
    end
end

function str = removenumalphabetelemtns(str)
    n = numel(str);
    chk = true(1,n);
    for i = 1:n
        chking = regexp(str(i),'\w','match');
        if isempty(chking)
            chk(i) = false;
        end
    end
    str = str(chk);

% --- Executes on button press in GaussianCheck.
function GaussianCheck_Callback(hObject, eventdata, handles)
% hObject    handle to GaussianCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of GaussianCheck



function GaussianNumEdit_Callback(hObject, eventdata, handles)
% hObject    handle to GaussianNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GaussianNumEdit as text
%        str2double(get(hObject,'String')) returns contents of GaussianNumEdit as a double


% --- Executes during object creation, after setting all properties.
function GaussianNumEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GaussianNumEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PixSizeEdit.
function PixSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PixSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function PixSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PixSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in SRCountsCheck.
function SRCountsCheck_Callback(hObject, eventdata, handles)
% hObject    handle to SRCountsCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SRCountsCheck


% --- Executes on selection change in Stylepopupmenu.
function Stylepopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Stylepopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Stylepopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Stylepopupmenu


% --- Executes during object creation, after setting all properties.
function Stylepopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stylepopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PlotSizeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PlotSizeEdit as text
%        str2double(get(hObject,'String')) returns contents of PlotSizeEdit as a double


% --- Executes during object creation, after setting all properties.
function PlotSizeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotSizeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
