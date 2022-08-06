function varargout = selectSaccades2(varargin)
% SELECTSACCADES2 MATLAB code for selectSaccades2.fig
%      SELECTSACCADES2, by itself, creates a new SELECTSACCADES2 or raises the existing
%      singleton*.
%
%      H = SELECTSACCADES2 returns the handle to a new SELECTSACCADES2 or the handle to
%      the existing singleton*.
%
%      SELECTSACCADES2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTSACCADES2.M with the given input arguments.
%
%      SELECTSACCADES2('Property','Value',...) creates a new SELECTSACCADES2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectSaccades2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectSaccades2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectSaccades2

% Last Modified by GUIDE v2.5 09-Mar-2019 14:26:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectSaccades2_OpeningFcn, ...
                   'gui_OutputFcn',  @selectSaccades2_OutputFcn, ...
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


% --- Executes just before selectSaccades2 is made visible.
function selectSaccades2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectSaccades2 (see VARARGIN)

% get varargin
Saccade_timing = varargin{1};
Eye_per_trial = varargin{2};
str = varargin{3};

handles.fr = 200; % camera frame rate Hz
handles.azMs = 800; % length of eye trace used for each saccade in ms
numnan = 2; % maximum number of NaN values in the 800ms trace accepted
trunc = 0.96875; % ratio of existing samples to the total samples within 800ms accepted

for i=1:size(Saccade_timing,1)
    sacTim = Saccade_timing{i,1};
    if ~isempty(sacTim)
        EYE = Eye_per_trial{i,1};
        for j=1:size(sacTim,1)
%             az = EYE(1,EYE(4,:)>sacTim(j,2)-500 & EYE(4,:)<sacTim(j,2)+500); 
%             az = EYE(1,EYE(4,:)>sacTim(j,1)-400 & EYE(4,:)<sacTim(j,1)+400);
            az = EYE(1,EYE(4,:)>sacTim(j,1)-1000 & EYE(4,:)<sacTim(j,1)+1000);
            if sum(isnan(az))>numnan %any(isnan(az))
            else
                if length(az)<handles.fr*handles.azMs/1000*trunc
                else
%                         plot(EYE(4,EYE(4,:)>sacTim(j,2)-500 & EYE(4,:)<sacTim(j,2)+500)-sacTim(j,2),az) 
%                         plot(sacTim(j,1)-sacTim(j,2),EYE(1,EYE(4,:)==sacTim(j,1)),'og','MarkerFaceColor','g','MarkerSize',4)
%                         plot(0,EYE(1,EYE(4,:)==sacTim(j,2)),'or','MarkerFaceColor','r','MarkerSize',4)
                    plot(handles.rawax, EYE(4,EYE(4,:)>sacTim(j,1)-1000 & EYE(4,:)<sacTim(j,1)+1000)-sacTim(j,1),az)
%                         plot(0,EYE(1,EYE(4,:)==sacTim(j,1)),'og','MarkerFaceColor','g','MarkerSize',4)
%                         plot(sacTim(j,2)-sacTim(j,1),EYE(1,EYE(4,:)==sacTim(j,2)),'or','MarkerFaceColor','r','MarkerSize',4)
                    hold(handles.rawax, 'on');
                    plot(handles.alignedax, EYE(4,EYE(4,:)>sacTim(j,1)-1000 & EYE(4,:)<sacTim(j,1)+1000)-sacTim(j,1),az-EYE(1,EYE(4,:)==sacTim(j,1)))
%                     plot(0,0,'og','MarkerFaceColor','g','MarkerSize',4)
%                     plot(sacTim(j,2)-sacTim(j,1),EYE(1,EYE(4,:)==sacTim(j,2))-EYE(1,EYE(4,:)==sacTim(j,1)),'or','MarkerFaceColor','r','MarkerSize',4)
                    hold(handles.alignedax, 'on');
                end
            end
        end
    end
end

h=title(handles.rawax, str);
set(h,'interpreter','none');
hold(handles.rawax, 'off');
hold(handles.alignedax, 'off');

handles.Saccade_timing = Saccade_timing;
handles.Eye_per_trial = Eye_per_trial;
handles.str=str;

% Choose default command line output for selectSaccades2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes selectSaccades2 wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = selectSaccades2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Saccade_timing;



function minaz_Callback(hObject, eventdata, handles)
% hObject    handle to minaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minaz as text
%        str2double(get(hObject,'String')) returns contents of minaz as a double

minaz = str2double(get(hObject,'String'));
if isnan(minaz) || ~isreal(minaz)  
    % isdouble returns NaN for non-numbers and f1 cannot be complex
    % Disable the Plot button and change its string to say why
    set(handles.changeButton,'String','Value invalid')
    set(handles.changeButton,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Plot button with its original name
    set(handles.changeButton,'String','Change')
    set(handles.changeButton,'Enable','on')
end


% --- Executes during object creation, after setting all properties.
function minaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxaz_Callback(hObject, eventdata, handles)
% hObject    handle to maxaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxaz as text
%        str2double(get(hObject,'String')) returns contents of maxaz as a double

maxaz = str2double(get(hObject,'String'));
if isnan(maxaz) || ~isreal(maxaz)  
    % isdouble returns NaN for non-numbers and f1 cannot be complex
    % Disable the Plot button and change its string to say why
    set(handles.changeButton,'String','Value invalid')
    set(handles.changeButton,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Plot button with its original name
    set(handles.changeButton,'String','Change')
    set(handles.changeButton,'Enable','on')
end


% --- Executes during object creation, after setting all properties.
function maxaz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxaz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in changeButton.
function changeButton_Callback(hObject, eventdata, handles)
% hObject    handle to changeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% First read in the new parameter values
minaz = str2double(get(handles.minaz,'String'));
maxaz = str2double(get(handles.maxaz,'String'));
numnan = str2double(get(handles.numnan,'String'));
trunc = str2double(get(handles.trunc,'String'));

% clear figures
% clf(handles.rawax,'reset');
% clf(handles.alignedax,'reset');

Saccade_timing = handles.Saccade_timing;
Eye_per_trial = handles.Eye_per_trial;

for i=1:size(Saccade_timing,1)
    sacTim = Saccade_timing{i,1};
    if ~isempty(sacTim)
        EYE = Eye_per_trial{i,1};
        for j=1:size(sacTim,1)
%             az = EYE(1,EYE(4,:)>sacTim(j,2)-500 & EYE(4,:)<sacTim(j,2)+500); 
            az = EYE(1,EYE(4,:)>sacTim(j,1)-400 & EYE(4,:)<sacTim(j,1)+400);
            if sum(isnan(az))>numnan %any(isnan(az))
            else
                if min(az)<minaz || max(az)>maxaz
                else
                    if length(az)<handles.fr*handles.azMs/1000*trunc
                    else
%                         plot(EYE(4,EYE(4,:)>sacTim(j,2)-500 & EYE(4,:)<sacTim(j,2)+500)-sacTim(j,2),az) 
%                         plot(sacTim(j,1)-sacTim(j,2),EYE(1,EYE(4,:)==sacTim(j,1)),'og','MarkerFaceColor','g','MarkerSize',4)
%                         plot(0,EYE(1,EYE(4,:)==sacTim(j,2)),'or','MarkerFaceColor','r','MarkerSize',4)
                        plot(handles.rawax, EYE(4,EYE(4,:)>sacTim(j,1)-400 & EYE(4,:)<sacTim(j,1)+400)-sacTim(j,1),az)
%                         plot(0,EYE(1,EYE(4,:)==sacTim(j,1)),'og','MarkerFaceColor','g','MarkerSize',4)
%                         plot(sacTim(j,2)-sacTim(j,1),EYE(1,EYE(4,:)==sacTim(j,2)),'or','MarkerFaceColor','r','MarkerSize',4)
                        hold(handles.rawax, 'on');
                        plot(handles.alignedax, EYE(4,EYE(4,:)>sacTim(j,1)-400 & EYE(4,:)<sacTim(j,1)+400)-sacTim(j,1),az-EYE(1,EYE(4,:)==sacTim(j,1)))
    %                     plot(0,0,'og','MarkerFaceColor','g','MarkerSize',4)
    %                     plot(sacTim(j,2)-sacTim(j,1),EYE(1,EYE(4,:)==sacTim(j,2))-EYE(1,EYE(4,:)==sacTim(j,1)),'or','MarkerFaceColor','r','MarkerSize',4)
                        hold(handles.alignedax, 'on');
                    end
                end
            end
        end
    end
end

h=title(handles.rawax, handles.str);
set(h,'interpreter','none');
hold(handles.rawax, 'off')
hold(handles.alignedax, 'off')



% --- Executes on button press in confirmButton.
function confirmButton_Callback(hObject, eventdata, handles)
% hObject    handle to confirmButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% First read in the new parameter values
minaz = str2double(get(handles.minaz,'String'));
maxaz = str2double(get(handles.maxaz,'String'));
numnan = str2double(get(handles.numnan,'String'));
trunc = str2double(get(handles.trunc,'String'));

Saccade_timing = handles.Saccade_timing;
Eye_per_trial = handles.Eye_per_trial;

for i=1:size(Saccade_timing,1)
    sacTim = Saccade_timing{i,1};
    directions = Saccade_timing{i,2};
    amplitudes = Saccade_timing{i,3};
    positions = Saccade_timing{i,4};
    eyeTrace = Saccade_timing{i,5};
    if ~isempty(sacTim)
        EYE = Eye_per_trial{i,1};
        del = [];
        for j=1:size(sacTim,1)
%             az = EYE(1,EYE(4,:)>sacTim(j,2)-500 & EYE(4,:)<sacTim(j,2)+500); 
            az = EYE(1,EYE(4,:)>sacTim(j,1)-400 & EYE(4,:)<sacTim(j,1)+400);
            if sum(isnan(az))>numnan %any(isnan(az))
                del = [del j];
            else
                if min(az)<minaz || max(az)>maxaz
                    del = [del j];
                else
                    if length(az)<handles.fr*handles.azMs/1000*trunc
                        del = [del j];
                    else
                    end
                end
            end
        end
        sacTim(del,:)=[]; % delete the flagged saccades
        directions(del,:)=[];
        amplitudes(del,:)=[];
        positions(del,:)=[];
        eyeTrace(del,:)=[];
        Saccade_timing{i,1} = sacTim; % update the timing
        Saccade_timing{i,2} = directions; % update the timing
        Saccade_timing{i,3} = amplitudes; % update the timing
        Saccade_timing{i,4} = positions;
        Saccade_timing{i,5} = eyeTrace;
    end
end

handles.Saccade_timing = Saccade_timing; % update handles

% Update handles structure
guidata(hObject, handles);

% Resume
hFig = ancestor(hObject,'Figure');
if isequal(get(hFig, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT
    uiresume(handles.figure1);
end




function numnan_Callback(hObject, eventdata, handles)
% hObject    handle to numnan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numnan as text
%        str2double(get(hObject,'String')) returns contents of numnan as a double

numnan = str2double(get(hObject,'String'));
if isnan(numnan) || ~isreal(numnan)  
    % isdouble returns NaN for non-numbers and f1 cannot be complex
    % Disable the Plot button and change its string to say why
    set(handles.changeButton,'String','Value invalid')
    set(handles.changeButton,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Plot button with its original name
    set(handles.changeButton,'String','Change')
    set(handles.changeButton,'Enable','on')
end


% --- Executes during object creation, after setting all properties.
function numnan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numnan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trunc_Callback(hObject, eventdata, handles)
% hObject    handle to trunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trunc as text
%        str2double(get(hObject,'String')) returns contents of trunc as a double

trunc = str2double(get(hObject,'String'));
if isnan(trunc) || ~isreal(trunc)  
    % isdouble returns NaN for non-numbers and f1 cannot be complex
    % Disable the Plot button and change its string to say why
    set(handles.changeButton,'String','Value invalid')
    set(handles.changeButton,'Enable','off')
    % Give the edit text box focus so user can correct the error
    uicontrol(hObject)
else 
    % Enable the Plot button with its original name
    set(handles.changeButton,'String','Change')
    set(handles.changeButton,'Enable','on')
end


% --- Executes during object creation, after setting all properties.
function trunc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trunc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Units','points');
set(handles.figure1,'PaperUnits','points');
size = get(handles.figure1,'Position');

size = size(3:4);
set(handles.figure1,'PaperSize',size);
set(handles.figure1,'PaperPosition',[0,0,size(1),size(2)]);

exp_id = handles.str;
folder = strcat('\',exp_id,'_ephys');
fname = strcat('\',exp_id,'_selSaccades');
% print(handles.figure1, strcat('D:\Satoru\Documents\MATLAB\Ephys data',folder,fname),'-depsc');
% print(handles.figure1, strcat('E:\MATLAB2\Ephys data',folder,fname),'-depsc');
print(handles.figure1, strcat('F:\MATLAB3\Ephys data',folder,fname),'-depsc');
