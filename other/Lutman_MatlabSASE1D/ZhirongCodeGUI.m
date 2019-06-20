function varargout = ZhirongCodeGUI(varargin)
% ZHIRONGCODEGUI MATLAB code for ZhirongCodeGUI.fig
%      ZHIRONGCODEGUI, by itself, creates a new ZHIRONGCODEGUI or raises the existing
%      singleton*.
%
%      H = ZHIRONGCODEGUI returns the handle to a new ZHIRONGCODEGUI or the handle to
%      the existing singleton*.
%
%      ZHIRONGCODEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZHIRONGCODEGUI.M with the given input arguments.
%
%      ZHIRONGCODEGUI('Property','Value',...) creates a new ZHIRONGCODEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ZhirongCodeGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ZhirongCodeGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ZhirongCodeGUI

% Last Modified by GUIDE v2.5 05-Mar-2014 10:21:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ZhirongCodeGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ZhirongCodeGUI_OutputFcn, ...
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


% --- Executes just before ZhirongCodeGUI is made visible.
function ZhirongCodeGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ZhirongCodeGUI (see VARARGIN)

T1_data{1,1}=2;T1_data{2,1}=6000;T1_data{3,1}=1024;T1_data{4,1}=0;
set(handles.T1,'data',T1_data)

T2_data{1,1}=5.810595384217249e+03;T2_data{2,1}=0.5e-4;T2_data{3,1}= 0.4e-6;T2_data{4,1}=3000.0;
T2_data{5,1}=25.0; T2_data{6,1}=50; T2_data{7,1}=0;T2_data{8,1}=0;
set(handles.T2,'data',T2_data)

T3_data{1,1}=0.03;T3_data{2,1}=40;T3_data{3,1}=1;T3_data{4,1}=3.5;
set(handles.T3,'data',T3_data)

T4_data{1,1}=16; T4_data{2,1}=200;T4_data{3,1}=40;T4_data{4,1}='[0,0,4]';
T4_data{5,1}=100; T4_data{6,1}=56; T4_data{7,1}=0;
set(handles.T4,'data',T4_data)

T5_data{1,1}=0; T5_data{2,1}=1;T5_data{3,1}=1230;T5_data{4,1}=30;
set(handles.T5,'data',T5_data)

T6_data{1,1}=1; T6_data{2,1}=0;T6_data{3,1}=0;
set(handles.T6,'data',T6_data)

T7_data{1,1}=0; T7_data{2,1}=0;T7_data{3,1}=1;T7_data{4,1}=0;
T7_data{5,1}=0;T7_data{6,1}=0;T7_data{7,1}=0;
set(handles.T7,'data',T7_data)

T8_data{1,1}=1; T8_data{2,1}=1;T8_data{3,1}=1;T8_data{4,1}='Sim';
set(handles.T8,'data',T8_data)

set(handles.ManualSetupPanel,'visible','off')
handles.ColorIDLE=get(handles.CT,'BackgroundColor');
handles.ColorON=[0,1,0];
handles.ColorRED=[1,0,0];
handles.ColorWAIT=[1,1,0];

handles.BranchablesDivision=2;


% Choose default command line output for ZhirongCodeGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
ST_Callback(hObject, eventdata, handles)
% UIWAIT makes ZhirongCodeGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ZhirongCodeGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ManualK.
function ManualK_Callback(hObject, eventdata, handles)
set(handles.ManualSetupPanel,'visible','on')
for II=1:10
   varname=eval(['handles.uipanel',num2str(II)]);
   set(varname,'visible','off') 
end
str=get(handles.text20,'String')
if(strcmp(str,'No Manual Configuration Stored'))    
    FillTable(handles,1)
else
    Current=get(handles.text20,'Userdata');
    set(handles.ManualConf,'data',Current);
end
% --- Executes on selection change in EBunchShape.
function EBunchShape_Callback(hObject, eventdata, handles)
% hObject    handle to EBunchShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns EBunchShape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from EBunchShape


% --- Executes during object creation, after setting all properties.
function EBunchShape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EBunchShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EBunchName_Callback(hObject, eventdata, handles)
% hObject    handle to EBunchName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EBunchName as text
%        str2double(get(hObject,'String')) returns contents of EBunchName as a double


% --- Executes during object creation, after setting all properties.
function EBunchName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EBunchName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CrystalIN.
function CrystalIN_Callback(hObject, eventdata, handles)
if(strcmp(get(hObject,'String'),'Crystal Inserted'))
    set(hObject,'String','No Crystal Inserted')
else
    set(hObject,'String','Crystal Inserted')
end


% --- Executes on button press in ManualD.
function ManualD_Callback(hObject, eventdata, handles)
set(handles.ManualSetupPanel,'visible','on')
for II=1:10
   varname=eval(['handles.uipanel',num2str(II)]);
   set(varname,'visible','off') 
end
FillTable(handles,1)

% --- Executes on button press in ManualB.
function ManualB_Callback(hObject, eventdata, handles)
set(handles.ManualSetupPanel,'visible','on')
for II=1:10
   varname=eval(['handles.uipanel',num2str(II)]);
   set(varname,'visible','off') 
end
FillTable(handles,1)

% --- Executes on selection change in SeedShape.
function SeedShape_Callback(hObject, eventdata, handles)
% hObject    handle to SeedShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SeedShape contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SeedShape


% --- Executes during object creation, after setting all properties.
function SeedShape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SeedShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SAVESETUP.
function SAVESETUP_Callback(hObject, eventdata, handles)
[FILENAME, PATHNAME]=uiputfile('*.mat','give filename');
ZhirongCodeConf.T1=get(handles.T1,'data');
ZhirongCodeConf.T2=get(handles.T2,'data');
ZhirongCodeConf.T3=get(handles.T3,'data');
ZhirongCodeConf.T4=get(handles.T4,'data');
ZhirongCodeConf.T5=get(handles.T5,'data');
ZhirongCodeConf.T6=get(handles.T6,'data');
ZhirongCodeConf.T7=get(handles.T7,'data');
ZhirongCodeConf.T8=get(handles.T8,'data');
ZhirongCodeConf.SeedShape=get(handles.SeedShape,'value');
ZhirongCodeConf.CrystalInString=get(handles.CrystalIN,'String');
ZhirongCodeConf.BunchShape=get(handles.EBunchShape,'value');
ZhirongCodeConf.NameShape=get(handles.EBunchName,'string');
ZhirongCodeConf.UndulatorConfSet=get(handles.text20,'string');
ZhirongCodeConf.Scan=get(handles.Scan,'value');
ZhirongCodeConf.ScanOption=get(handles.ScanOption,'value');
ZhirongCodeConf.FROM=get(handles.ScanFROM,'string');
ZhirongCodeConf.TO=get(handles.ScanTO,'string');
ZhirongCodeConf.STEPS=get(handles.ScanSTEPS,'string');
ZhirongCodeConf.ENERGY_JITTER=get(handles.E_Jitter,'string');
ZhirongCodeConf.SCANAFTERMODULE=get(handles.edit13,'string');
ZhirongCodeConf.UndulatorManualConfiguration=get(handles.ManualConf,'data');
ZhirongCodeConf.UndulatorManualConfigurationRows=get(handles.ManualConf,'rowname');
save([PATHNAME,FILENAME],'ZhirongCodeConf');




% --- Executes on selection change in Scan.
function Scan_Callback(hObject, eventdata, handles)
% hObject    handle to Scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Scan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Scan


% --- Executes during object creation, after setting all properties.
function Scan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ScanOption.
function ScanOption_Callback(hObject, eventdata, handles)
% hObject    handle to ScanOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ScanOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ScanOption


% --- Executes during object creation, after setting all properties.
function ScanOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScanOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScanFROM_Callback(hObject, eventdata, handles)
% hObject    handle to ScanFROM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScanFROM as text
%        str2double(get(hObject,'String')) returns contents of ScanFROM as a double


% --- Executes during object creation, after setting all properties.
function ScanFROM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScanFROM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScanTO_Callback(hObject, eventdata, handles)
% hObject    handle to ScanTO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScanTO as text
%        str2double(get(hObject,'String')) returns contents of ScanTO as a double


% --- Executes during object creation, after setting all properties.
function ScanTO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScanTO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function E_Jitter_Callback(hObject, eventdata, handles)
% hObject    handle to E_Jitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Jitter as text
%        str2double(get(hObject,'String')) returns contents of E_Jitter as a double


% --- Executes during object creation, after setting all properties.
function E_Jitter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Jitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScanSTEPS_Callback(hObject, eventdata, handles)
% hObject    handle to ScanSTEPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScanSTEPS as text
%        str2double(get(hObject,'String')) returns contents of ScanSTEPS as a double


% --- Executes during object creation, after setting all properties.
function ScanSTEPS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScanSTEPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
T1=get(handles.T1,'data');
T2=get(handles.T2,'data');
T3=get(handles.T3,'data');
T4=get(handles.T4,'data');
T5=get(handles.T5,'data');
T6=get(handles.T6,'data');
T7=get(handles.T7,'data');
T8=get(handles.T8,'data');
SeedShape=get(handles.SeedShape,'value');
CrystalInString=get(handles.CrystalIN,'String');
BunchShape=get(handles.EBunchShape,'value');
NameShape=get(handles.EBunchName,'string');
ManualUndConf=get(handles.text20,'string');
Scan=get(handles.Scan,'value');
ScanOption=get(handles.ScanOption,'value');
FROM=get(handles.ScanFROM,'string');
TO=get(handles.ScanTO,'string');
STEPS=get(handles.ScanSTEPS,'string');
ENERGY_JITTER=get(handles.E_Jitter,'string');
EJ=str2num(ENERGY_JITTER);
SCANAFTERMODULE=get(handles.edit13,'string');
UndulatorManualConfiguration=get(handles.ManualConf,'data');
UndulatorManualConfigurationRows=get(handles.ManualConf,'rowname');

if((Scan>=handles.BranchablesDivision) && (ScanOption==2)) 
    BRANCHINGCYCLE=1;
else
    BRANCHINGCYCLE=0;
end

InitializationScript
UndulatorConfigurationScript
ENERGIA=energy;
ElectronBunchConfigurationScript
StartingFieldConfigurationScript

%save TEMP

if BRANCHINGCYCLE
   for FileID=1:T8{2}    
       for StartID=1:T8{1} %Simulate Up to branch
           
           for BranchID=1:STEPS %Restart from Branch after doing something
                RunID=(StartID-1)*STEPS+BranchID;
       
           end
       end      
   end
else
    for FileID=1:T8{2}  
       for RunID=1:T8{1} %No branching every run is a new run
            if(Scan==2)
                energyV=linspace(FROM,TO,STEPS);
                energy=energyV(RunID);
                ElectronBunchConfigurationScript
            end
            if(EJ>0)
                energy=energy + randn(1)*EJ*ENERGIA;
                ElectronBunchConfigurationScript
            end
            if(Scan==3) %scan undulator k
            
            end
            if(Scan==4) %delay (phase shift)
            
            end
            if(Scan==5) %taper strenght
            
            end
            if(Scan==6) %taper start
            
            end
            power_z=[];
            power_z_smart=[];
            full_a2=[];
            cumulative_z=0;
            ENERGIA=T2{1};
            
            ElectronBunchConfigurationScript
            StartingFieldConfigurationScript
            
            for ElID=1:numel(Element);
               switch(Element{ElID}.type)
                   case 'U' %FEL evolution in an undulator
                       %GammaCValues=sqrt(unduPeriod/2/radWavelength*(1+Element{ElID}.Kvalue.^2/2))
                       
                       GammaFake=sqrt( unduPeriod /radWavelength/2*(1+T3{4}^2/2) );
                       MevFake=GammaFake*mc2;
                       MevEq=(1+Element{ElID}.Kvalue.^2/2)/(1+T3{4}^2/2)*MevFake;
                       Egain=MevEq-energy;
                       Egain=-Egain;
                       UndulatorLength=Element{ElID}.length;
                       cumulative_z=cumulative_z+UndulatorLength;
                       disp(['Undulator section m= ',num2str(UndulatorLength)]);
                       FEL_Code_Script;
                       a2 = ar.^2+ai.^2;
                       full_a2=[full_a2,ar(:,2:end) + 1i*ai(:,2:end)];
                       power_z_smart=[power_z_smart,mean(a2(find(shape>0.05),:))*rhoPbeam];
                       power_z=[power_z,mean(a2)*rhoPbeam];
                   case 'K' %destroy bunching
                       for k = 1:s_steps-1
                            [thet0,gam0] = load_bucket(npart,mean(CurrentPhaseSpace.theta(k,:)),std(CurrentPhaseSpace.theta(k,:)),5-T1{4},Ns);	% load each bucket
                            CurrentPhaseSpace.theta(k,:)=thet0;
                            CurrentPhaseSpace.gamma(k,:)=gam0;
                        end
                   case 'D' %apply delay to field
                       disp(['apply delay ',num2str(Element{ElID}.length)]);
                       Final_field=Final_field*exp(1i*Element{ElID}.length*2*pi);
                   case 'C' %do self seeding crystal stuff
               end 
            end
            if(T6{3}), T6{1}=1;, end %If need to separate 2 colors, need to calculate spectrum
            % Evaluate Spectrum and store as output
            [SA,SB]=size(full_a2);
            OUTPUT{RunID}.time=s_step_safe/c*(1:SA);
            OUTPUT{RunID}.full_field=full_a2*sqrt(rhoPbeam);
            OUTPUT{RunID}.power_z_smart=power_z_smart;
            OUTPUT{RunID}.power_z=power_z;
            OUTPUT{RunID}.z=linspace(0,cumulative_z,length(OUTPUT{RunID}.power_z));
            %save TEMP3s
            if(T6{1})
                fieldFFT_z=zeros(SA+2*T6{2},SB);
                for SpecID=1:SB
                    fieldFFT_z(:,SpecID)=s_step_safe*circshift(fft([zeros(T6{2},1);OUTPUT{RunID}.full_field(:,SpecID);zeros(T6{2},1)]),round(length(Final_field)/2+T6{2}));
                end
                OrelMin=unduPeriod*T1{1}*(-1/2);
                OrelMax=unduPeriod*T1{1}*(+1/2);
                OUTPUT{RunID}.orel=linspace(OrelMin,OrelMax,SA+2*T6{2});
                OUTPUT{RunID}.fieldFFT_z=fieldFFT_z;
            end
            % Evaluate Two Colors (when it will be added... in the future)
            if(T6{3})
                OUTPUT{RunID}.Color1=0*OUTPUT{RunID}.full_field; OUTPUT{RunID}.Color2 = OUTPUT{RunID}.Color1;
                for SpecID=1:SB
                    TempColor1=zeros(size(OUTPUT{RunID}.fieldFFT_z(:,SpecID)));
                    TempColor2=zeros(size(OUTPUT{RunID}.fieldFFT_z(:,SpecID)));
                    TempColor1(1:(SA+T6{2}))=OUTPUT{RunID}.fieldFFT_z(1:(SA+T6{2}),SpecID);
                    TempColor2((SA+T6{2}+1):end)=OUTPUT{RunID}.fieldFFT_z((1+SA+T6{2}):end,SpecID);
                    TempColor1=1/(s_step_safe)*circshift(ifft(TempColor1),-round(length(Final_field)/2+T6{2}));
                    TempColor2=1/(s_step_safe)*circshift(ifft(TempColor2),-round(length(Final_field)/2+T6{2}));
                    OUTPUT{RunID}.Color1(:,SpecID)=TempColor1((T6{2}+1):(T6{2}+SA));
                    OUTPUT{RunID}.Color2(:,SpecID)=TempColor2((T6{2}+1):(T6{2}+SA));
                end 
            end
            % Store in output
            OUTPUT{RunID}.Energy=ENERGIA;
       end  
       eval(['save ',T8{4},num2str(FileID),' OUTPUT -v7.3'])
    end
end

%save TEMP2s


% --- Executes on button press in LoadSetup.
function LoadSetup_Callback(hObject, eventdata, handles)
[FILENAME, PATHNAME]=uigetfile('*.mat','give filename');
load([PATHNAME,FILENAME]);
set(handles.T1,'data',ZhirongCodeConf.T1);
set(handles.T2,'data',ZhirongCodeConf.T2);
set(handles.T3,'data',ZhirongCodeConf.T3);
set(handles.T4,'data',ZhirongCodeConf.T4);
set(handles.T5,'data',ZhirongCodeConf.T5);
set(handles.T6,'data',ZhirongCodeConf.T6);
set(handles.T7,'data',ZhirongCodeConf.T7);
set(handles.T8,'data',ZhirongCodeConf.T8);
set(handles.SeedShape,'value',ZhirongCodeConf.SeedShape);
set(handles.CrystalIN,'String',ZhirongCodeConf.CrystalInString);
set(handles.EBunchShape,'value',ZhirongCodeConf.BunchShape);
set(handles.EBunchName,'string',ZhirongCodeConf.NameShape);
set(handles.text20,'string',ZhirongCodeConf.UndulatorConfSet);
set(handles.Scan,'value',ZhirongCodeConf.Scan);
set(handles.ScanOption,'value',ZhirongCodeConf.ScanOption);
set(handles.ScanFROM,'string',ZhirongCodeConf.FROM);
set(handles.ScanTO,'string',ZhirongCodeConf.TO);
set(handles.ScanSTEPS,'string',ZhirongCodeConf.STEPS);
set(handles.E_Jitter,'string',ZhirongCodeConf.ENERGY_JITTER);
set(handles.edit13,'string',ZhirongCodeConf.SCANAFTERMODULE);
set(handles.ManualConf,'data',ZhirongCodeConf.UndulatorManualConfiguration);
set(handles.ManualConf,'rowname',ZhirongCodeConf.UndulatorManualConfigurationRows);


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ApplyAndClose.
function ApplyAndClose_Callback(hObject, eventdata, handles)
set(handles.text20,'String','Manual Configuration Stored')
Current=get(handles.ManualConf,'data');
set(handles.text20,'Userdata',Current);
set(handles.ManualSetupPanel,'visible','off')
for II=1:10
   varname=eval(['handles.uipanel',num2str(II)]);
   set(varname,'visible','on') 
end

% --- Executes on button press in DiscardandClose.
function DiscardandClose_Callback(hObject, eventdata, handles)
set(handles.ManualSetupPanel,'visible','off')
for II=1:10
   varname=eval(['handles.uipanel',num2str(II)]);
   set(varname,'visible','on') 
end

function FillTable(handles,FROM)
T4=get(handles.T4,'data');
T3=get(handles.T3,'data')
if(strcmp(get(handles.CrystalIN,'String'),'Crystal Inserted'))
    Crystal=1;
else
    Crystal=0;
end
RowName={};
inserted=0;
for II=1:T3{3} 
    inserted=inserted+1;
    RowName{inserted}=['Und #',num2str(II)];
    data{inserted,1}=T3{4};
    data{inserted,2}=0;
    data{inserted,3}=false;
    data{inserted,4}=T3{2};
    data{inserted,5}=true;
    data{inserted,6}='';
    if Crystal
        if(II==T4{1})
            RowName{inserted}=['Crystal'];
            inserted=inserted+1;
        end
    end
end

set(handles.ManualConf,'RowName',RowName)
set(handles.ManualConf,'data',data)


% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
C=get(handles.CT,'backgroundcolor');
V1=(get(handles.edit17,'String'));
V2=(get(handles.edit18,'String'));
V3=(get(handles.edit19,'String'));
V4=(get(handles.edit20,'String'));
V5=str2num(get(handles.edit21,'String'));
if(sum(C==handles.ColorON)==3) %Cont Taper
    Current=get(handles.ManualConf,'data');
    Current{V5,6}=['[',V1,',',V2,',',V3,',',V4,']'];
    Current{V5,1}=NaN;
else %Discrete Taper
    V1=str2num(V1);
    V2=str2num(V2);
    V3=str2num(V3);
    V4=str2num(V4);    
    inserted=0;
    quadratic=0;
    Current=get(handles.ManualConf,'data');
    ItemNames=get(handles.ManualConf,'rowname');
    for II=1:numel(ItemNames)
       if(ItemNames{II}(1)=='U')
           if(Current{II,5}) 
              inserted=inserted+1;
              if(II>=V3)
                  quadratic=quadratic+1;
              end
              Kvalue= V1 + (inserted-1)*V2 + quadratic*(quadratic+1)/2*V4;
              Current{II,1}=Kvalue;
              Current{II,6}='';
           end
       end
    end
end
set(handles.ManualConf,'data',Current);

% --- Executes during object creation, after setting all properties.
function ManualD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ManualD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in ContTaper.
function ContTaper_Callback(hObject, eventdata, handles)
% hObject    handle to ContTaper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in CT.
function CT_Callback(hObject, eventdata, handles)
set(handles.CT,'BackgroundColor',handles.ColorON);
set(handles.ST,'BackgroundColor',handles.ColorIDLE);
set(handles.text22,'String','Linear Taper [/m]')
set(handles.text24,'String','Quadratic Start [m]')
set(handles.text23,'String','Quadratic Taper [/m]')
set(handles.UndSelect,'visible','off')
set(handles.edit21,'visible','on')
set(handles.text27,'visible','on')


% --- Executes on button press in ST.
function ST_Callback(hObject, eventdata, handles)
set(handles.ST,'BackgroundColor',handles.ColorON);
set(handles.CT,'BackgroundColor',handles.ColorIDLE);
set(handles.text22,'String','Linear Taper /U')
set(handles.text24,'String','Quadratic Start U#')
set(handles.text23,'String','Quadratic Taper /U')
set(handles.UndSelect,'visible','on')
set(handles.edit21,'visible','off')
set(handles.text27,'visible','off')


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');
Delta=str2num(get(handles.edit14,'String'));
for II=1:numel(ItemNames)
 %  if(ItemNames{II}(1)=='U')
        if(sum(find(Delta==II)))
            Current{II,5}=true;
        else
            Current{II,5}=false;
        end
  % end
end
set(handles.ManualConf,'data',Current);


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');
Delta=str2num(get(handles.edit15,'String'))
counter=0;
putting=true;
for II=1:numel(ItemNames)
   if(ItemNames{II}(1)=='U')
        counter=counter+1;
        Current{II,5}=putting;
        if(counter==Delta)
            putting=~putting;
            counter=0;
        end
   end
end
set(handles.ManualConf,'data',Current);

% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');

for II=1:numel(ItemNames)
   if(ItemNames{II}(1)=='U')
       if(~Current{II,5})
            Current{II,5}=true;
       end
   end
end
set(handles.ManualConf,'data',Current);

% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');

for II=1:numel(ItemNames)
   if(ItemNames{II}(1)=='U')
       if(Current{II,5})
            Current{II,5}=false;
       end
   end
end
set(handles.ManualConf,'data',Current);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');
Delta=str2num(get(handles.edit16,'String'));

for II=1:numel(ItemNames)
   if(ItemNames{II}(1)=='U')
       if(Current{II,5})
            Current{II,1}=Current{II,1}+Delta;
       end
   end
end
set(handles.ManualConf,'data',Current);



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
set(handles.text20,'String','No Manual Configuration Stored');



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
Current=get(handles.ManualConf,'data');
ItemNames=get(handles.ManualConf,'rowname');
for II=1:numel(ItemNames)
   if(ItemNames{II}(1)=='U')
        if(Current{II,5})
            Current{II,5}=false
        else
            Current{II,5}=true
        end
   end
end
set(handles.ManualConf,'data',Current);



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in T5.
function T5_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to T5 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
c   = 2.99792458E8;
e   = 1.60217733E-19;
h   = 4.13566751691*10^-15; % [eV s]
Current=get(hObject,'data');
if(eventdata.Indices(1)==2)
    Current{3}=h*c/(eventdata.NewData*10^-9);
    set(handles.T5,'data',Current);
end

if(eventdata.Indices(1)==3)
    Current{2}=h*c/eventdata.NewData;
    set(handles.T5,'data',Current);
end


% --- Executes when entered data in editable cell(s) in ManualConf.
function ManualConf_CellEditCallback(hObject, eventdata, handles)
