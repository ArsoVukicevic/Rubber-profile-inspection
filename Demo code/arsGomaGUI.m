function varargout = arsGomaGUI(varargin)
% ARSGOMAGUI MATLAB code for arsGomaGUI.fig
%      ARSGOMAGUI, by itself, creates a new ARSGOMAGUI or raises the existing
%      singleton*.
%
%      H = ARSGOMAGUI returns the handle to a new ARSGOMAGUI or the handle to
%      the existing singleton*.
%
%      ARSGOMAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARSGOMAGUI.M with the given input arguments.
%
%      ARSGOMAGUI('Property','Value',...) creates a new ARSGOMAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before arsGomaGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to arsGomaGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help arsGomaGUI

% Last Modified by GUIDE v2.5 20-Nov-2018 20:29:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @arsGomaGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @arsGomaGUI_OutputFcn, ...
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

%za resetovanje eventa
function praznaFja(hObject,eventdata)
a=1;
%do nothing fja

function resetClickEvents(hObject,handles)
try
    set(handles.axesKamera.Parent,'WindowButtonDownFcn'    , []);
    set(handles.axesKamera.Parent,'WindowButtonUpFcn'      , []);
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn'  , []);
catch
end
% guidata(hObject, handles);

% --- Executes just before arsGomaGUI is made visible.
function arsGomaGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to arsGomaGUI (see VARARGIN)
initAkviziciju(hObject, eventdata, handles, varargin);


 function initAkviziciju(hObject, eventdata, handles, varargin)
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
global arsData;

dateTime = clock;

arsData.vremeAkvizicije         = [date ' ' num2str(dateTime(4)) 'h ' num2str(dateTime(5)) 'min ' num2str(round(dateTime(6))) 'sec'];

arsData.debugMode               = 0;

arsData.dozvoljeniMaxZaPolozaje = [1.2, 1.3, 0.8, 1.7, 1.75, 1.7, 0.8, 1.3, 1.2];

for i = 1:9
   arsData.skicePozicijaKojeSeIspituju{i}=imread(['sketches of predefined positons\Skica0' mat2str(i) '.jpg']);
end

% Choose default command line output for arsGomaGUI
handles.output = hObject;

%postavi kameru na axes 1; strukture moraju da budu global
try
handles.vid = webcam (1);
axes (handles.axesKamera);
handles.hImage = image (zeros (333,333,3), 'Parent', handles.axesKamera);
preview (handles.vid, handles.hImage);
catch
arsData.imgPath     ='_Test images\STR1#1.JPG';
arsData.img         = imread(arsData.imgPath);
for i=1:9
    arsData.rezultatiMerenja{i}.imgVelika = [];
    arsData.rezultatiMerenja{i}.imgMala   = [];
end

axes(handles.axesKamera);
%#1 Korak - nadji krugove i izvrsi kalibraciju
[arsData.trenutniProfil.koordinateCentar arsData.trenutniProfil.xPravac arsData.trenutniProfil.yPravac arsData.trenutniProfil.scaleFaktor] = skripta1_nadjiKrugKalibrisiKooSistem(arsData.imgPath, arsData.debugMode);

%#2 Korak - Nadji konturu profila
[arsData.trenutniProfil.konturaProfila   arsData.trenutniProfil.imgFiltrirano                                          ] = skripta2_nadjiKontureProfila(arsData.trenutniProfil.koordinateCentar,...
                                                                                                                                        arsData.trenutniProfil.xPravac,.......
                                                                                                                                        arsData.trenutniProfil.yPravac,.......
                                                                                                                                        arsData.trenutniProfil.scaleFaktor,...
                                                                                                                                        arsData.imgPath,...
                                                                                                                                        arsData.debugMode);
%#3 Korak - Preklopi crtez preko nosacha
[arsData.postavljenaKontura              arsData.rezGTnoveKoordinate   ] = skripta3_PreklopiDetektovanuSaoriginalnomKonturom(arsData.trenutniProfil.koordinateCentar,....
                                                                                                             arsData.trenutniProfil.xPravac,.......
                                                                                                             arsData.trenutniProfil.yPravac,.......
                                                                                                             arsData.trenutniProfil.scaleFaktor,...
                                                                                                             arsData.imgPath,..........
                                                                                                             'orgKonturaMM',... 
                                                                                                             arsData.debugMode);
%#4 Korak - Odradi Euclide distance i nadji Skeletone
[arsData.trenutniProfil.euclidianDistanceImg....
 arsData.trenutniProfil.skeletonImg.............
 arsData.trenutniProfil.koordinateBifurkacija...
 arsData.trenutniProfil.koordianateCentralnihLinija] = skripta4_OdradiEuclidianDistanceNadjiSkeletoneKontura_b(arsData.trenutniProfil.konturaProfila, arsData.imgPath, arsData.trenutniProfil.scaleFaktor, arsData.debugMode);
end

popupMenuListaPolozaja_Callback(hObject, eventdata, handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes arsGomaGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = arsGomaGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnSacuvajReport.
function btnSacuvajReport_Callback(hObject, eventdata, handles)
% hObject    handle to btnSacuvajReport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sacuvajRezultateMerenja();

% --- Executes on selection change in popupMenuListaPolozaja.
function popupMenuListaPolozaja_Callback(hObject, eventdata, handles)
% hObject    handle to popupMenuListaPolozaja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global arsData;
arsData.idPozicijeKojaSeIspituje = handles.popupMenuListaPolozaja.Value;
iscrtaj(hObject, eventdata, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupMenuListaPolozaja contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupMenuListaPolozaja


% --- Executes during object creation, after setting all properties.
function popupMenuListaPolozaja_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupMenuListaPolozaja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnMoveUp.
function btnMoveUp_Callback(hObject, eventdata, handles)
% hObject    handle to btnMoveUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMoveLeft.
function btnMoveLeft_Callback(hObject, eventdata, handles)
% hObject    handle to btnMoveLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMoveRight.
function btnMoveRight_Callback(hObject, eventdata, handles)
% hObject    handle to btnMoveRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnMoveDown.
function btnMoveDown_Callback(hObject, eventdata, handles)
% hObject    handle to btnMoveDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnRotateLeft.
function btnRotateLeft_Callback(hObject, eventdata, handles)
% hObject    handle to btnRotateLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnRotateRight.
function btnRotateRight_Callback(hObject, eventdata, handles)
% hObject    handle to btnRotateRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txtPomerajTranslacija_Callback(hObject, eventdata, handles)
% hObject    handle to txtPomerajTranslacija (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPomerajTranslacija as text
%        str2double(get(hObject,'String')) returns contents of txtPomerajTranslacija as a double


% --- Executes during object creation, after setting all properties.
function txtPomerajTranslacija_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPomerajTranslacija (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPomerajRotacija_Callback(hObject, eventdata, handles)
% hObject    handle to txtPomerajRotacija (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPomerajRotacija as text
%        str2double(get(hObject,'String')) returns contents of txtPomerajRotacija as a double


% --- Executes during object creation, after setting all properties.
function txtPomerajRotacija_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPomerajRotacija (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function iscrtaj(hObject, eventdata, handles, rez)
global arsData;
global arsDataFake;
%koji se polozaj ispituje
axes(handles.axisSkica);
imshow(arsData.skicePozicijaKojeSeIspituju{arsData.idPozicijeKojaSeIspituje});

handles.textTrenutniPolozaj.String = ['Polozaj ' num2str(handles.popupMenuListaPolozaja.Value)];
%--
%#5 Uradi alignement za inspekciju pojedinih delova
if ~exist('rez')
    rez = skripta5_OdradiInspekcijuPojedinihDelova(arsData.trenutniProfil.koordinateBifurkacija(:,[2 1 3]),.........
                                         arsData.trenutniProfil.koordianateCentralnihLinija(:,[2 1 3]),...
                                         arsData.trenutniProfil.konturaProfila(:,[1 2 3]),................ NOVI PROFIL 
                                         arsData.rezGTnoveKoordinate,... referentni profil
                                         arsData.imgPath,.....
                                         arsData.idPozicijeKojaSeIspituje,....
                                         handles.axesKamera, false);
end                                      
global dragableLinije;
global centarRotacije;
dragableLinije = [];
axes(handles.axesKamera);
cla(handles.axesKamera,'reset')

if handles.checkboxKoristiEuclidianDistance.Value
    img           = arsData.trenutniProfil.euclidianDistanceImg;  
    img(:,:,2)    = arsData.trenutniProfil.euclidianDistanceImg;
    img(:,:,3)    = arsData.trenutniProfil.euclidianDistanceImg;
    GrayIndex = uint8(floor(arsData.trenutniProfil.euclidianDistanceImg*255/arsData.dozvoljeniMaxZaPolozaje(handles.popupMenuListaPolozaja.Value)));
    Map       = jet(255);
    RGB       = ind2rgb(GrayIndex, Map);
    RGB(isnan(img))=255;
    imshow(RGB);
else
    if handles.checkboxTransparentnaPozadina.Value
        img                 = arsData.trenutniProfil.euclidianDistanceImg; 
        img(~isnan(img))    = 100; 
        img(isnan(img))     = 222;%procenat transparentnosti
        img=double(img);
        img = img/100.0;
        slikaZaPrikazivanje = mat2gray(imread(arsData.imgPath));
        slikaZaPrikazivanje = slikaZaPrikazivanje(:,:,1);
        slikaZaPrikazivanje = slikaZaPrikazivanje .* img;
        imshow(slikaZaPrikazivanje);
    else
        imshow(arsData.imgPath);
    end
end
centarRotacije = mean(rez.TackaKojaSeTrazi{:});
dragableLinije{1}=line(rez.Dicp(:,1)   , rez.Dicp(:,2), 'color','blue');  hold on;
scatter(rez.tackaBifurkacija(:,1), rez.tackaBifurkacija(:,2));
% scatter(rez.pointCloud(:,1), rez.pointCloud(:,2));
for i = 1 : numel(rez.tolerancija)
    dragableLinije{1+i}=line(rez.tolerancija{i}(:,1), rez.tolerancija{i}(:,2), 'color','red');
end


axes(handles.axesScreenShot);
cla(handles.axesScreenShot,'reset')
if ~isempty(arsData.rezultatiMerenja{handles.popupMenuListaPolozaja.Value}.imgMala)
	imshow(arsData.rezultatiMerenja{handles.popupMenuListaPolozaja.Value}.imgMala);
end
        
txt = 'Uradjeni delovi:';
for i=1:9
    if ~isempty(arsData.rezultatiMerenja{i}.imgMala)
        imshow(arsData.rezultatiMerenja{handles.popupMenuListaPolozaja.Value}.imgMala);
        txt = [txt num2str(i) '  '];
    end
end
handles.textUradjeniPogledi.String = txt;
txt = ['"'];
arsDataFake = arsData;


% --------------------------------------------------------------------
function uitoggletoolPanCrtez_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolPanCrtez1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'pointer','hand');
%% dragable linija
    set(handles.axesKamera.Parent,'WindowButtonDownFcn', @startDragFcn);
    set(handles.axesKamera.Parent,'WindowButtonUpFcn'  , @endDragFcn);
%%    
    function startDragFcn(hObject,eventdata)
    handles = guidata(hObject);  
    global pocetnaPozicijaMisa;
    set(gcf,'pointer','fleur');
    pocetnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', @draggingFcn);
    global linija;
    linija = line(handles.axesKamera, [pocetnaPozicijaMisa(1,1), pocetnaPozicijaMisa(1,1)+10],....
                   [pocetnaPozicijaMisa(1,2), pocetnaPozicijaMisa(1,2)+10], 'color','red');
%%
function endDragFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    a=get(handles.axesKamera, 'CurrentPoint');
    set(gcf,'pointer','hand');  
    global linija;
    if ~isempty(linija)
        delete(linija);
        clear('linija');
    end
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', '');
%%
function draggingFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    global linija;
    global pocetnaPozicijaMisa;
    trenutnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');       
    linija.XData(end)=trenutnaPozicijaMisa(1,1);
    linija.YData(end)=trenutnaPozicijaMisa(1,2);
    global dragableLinije;
    vektorPomeraja = trenutnaPozicijaMisa-pocetnaPozicijaMisa;
    pocetnaPozicijaMisa = trenutnaPozicijaMisa;
    for i = 1:numel(dragableLinije)
        dragableLinije{i}.XData=dragableLinije{i}.XData+vektorPomeraja(1,1);
        dragableLinije{i}.YData=dragableLinije{i}.YData+vektorPomeraja(1,2);
    end
    global centarRotacije;
    centarRotacije = centarRotacije + vektorPomeraja(1,:);
    
    


% --------------------------------------------------------------------
function uitoolbar1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uitoolbar1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  ROTACIJA START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function uitoggletoolRotateCrtez_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolRotateCrtez (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% dragable linija
    set(handles.axesKamera.Parent,'WindowButtonDownFcn', @startRotateFcn);
    set(handles.axesKamera.Parent,'WindowButtonUpFcn'  , @endRotateFcn);
  
function startRotateFcn(hObject,eventdata)
    handles = guidata(hObject);  
    global pocetnaPozicijaMisa;
    global dragableLinije;
    global centarRotacije;
    set(gcf,'pointer','fleur');
    pocetnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', @rotatingFcn);
    global linija;
    linija = line(handles.axesKamera, [pocetnaPozicijaMisa(1,1), centarRotacije(1,1)],....
                                      [centarRotacije(1,2)     , centarRotacije(1,2)], 'color','red');

function endRotateFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    a=get(handles.axesKamera, 'CurrentPoint');
    set(gcf,'pointer','hand');  
    global linija;
    if ~isempty(linija)
        delete(linija);
        clear('linija');
    end
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', '');
%%
function rotatingFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    global linija;
    global pocetnaPozicijaMisa;
    global centarRotacije;
    trenutnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');   
    linija.XData(3)=trenutnaPozicijaMisa(1,1);
    linija.YData(3)=trenutnaPozicijaMisa(1,2);
    global dragableLinije;
    pom  = pocetnaPozicijaMisa(1,:)-centarRotacije(1,:);
    pom1 = trenutnaPozicijaMisa(1,:)-centarRotacije(1,:);
    pom2 = [    0.7071   -0.7071         0;....
                0.7071    0.7071         0;....
                     0         0    1.0000]*pom1'; 
    pom2 = pom2';       
    ugaoRotacijeStepeni = AngioIvusMath.getUgaoIzmedjuDvaVektora(pom,pom1);
    znak = -(-1)^ AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB([0 0 0],pom, pom1, pom2);
    matricaRotacije = [cos(deg2rad(ugaoRotacijeStepeni)*znak), -sin(deg2rad(ugaoRotacijeStepeni)*znak);...
                       sin(deg2rad(ugaoRotacijeStepeni)*znak),  cos(deg2rad(ugaoRotacijeStepeni)*znak)];
    pocetnaPozicijaMisa = trenutnaPozicijaMisa;               
    for i = 1:numel(dragableLinije)
        XY = [dragableLinije{i}.XData-centarRotacije(1,1);...
              dragableLinije{i}.YData-centarRotacije(1,2)];
        XY = matricaRotacije * XY;
        dragableLinije{i}.XData=XY(1,:)+centarRotacije(1,1);
        dragableLinije{i}.YData=XY(2,:)+centarRotacije(1,2);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  ROTACIJA ENDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
               
               


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

buttonId = double(get(gcf,'CurrentCharacter'));

if isempty(buttonId)
    return;
end

switch buttonId    
   case 115
        global arsData;
        imgCela = getimage(handles.axesKamera); %napravi img od cele slike        
        frame   = getframe(handles.axesKamera);
        imgMala = frame2im(frame);
        axes(handles.axesScreenShot)
        imshow(imgMala);
        arsData.rezultatiMerenja{handles.popupMenuListaPolozaja.Value}.imgVelika = imgCela;
        arsData.rezultatiMerenja{handles.popupMenuListaPolozaja.Value}.imgMala   = imgMala;
   case 49 
       handles.popupMenuListaPolozaja.Value = 1; ...Taster 1  Pozicija br 1-Levo spoljasnje perce        
   case 50 
       handles.popupMenuListaPolozaja.Value = 2; ...Taster 2
   case 51 
       handles.popupMenuListaPolozaja.Value = 3;    
   case 52 
       handles.popupMenuListaPolozaja.Value = 4;
   case 53 
       handles.popupMenuListaPolozaja.Value = 5;
   case 54 
       handles.popupMenuListaPolozaja.Value = 6;
   case 55 
       handles.popupMenuListaPolozaja.Value = 7;
   case 56 
       handles.popupMenuListaPolozaja.Value = 8;
   case 57 
       handles.popupMenuListaPolozaja.Value = 9;    
   case 114
       ...Rotacija
       resetClickEvents(hObject, handles);   
       uitoggletoolRotateCrtez_ClickedCallback(hObject, eventdata, handles);
   case 116
       resetClickEvents(hObject, handles) ; 
       ...T Translacija 
       uitoggletoolPanCrtez_ClickedCallback(hObject, eventdata, handles);
   case 121        
%        zoom on; ...Y Zoom 
   case 100
%         pan on; ...T Pad  
   case 101
        ...E Merenje-kotiranje  
       uitoggletoolMerenjeRastojanja_ClickedCallback(hObject, eventdata, handles);
   case 97
        ...A Automatsko poravnanje
        resetClickEvents(hObject, handles);
        uitoggletoolLocalFitting_ClickedCallback(hObject,eventdata,handles);
   otherwise
      ...fprintf('Invalid grade\n' );
      aaaa=1;  
end ...switch

if(buttonId>48 && buttonId<58) 
       popupMenuListaPolozaja_Callback(hObject, eventdata, handles);
end


% --- Executes on button press in checkboxKoristiEuclidianDistance.
function checkboxKoristiEuclidianDistance_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxKoristiEuclidianDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxKoristiEuclidianDistance


function sacuvajRezultateMerenja()
global arsData;
folderName = ['rezultati/' arsData.vremeAkvizicije];
if ~exist(folderName, 'dir')
	mkdir(folderName)
end
for i = 1:numel(arsData.rezultatiMerenja)
    try
        imwrite(arsData.rezultatiMerenja{i}.imgVelika, [folderName '/imgVelika.jpg']);
        imwrite(arsData.rezultatiMerenja{i}.imgMala,   [folderName '/imgMala_Pozicija#' num2str(i) '.jpg']);
    catch
    end
end


% --- Executes on button press in checkboxTransparentnaPozadina.
function checkboxTransparentnaPozadina_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTransparentnaPozadina (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTransparentnaPozadina


%----------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------
function uitoggletoolMerenjeRastojanja_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolMerenjeRastojanja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% dragable linija
    set(handles.axesKamera.Parent,'WindowButtonDownFcn', @startMerenjeFcn);
    set(handles.axesKamera.Parent,'WindowButtonUpFcn'  , @endRMerenjeFcn);
% --------------------------------------------------------------------
function uitoggletoolMerenjeRastojanja_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolMerenjeRastojanja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
resetClickEvents(handles);
    
function startMerenjeFcn(hObject,eventdata)
    handles = guidata(hObject);  
    global pocetnaPozicijaMisa;
    global dragableLinije;
    set(gcf,'pointer','fleur');
    pocetnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', @merenjeFcn);
    global linija; global textRastojanje;
    linija = line(handles.axesKamera, [pocetnaPozicijaMisa(1,1), pocetnaPozicijaMisa(1,1)+0.1],...................
                                      [pocetnaPozicijaMisa(1,2), pocetnaPozicijaMisa(1,2)+0.1],'color','red');....  
    textRastojanje       = text(pocetnaPozicijaMisa(1,1),pocetnaPozicijaMisa(1,2),'0','color','red');............                                  

function endRMerenjeFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    a=get(handles.axesKamera, 'CurrentPoint');
    set(gcf,'pointer','hand');  
    global linija;
    if ~isempty(linija)
%         delete(linija);
        clear('linija');
    end
    set(handles.axesKamera.Parent,'WindowButtonMotionFcn', '');
%%
function merenjeFcn(hObject,eventdata,handles)
    handles = guidata(hObject);
    global linija; global textRastojanje; global arsData;
    global pocetnaPozicijaMisa;
    trenutnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint');   
    linija.XData(2)      = trenutnaPozicijaMisa(1,1);
    linija.YData(2)      = trenutnaPozicijaMisa(1,2);
    textRastojanje.Position = [mean(linija.XData),mean(linija.YData),0];
    textRastojanje.String   = num2str(arsData.trenutniProfil.scaleFaktor*norm([linija.XData(1), linija.YData(1)]-[linija.XData(2), linija.YData(2)]));
    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  ROTACIJA ENDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%vrsi lokalno poravnanje sto najblizih tacki referentne i trenutne konture
%u odnosu na trenutni polozaj misa
% --------------------------------------------------------------------
function uitoggletoolLocalFitting_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolLocalFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
set(handles.axesKamera.Parent,'WindowButtonDownFcn'    , @uitoggletoolLocalFittingButtonDown);


function uitoggletoolLocalFittingButtonDown(hObject,eventdata,handles)
global arsData;
global dragableLinije;
global centarRotacije;
handles = guidata(hObject);
trenutnaPozicijaMisa = get(handles.axesKamera, 'CurrentPoint'); 
trenutnaPozicijaMisa = [trenutnaPozicijaMisa(1,1:2) 0];
%nadji najblizu tacku GT bifurkacija trenuntoj poziciji misa


rez = skripta5_OdradiInspekcijuPojedinihDelova(arsData.trenutniProfil.koordinateBifurkacija(:,[2 1 3]),.........
                                         arsData.trenutniProfil.koordianateCentralnihLinija(:,[2 1 3]),...
                                         arsData.trenutniProfil.konturaProfila(:,[1 2 3]),................ NOVI PROFIL 
                                         arsData.rezGTnoveKoordinate,... referentni profil
                                         arsData.imgPath,.....
                                         arsData.idPozicijeKojaSeIspituje,....
                                         handles.axesKamera, false,....
                                         trenutnaPozicijaMisa);
                                     
% arsDataFake.trenutniProfil.konturaProfila             =  rez.Dicp       ;
% arsDataFake.rezGTnoveKoordinate.tolerancije                =  rez.tolerancija;
% arsDataFake.rezGTnoveKoordinate.koordinateCentralnihLinija =  rez.a;
% arsDataFake.rezGTnoveKoordinate.TackeZaMatching            =  rez.kontura;

% % centarRotacije = mean(rez.TackaKojaSeTrazi{:});
% % dragableLinije{1}.XData = rez.Dicp(:,1);
% % dragableLinije{1}.YData = rez.Dicp(:,2);
% % for i = 1 : numel(rez.tolerancija)
% %     dragableLinije{1+i}.XData = rez.tolerancija{i}(:,1);
% %     dragableLinije{1+i}.YData = rez.tolerancija{i}(:,2);
% % end
% iscrtaj(hObject, eventda1ta, handles, rez);
