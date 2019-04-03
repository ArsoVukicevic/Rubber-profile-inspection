%Skripta nalazi konturu profila na osnovu input img
%INPUTS
    %imgPath - slika sa fotoaparata
%OUTPUTS
    %kontura       - 2D kontura
    %imgFiltrirano - filtrirana sika bez sumova
function [kontura, imgFiltrirano] = skripta2_nadjiKontureProfila(koordinateCentar, xPravac, yPravac, scaleFaktor, imgPath, debugMode)
if ~exist('debugMode')
    debugMode = 0;
end
img = imread(imgPath);
img = mat2gray(img(:,:,1));

%nadji konture
img(img>0.2)=1;
img(img<=0.2)=0;
%odseci ivice
% % img(1:10,:)=1; img(end-10:end,:)=1;
% % img(:,1:10)=1; img(:,end-10:end)=1;
%obrisi krugove
%dimenzije nalepnice su 10x10mm
% % img(round(koordinateCentar(2)-6/scaleFaktor):round(koordinateCentar(2)+4/scaleFaktor),:)=1;

%% OBRIShI ShUMOVE
img = ~img;
BW = bwareaopen(img, 1000);
if debugMode
    figure; imshowpair(img,BW,'montage');
end

%% NADJI IVICE
img = mat2gray(BW);
c=contourc(img,[.5 .5]);
c=c'; c(:,3)=0;  
c = c(2:end-1,:);

if debugMode
    figure; imshow(img); hold on; plotLine(c);
end
%OUTPUTS
kontura       =   c;
imgFiltrirano = img;
end