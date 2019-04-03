%Skripta vrsi kalibraciju - vraca centar koordinatnog pocetka i pravce X Y

%INPUTS
    %imgPath - putanja ka slici
%OUTPUTS
    %koordinateCentar - centar u kome se postavlja kontura za preklapanje
    %X                - pravac od najveceg kruga ka srednjem
    %Y                - pravac od najveceg kruga ka najmanjem
    %scaleFaktor      - koliko mm ima jedan pixel
function [koordinateCentar xPravac yPravac scaleFaktor] = skripta1_nadjiKrugKalibrisiKooSistem(imgPath, dubugMode)
if nargin ==0
    imgPath   = 'arsIMG2.JPG';
    dubugMode = 0;
elseif nargin == 1
    dubugMode = 0;
end
%UCITAJ SLIKU
img = imread(imgPath);
R   = img(:,:,1) > 200;
G   = img(:,:,2) < 100;
img = R .* G; 
% img = gradient(img);
img=mat2gray(~img);
if dubugMode
    imshow(img);
end
%nadji krugove
radiusRange     = [44  55];
[centers1,radii1] = imfindcircles(img,[50  66], 'ObjectPolarity','dark', 'Method','TwoStage',  'Sensitivity',0.97);
[centers2,radii2] = imfindcircles(img,[26  33], 'ObjectPolarity','dark', 'Method','TwoStage',  'Sensitivity',0.98);
[centers3,radii3] = imfindcircles(img,[70  88], 'ObjectPolarity','dark', 'Method','TwoStage',  'Sensitivity',0.97);
radii   = [radii1;radii2;radii3];
centers = [centers1;centers2;centers3];

for i = numel(radii):-1:1
    if ~proveriDaLiJeKrugValidan(img, centers(i,:), radii(i) )
        radii(i)    =[];  
        centers(i,:)=[];
    end
end
% dubugMode=1;
if dubugMode
    imshow(img);
%     imshow(imgPath);
    for i = 1:numel(radii)
        iscrtajKrug(centers(i,1), centers(i,2), radii(i));
	end
end

[radii id] = sort(radii);
centers    = centers(id,:);
%nadji centar
radiusCentar      = radii(3,:)  ;
koordinateCentar  = round(centers(3,:)); koordinateCentar(3)=0;
%nadji x koordinatu
xPravac           = centers(1,:)-centers(3,:);
xPravac           = xPravac/norm(xPravac)    ;  xPravac(3)=0;
%nadji x koordinatu
yPravac           = centers(3,:)-centers(2,:);
yPravac           = yPravac/norm(yPravac)    ;  yPravac(3)=0;
%scale faktor       %br mm  / brpixela
% rastojanje izmedju centra najveceg i najmanjeg je 5mm, izmedju najveceg srednjeg 6mm
scaleFaktor       =  mean([5/norm(centers(3,:)-centers(1,:)), 6/norm(centers(3,:)-centers(2,:))]);
% scaleFaktor       = mean([1 2 3] ./ radii');
end
function iscrtajKrug(x,y,r)
    th    = 0:pi/50:2*pi   ;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plotLine([xunit(:), yunit(:)],2);
end
%INPUTS
    %centar  - koordinate centra kruga
    %radijs  - POLUprecnik kruga
%OUTPUTS
    %validan - true/false
% proverava se da li su svi pikseli kvadratica unutar kruga 
% (lakse je za selektovanje nego da pravim konture kruga pa da selektujem pixele)
% ako su skoro svi pixeli unutar kvadracitca crni to je validan krug
function validan = proveriDaLiJeKrugValidan(img, centar, radijus)
    %
%     imshow(img);
%     iscrtajKrug(centar(1),centar(2),radijus);
    a = round(radijus/sqrt(2));
    b = round(radijus*0.2);
    pixeli1 =  img(centar(2)-a:centar(2)+a, centar(1)-a:centar(1)+a);
    %% takocje svi pixeli na coskovima moraju biti beli
    pixeli2 = img(centar(2)-(radijus+b):centar(2)+(radijus+b), centar(1)-(radijus+b):centar(1)+(radijus+b));
    if ((mean(pixeli1(:)) > 0.05) || (mean(pixeli2(:)) < 0.25))
        validan = false;
    else
        validan = true;
    end
end