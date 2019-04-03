%Skripta detektuje centralne linije i 
%INPUTS
    %detektovanaKontura - detektovane ivice profila
    %img                - img, sluzi da se zna dimenzija slike
    %debugMode          - da li treba stampati medjurezultate ili ne
%OUTPTUS
    %skeletonImg        - centralne linije
    %geodesicImg        - matrica rastojanja pixela od centralnih linija
function [geodesicImg skeletonImg] = skripta4_OdradiEuclidianDistanceNadjiSkeletoneKontura(detektovanaKontura, img, debugMode)
    if isstr(img)
        img = imread(img);
        img = img(:,:,1);
    end
    %#1 Napravi poly mask od konture
    [x y]=size(img);
    imgPolyMask = poly2mask(detektovanaKontura(:,1), detektovanaKontura(:,2), x, y);
    imgPolyMask = int16(~imgPolyMask);
    if debugMode
        figure;
        imshow(mat2gray(imgPolyMask));
    end
    
    %#2 Nadji Centralne Linije
    skeletonImg = bwmorph(~imgPolyMask,'skel',inf);    
    centralneLinijePixelsID  = find(skeletonImg==1);
    
    % Prebaci iz ID u koordinate na slici
    [rez] = arsIMG.pxelIDS2Koordinate(img,centralneLinijePixelsID);
    
    %#3 Nadji geodesic distance         
    [geodesicImg] = bwdistgeodesic(logical(~imgPolyMask), logical(skeletonImg)); 
    geodesicImg(isnan(geodesicImg))=0;
    %Vizuelizacija
    if debugMode
        figure, imshow(mat2gray(geodesicImg)), title('Distance transform of ~bw');
        figure; imshow(skeletonImg)
    end
end