classdef arsIMG
    properties(SetAccess = public, GetAccess = public)
        points;
    end%properties
    methods(Static=true)
        %% 
        function [rez] = nadjiUnutrasnjuSpoljanjuKonturuOdBinarneSlikeIvica(bwImg, komponenteKonture)
            debugMode = 0;
            bwImg = arsIMG.izbrisiKratkeDeloveBinarneSlike(bwImg,40);
%             BW3 = bwmorph(bb,'skel',inf)        ;
            komponenteKonture  = bwconncomp(bwImg,4)  ;
            brKontura = komponenteKonture.NumObjects;
            imgSize   = komponenteKonture.ImageSize;
            maxPix=[0 0 ];
            minPix = maxPix;
            for iColumn = 1:1:imgSize(2)
                kol = bwImg(:,iColumn);
                  kol = find(kol==1);
                  if ~isempty(kol)% && (numel(kol)>1)
                    maxPix = [maxPix; [max(kol) iColumn]];
                    minPix = [minPix; [min(kol) iColumn]];
                  end
            end
            %
            if debugMode
                xMin = minPix(:,1);
                yMin = minPix(:,2);
                xMax = maxPix(:,1);
                yMax = maxPix(:,2);
                figure; imshow(bwImg);
                line(yMin,xMin);
                figure; imshow(bwImg);
                line(yMax,xMax);
            end
            presek = intersect([maxPix(:,2) maxPix(:,1)],[minPix(:,2) minPix(:,1)],'rows');
            if debugMode
                figure;
                xPresek = presek(:,1); yPresek = presek(:,2);
                imshow(bwImg); line(xPresek, yPresek);
            end
            %brisanje duplih elemenata
            pom1 = [ 0 0];
            for i=1:numel(presek(:,1))
                rezMaxPix = find(maxPix(:,2)==presek(i,1));
                maxPix(numel(maxPix(:,1))+rezMaxPix) = 0;
                rezMinPix = find(minPix(:,2)==presek(i,1));
                minPix(numel(minPix(:,1))+rezMinPix) = 0;
            end
            rezMaxPix = [ 0 0];
            for i=1:numel(maxPix(:,1))
                if maxPix(i,2)~=0
                    rezMaxPix = [rezMaxPix;maxPix(i,:)];
                end
            end
            rezMinPix = [ 0 0];
            for i=1:numel(minPix(:,1))
                if minPix(i,2)~=0
                    rezMinPix = [rezMinPix;minPix(i,:)];
                end
            end
            maxPix = rezMaxPix(2:end,:);
            minPix = rezMinPix(2:end,:);
            if debugMode
                rez    = maxPix;
                x = rez(:,1); y = rez(:,2);
                figure;imshow(bwImg);   line(y,x);
                rez    = minPix;
                x = rez(:,1); y = rez(:,2);
                figure;imshow(bwImg);   line(y,x);
            end
            spoljasnjeKonture=0;
            imgSpoljasnjeKonture = zeros(imgSize);
            unutrasnjeKonture=0;
            imgUnutrasnjeKonture = zeros(imgSize);
            for i=1:brKontura
                fff= arsIMG.daLiKonturaPripadaUnutrasnjojIliSpoljasnjoj(bwImg,komponenteKonture.PixelIdxList{i},maxPix,minPix);
                if fff==1
                    spoljasnjeKonture=[spoljasnjeKonture; i];
                    imgSpoljasnjeKonture(komponenteKonture.PixelIdxList{i})=1;
                elseif fff==2
                    unutrasnjeKonture=[unutrasnjeKonture; i];
                    imgUnutrasnjeKonture(komponenteKonture.PixelIdxList{i})=1;
                end
            end
%             figure; imshow(imgUnutrasnjeKonture);
%             figure; imshow(imgSpoljasnjeKonture);
            koordinateSpoljasnje  = arsIMG.pxelIDS2Koordinate(imgUnutrasnjeKonture,find(imgUnutrasnjeKonture>0));
            koordinateUnutrasnje = arsIMG.pxelIDS2Koordinate(imgSpoljasnjeKonture,find(imgSpoljasnjeKonture>0));
            rez.spoljasnjaKontura    = koordinateSpoljasnje;
            rez.imgUnutrasnjeKonture = imgUnutrasnjeKonture;
            rez.unutrasnjaKontura    = koordinateUnutrasnje;
            rez.imgSpoljasnjeKonture = imgSpoljasnjeKonture;
        end
        %%
        function [rez] = pxelIDS2Koordinate(img,PixelsID)
            [imgSizeX imgSizeY] = size(img);
            rezKontura1         = [0 0];
            rez(:,1)            = mod(PixelsID,imgSizeX);
            rez(rez==0)         = imgSizeX;
            rez(:,2)            = ceil(PixelsID/imgSizeX);
            rez(:,3)            = 0;
        end
        %%
        function [rez] = pixelsKoordinate2IDS(img,PixelsKoordinate)
            imgSize= size(img);
            rezKontura1=1;
            for i = 1:numel(PixelsKoordinate(:,1))
                rezKontura1=[rezKontura1; PixelsKoordinate(i,1)+(PixelsKoordinate(i,2)-1)*imgSize(1)];
            end
            rez = rezKontura1(2:end);
        end
        %%
        function [rez] = daLiKonturaPripadaUnutrasnjojIliSpoljasnjoj(img, kontura, kontura1, kontura2)
            imgSize= size(img);
            rezKontura1 = arsIMG.pixelsKoordinate2IDS(img,kontura1);
            rezKontura2 = arsIMG.pixelsKoordinate2IDS(img,kontura2);
            unionKonturaRezKontura1 = unique([kontura;rezKontura1]);
            unionKonturaRezKontura2 = unique([kontura;rezKontura2]);
%             figure;
%             img = zeros(imgSize);
%             img(kontura) = 1;
%             imshow(img);
%             x = kontura1(:,2); y = kontura1(:,1); line(x,y);
%             x = kontura2(:,2); y = kontura2(:,1); line(x,y);
            if (numel(unionKonturaRezKontura1)-numel(rezKontura1)==0 &&...
                numel(unionKonturaRezKontura2)-numel(rezKontura2)==0)
                rez = 0;
            elseif(numel(unionKonturaRezKontura1)-numel(rezKontura1)>...
                    numel(unionKonturaRezKontura2)-numel(rezKontura2))
%                 rez = 'pripada konturi 1-spoljasnjem zidu';
%                 title(rez);
                rez = 1;
            else
%                 rez = 'pripada unutrasnjem zidu';
%                 title(rez);
                rez = 2;
            end
        end
        %%
        function [rez]= iscrtajKonturuUDekartovomKooSistemu(imgPath, uglovi,precnici)
            %interpoliraj konturu(polarni kooordinatni sistem)
%             precnikUgao = obj.BSpline.interpolateFromAtoB(obj.BSpline, 0.01,0.01, 0.99);
            if isstr(imgPath)
                imgData = imread(imgPath);
            else
                imgData = imgPath;
            end            
            imgData= imgData(:,:,1);            
            imgSize = size(imgData);
            precnikUgao = [ uglovi precnici ];
            precnik = precnikUgao(2:end-1,2)       ;%izdvoj koordinate precnika
            radijus = imgSize(1)/2-precnik         ;%koriguj precnik, u matlabu Height=0 je na vrhu slike
            ugao    = precnikUgao(2:end-1,1)       ;%izdvoj ugao
            ugao    = (ugao/180) * pi              ;%konvertuj u radijane
            cosUgao = cos(ugao)                    ;%izracunaj kosinus uglova
            sinUgao = sin(ugao)                    ;%izracunaj sinus uglova
            %izracunaj X-Y koordinate i zatvori konturu
            x       = cosUgao .* radijus           ; 
            x       = [ x ; x(1)]                  ;
            y       = sinUgao .* radijus           ;
            y       = [ y ; y(1)]                  ;
            %prikazivanje rezultate
%             s       = imread('img40frame.jpg')     ;
            s       = imgSize
            x       =  x + s(2)/2                  ;%centriraj na sliku
            y       = -y + s(1)/2                  ;
            %imshow('img40frame.jpg'); hold on;
            f = figure;
            imshow(imgData);
            hold on; line(x,y,'LineWidth',1); axis equal; z=x; z(z~=0)=0;
%             BSpline = arsBSplineOpenNonuniform;
%             BSpline = BSpline.set(BSpline, [x,y,z], 3);
%             BSpline.test(BSpline);
            fit_ellipse(x,y,f);
            rez = [x,y,z];
        end  
        %% gausovo smutovanje slike
            %INPUTS
                %img   - slika koja se smutuje
                %hsize - [x y] velicina filtera
                %sigma - standardna devijacija
        function [rez] = gaussianSmoothing(img, hsize, sigma)
            filter = fspecial('gaussian',hsize,sigma);
            rez    = imfilter(img,filter,'same')     ;
        end        
        %% Lokalni Binarni Patern
        % ref http://www.cse.oulu.fi/CMV
        function [rez] = LocalBinaryPattern(img, R, N)
%          imgPolar  = arsIMG.dekartToPolar(img);
           imgPolar  = img             ;
           %jedinicni krug
           ugao      = 0:pi/N:pi-0.1   ;
           xKruga    = cos(ugao)*R     ;
           yKruga    = sin(ugao)*R     ;
           %dimenzija polarne slike
           dimImg    = size(imgPolar)  ;
           LBP_image = zeros(dimImg)   ;            
            for i = R+1:dimImg(1)-R
                for j= R+1:dimImg(2)-R
                pozicijaTrenutnogCentra= [i,j];
                    for k = 1:numel(xKruga)
                        LBP_razlika(k) = arsIMG.getPixelBipolarnaInterpolacija(...
                                         imgPolar,pozicijaTrenutnogCentra + [xKruga(k) yKruga(k)])...
                                         -imgPolar(i,j);
                    end%                                      
                LBP_razlika(LBP_razlika>0) = 1                                    ;
                LBP_razlika(LBP_razlika<0) = 0                                    ;
                LBP_image(i,j)             = AngioIvusMath.bin2Double(LBP_razlika);
                end
            end 
          rez = LBP_image;
        end
        %%
        function [dominantniUgao]= getDominantniGradijentUgao(img)
            [xGrad, yGrad]    = arsIMG.izracunajGradient(img)              ;
            ugao              = atan(yGrad ./ xGrad)                       ;
%             ugao = rad2deg(ugao); za preveru konvertovano u stepene
            ugao(isnan(ugao)) = 666                                        ;
            %
            jedinstveniUglovi = ugao                                       ;
            intervali         = [-pi/2:0.2:pi/2]                           ;
            for i=1:numel(intervali)
                jedinstveniUglovi(abs(jedinstveniUglovi-intervali(i))<0.1)=intervali(i) ;
                pom   = jedinstveniUglovi==intervali(i)                                 ;                
                br(i) = sum(sum(pom))                                                   ;
            end
            x = max(br);
            %nadji maksimalan interval
            maxS = 0;
            for i=2:numel(br)-1
                s = sum(br(i-1:i+1));
                if s>maxS
                    maxS = s;
                    indexMax= i;
                end
            end
%             dominantniUgao = intervali(find(br==x));
            brUglova = br(indexMax-1:indexMax+1)       ;
            ugovi    = intervali(indexMax-1:indexMax+1);
            dominantniUgao = sum( brUglova .* ugovi) / sum(brUglova);
            dominantniUgao = rad2deg(mean(dominantniUgao));
        end
        %%
        function pripremiSlikeZaIVUS(osnovaSlike)
            for i=1:500
                img1Path = [osnovaSlike '_'  num2str(i,'%.4i') '.jpg']    ;
                img1     = imread(img1Path)                               ;
                img2Path = [osnovaSlike '_'  num2str(i+1,'%.4i') '.jpg']  ;
                img2     = imread(img2Path)                               ;
                img      = img1+img2                                      ;
                %
                imwrite(img, ['img_'  num2str(i,'%.4i')  '.jpg'] );
            end
        end
        %% konvertovanje slike u polarni koo.sistem
        function [imgPolar]= dekartToPolar(img,O)
            if isstr(img)
                img = imread(img);
            end
            img     = mat2gray(img);
            imgSize = size(img)    ;
            if(numel(imgSize)==3) %ima trecu komponentu
                img = ( img(:,:,1) + img(:,:,2) + img(:,:,3) )/3; 
            end
            %img         = img(30:500,80:550);
            if nargin ==1                
                O           = round(size(img) ./ 2)  ; 
                totalRadius = round(min(size(img))/2);
            else 
                O           = round(O);
                totalRadius = min(O(1),imgSize(1)-O(1));
                totalRadius = min(totalRadius,min(O(2),imgSize(2)-O(2)));
            end                        
            imgPolar    = zeros(totalRadius,360) ;
            
            for radius = 1:totalRadius-1
                for ugao = 1:360
                    ugaoRadiani = (ugao-1)*pi/180;
                        %width
                        x = O(1)-round(radius*sin(ugaoRadiani));
                        %height
                        y = O(2)+round(radius*cos(ugaoRadiani));
                        newPixel = img(x,y);
                    imgPolar(totalRadius-radius,ugao) = newPixel;                    
                end
            end%for
            %imshow(imgPolar);
        end%function set newPoints
        %% konvertuje sliku iz polarnih u dekartov koo. sistem
        function [rez] = polar2dekart(img)
            imgSize = size(img);            
            rez     = zeros(2*imgSize(1)+1,2*imgSize(1)+1);
            for i = 1:imgSize(1)%radius
                for j = 1:imgSize(2)%theta
                    %[X,Y] = pol2cart(THETA,RHO)
            %arsIMG.polarToDekartKoordinate(koordinate, centarSlike)
            %koordinate [ radius, ugao]
                    xy = arsIMG.polarToDekartKoordinate([imgSize(1)-i, j], [imgSize(1), imgSize(1)]);
                    x = round(xy(1)); y = round(xy(2));
                    if(x==0)
                        x=1;
                    end
                    if(y==0)
                        y=1;
                    end
                    try
                    rez(x,y) = img(i,j);
                    catch
                        [i j] 
                        [ x y]
                    end
                end
            end
            imshow(rez);
        end
        %% uzima piksele iz slike pri cemu polozaj piksela moze biti realan broj
        %interpolacija piksela se vrsi preko bipolarne interpolacije
        %ref:http://cboard.cprogramming.com/game-programming/87890-bilinear-interpolation.html
        %INPUTS:
            %img - slika sa koje se uzima interpolirani pixel
            %pointSelektor - koordinate piksela koji se trazi
        %OUTPUTS:
            %rezPiksel
        function [rezPiksel] = getPixelBipolarnaInterpolacija(img,pointSelektor)
            imgSize = size(img);
            if (min(pointSelektor)<1 || max(pointSelektor)>imgSize(1) || max(pointSelektor)>imgSize(1))
                rezPiksel=0;
                return;
            end
            xUpper = ceil(pointSelektor(1)) ;
            xLower = floor(pointSelektor(1));
            yUpper = ceil(pointSelektor(2)) ;
            yLower = floor(pointSelektor(2));
            % O(x=0,y=0)img
            %   (xLower,yLower)A........B(xUpper,yLower)   ----->Xosa(br.kolone  
            %                   :      :                   |
            %                   :   O  :                   |
            %   (xLower,yUpper)C........D(xUpper,yUpper)   V Yosa(br.reda
            A = img(xLower,yLower);
            B = img(xUpper,yLower);
            C = img(xLower,yUpper);
            D = img(xUpper,yUpper);
            p = (pointSelektor(1)-xLower);
            q = (pointSelektor(2)-yLower);
            rezPiksel = (1-p)*(1-q)*A... 
                         +(1-p)*q*C...
                         +p*(1-q)*B...
                         +p*q*D;  
        end        
        %% %%%%% The functions used in the main.m file %%%%%%%
        % Function "d2gauss.m":
        % This function returns a 2D Gaussian filter with size n1*n2; theta is 
        % the angle that the filter rotated counter clockwise; and sigma1 and sigma2
        % are the standard deviation of the gaussian functions.
        function h = d2gauss(n1,std1,n2,std2,theta)
        r=[cos(theta) -sin(theta);
           sin(theta)  cos(theta)];
        for i = 1 : n2 
            for j = 1 : n1
                u = r * [j-(n1+1)/2 i-(n2+1)/2]';
                h(i,j) = arsIMG.gauss(u(1),std1)*arsIMG.gauss(u(2),std2);
            end
        end
        h = h / sqrt(sum(sum(h.*h)));
        end
        %% Function "gauss.m":
        function y = gauss(x,std)
        y = exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));
        end
        %% smoothuj sliku na pravcu
        function [rez]=smoothDeoSlikuPoPravcu(img, pravac, tackaSmoothovanja, dimenzijeSmoothRegiona)
            %jedinicni vektori pravca selekcije
            pravacRadiani   = deg2rad(pravac)                              ;
            matricaRotacije = [cos(pravacRadiani),-sin(pravacRadiani);...
                               sin(pravacRadiani),cos(pravacRadiani)];
            x = [matricaRotacije*[1;0]]'                                   ;
            y = [matricaRotacije*[0;1]]'                                   ;
%             %prikazi sliku i region koji je selektovan
%             imshow(img);
%             A = tackaSmoothovanja - x*dimenzijeSmoothRegiona(1) + y*dimenzijeSmoothRegiona(2);
%             B = tackaSmoothovanja + x*dimenzijeSmoothRegiona(1) + y*dimenzijeSmoothRegiona(2);
%             C = tackaSmoothovanja - x*dimenzijeSmoothRegiona(1) - y*dimenzijeSmoothRegiona(2);
%             D = tackaSmoothovanja + x*dimenzijeSmoothRegiona(1) - y*dimenzijeSmoothRegiona(2);
%             plotLine([A,0],[B,0]);
%             plotLine([A,0],[C,0]);
%             plotLine([C,0],[D,0]);
%             plotLine([D,0],[B,0]);  
            xDim    = 2*dimenzijeSmoothRegiona(1)*cos(pravacRadiani) + 2*dimenzijeSmoothRegiona(2)*sin(pi/2-pravacRadiani)      ;
            xDim    = round(xDim)                                                                                               ;
            yDim    = 2*dimenzijeSmoothRegiona(1)*sin(pravacRadiani) + 2*dimenzijeSmoothRegiona(2)*cos(pi/2-pravacRadiani)      ;
            yDim    = round(yDim)                                                                                               ;
            matricaRotiraniGauss2D = arsIMG.d2gauss(xDim,dimenzijeSmoothRegiona(2),yDim,dimenzijeSmoothRegiona(1),pravacRadiani);            
        end
        %% Racuna gradijent dela slike selektujuci piksele po zadatom
        %  pracvcu. Funkcija se nadovezuje na arsIMG.getDominantniGradijentUgao
            %INPUTS
                %img                    - slika
                %pravac                 - ugao(degree)za koji je region selekcije rotiran od x ka y   Y|-->x,  
                %tackaSmoothovanja      - [x-br reda V,y-br kolone-->] centar pravougonika smoothovanja
                %dimenzijeSmoothRegiona - [xWidth-po kolonama,yHeight-po redovima] dimenzije pravougaonika smoothovanja
        %          X1
        %    y|   /^.
        %    /|  /   )pravac
        %   / |\/___/_
        %  / A|/\      gradient =  A - B
        %  \ c!--\---------------------------------->x0
        %   \/ B /    
        %   /\  /    c(x,y) = tackaSmoothovanja
        %  /  \/
        function [rez]=izracunajGradientDelaSlikePoPravcu(img, pravac, tackaSmoothovanja, dimenzijeSmoothRegiona)
            %jedinicni vektori pravca selekcije
            pravacRadiani   = deg2rad(pravac)                              ;
            matricaRotacije = [cos(pravacRadiani),-sin(pravacRadiani);...
                               sin(pravacRadiani),cos(pravacRadiani)]      ;
            x               = [matricaRotacije*[1;0]]'                     ;
            x               = [-x(2), x(1)]                                ;
            y               = [matricaRotacije*[0;1]]'                     ;
            y               = [-y(2), y(1)]                                ;
            %selektuj tacke i smoothuj
            smoothRegion    = zeros([dimenzijeSmoothRegiona(2),dimenzijeSmoothRegiona(1)]*2+[1,1])        ;
            for i=-dimenzijeSmoothRegiona(1):dimenzijeSmoothRegiona(1)
                for j=-dimenzijeSmoothRegiona(2):dimenzijeSmoothRegiona(2)
                    iX                =  i*x                               ;
                    jY                = -j*y                               ;
                    koordinatePiksela = iX+jY+tackaSmoothovanja            ;
                    piksel = arsIMG.getPixelBipolarnaInterpolacija(img, koordinatePiksela);
                    smoothRegion(dimenzijeSmoothRegiona(2)+j+1,dimenzijeSmoothRegiona(1)+i+1) = piksel;
                end
            end
            %prikazi sliku i region koji je selektovan
            f = figure;
            imshow(img);
            x=[-x(2), x(1)];y=[-y(2), y(1)];
            A = tackaSmoothovanja + x*dimenzijeSmoothRegiona(2) + y*dimenzijeSmoothRegiona(1); 
            A = fliplr(A);
            B = tackaSmoothovanja - x*dimenzijeSmoothRegiona(2) + y*dimenzijeSmoothRegiona(1);
            B = fliplr(B);
            C = tackaSmoothovanja + x*dimenzijeSmoothRegiona(2) - y*dimenzijeSmoothRegiona(1);
            C = fliplr(C);
            D = tackaSmoothovanja - x*dimenzijeSmoothRegiona(2) - y*dimenzijeSmoothRegiona(1);
            D = fliplr(D);
            plotLine([A,0],[B,0]);
            plotLine([A,0],[C,0]);
            plotLine([C,0],[D,0]);
            plotLine([D,0],[B,0]);  
            f1 = figure;
            pom = transp(smoothRegion);pom = pom(end:-1:1,:);
            imshow(pom);
            close(f); close(f1);
            %nadji gradijente 
            razlika = 0;%;smoothRegion(1,:);
            for i = 1: dimenzijeSmoothRegiona(1)-1
                if pravac>0
                    pom     = smoothRegion(:,dimenzijeSmoothRegiona(1)+i)-smoothRegion(:,dimenzijeSmoothRegiona(1)-i);
                else
                    pom     = smoothRegion(:,dimenzijeSmoothRegiona(1)-i)-smoothRegion(:,dimenzijeSmoothRegiona(1)+i);
                end
                pom     = sum(sum(pom));
                razlika = razlika+pom;
            end
            % totalni gradijent, oduzima sumu jednog bloka od drugoga
%             a=smoothRegion(1:dimenzijeSmoothRegiona(2),1:dimenzijeSmoothRegiona(1));
%             b=smoothRegion(dimenzijeSmoothRegiona(2)+1:2*dimenzijeSmoothRegiona(2),dimenzijeSmoothRegiona(1)+1:2*dimenzijeSmoothRegiona(1));
%             gradient = sum(sum(a-b));
            rez = razlika;
        end
        %% konvertovanje slike u polarni koo.sistem
        % INPUTS:
        %   img : slika, moze putanja moze i data
        %   CentarSlike : centarSlike-koo pocetak za polarnu slikuu
        %   Radijus : radijus slike
        %   rezolucijaPiksela : normalno je 1, moze i <1 i >1
        function [imgPolar]= dekartToPolarBipolarnaInterpolacija(img,CentarSlike, Radius, rezolucijaPiksela,rezolucijaUglova)
            if isstr(img)
                img = mat2gray(imread(img));
            end
            img     	   = mat2gray(img)                                   ;
            imgSize        = size(img)                                       ;
            if(numel(imgSize)==3) %ima trecu komponentu
                img = ( img(:,:,1) + img(:,:,2) + img(:,:,3) )/3; 
                imgSize        = size(img);
            end
            %
            if nargin==1 %ako nisu prosledjeni svi parametri radi default
                CentarSlike       = round(imgSize/2)     ;
                Radius            = min(CentarSlike)     ;
                rezolucijaPiksela = 1                    ;
                rezolucijaUglova  = 1                    ;
                CentarSlike       = imgSize/2            ;
            end            
            %
            radiusSelector = [1:rezolucijaPiksela:Radius]                    ; 
            ugaoSelector   = [1:rezolucijaUglova:360]                        ;
            imgPolar       = zeros(numel(radiusSelector),numel(ugaoSelector));
            O              = CentarSlike                                     ;     
            %brojaci za kolone/redove za slike
            iRadius =0;
            iUgao   =0;
            for radius = 1:numel(radiusSelector)-1
                iRadius  = iRadius+1                                      ;
                iUgao    = 0                                              ;
                for ugao = 1:rezolucijaUglova:360-1
                 iUgao   = iUgao+1                                        ;
                 %bilinearna interpolacija, funkcija arsIMG.getPixelBipolarnaInterpolacij
                 %http://cboard.cprogramming.com/game-programming/87890-bilinear-interpolation.html
                 ugaoRadiani   = deg2rad(ugao);
                 pointSelektor = [ radius*cos(ugaoRadiani) ...%x koordinate
                                   radius*sin(ugaoRadiani)];  %y
                 pointSelektor = [-pointSelektor(2)  pointSelektor(1)] + O;
                 imgPolar(iRadius,iUgao) = arsIMG.getPixelBipolarnaInterpolacija(img, pointSelektor);                                          
                end
            end%for
            imgPolar = imgPolar(end:-1:1,:);
            imshow(imgPolar);
        end%function set newPoints
        %%
        %prebacuje koordinate iz dekartovog koo sistema u polarni
        %za centar se uzima centar slike
        %koordinate=[VERTIAKLNA, HORIZONTALNA]
        function [rez] = dekartToPolarKoordinate(koordinate, centarSlike)
            if nargin==1
                centarSlike = [0 0 0];
            end
            brKoordinata = size(koordinate(1,:,:));
            if brKoordinata(2)== 2
                rez = [0 0];
                for i = 1:numel(koordinate(:,1))
                    koordinate(i,:) = [centarSlike(1)-koordinate(i,1),centarSlike(2)-koordinate(i,2)] ;
                    [ugao, radius]  = cart2pol(koordinate(i,1),koordinate(i,2));  
                    ugao            = rad2deg(ugao)              ;
                    rez             = [ rez; radius, ugao]       ; % proveri nije tacno
                end
                rez = rez(2:end,:);
            else
                rez = [0 0 0];
                for i = 1:numel(koordinate(:,1,:))
                    koordinate(i,:,:) = [centarSlike(1)-koordinate(i,1,:),centarSlike(2)-koordinate(i,2,:), 0] ;
                    [ugao, radius]  = cart2pol(koordinate(i,1,:),koordinate(i,2,:));  
                    ugao            = rad2deg(ugao)              ;
                    rez             = [ rez; radius, ugao,koordinate(i,3,:)]    ; % proveri nije tacno
                end
                rez = rez(2:end,:,:);
            end
            id = find(rez(:,2)<0);
            rez(id,2)=360-abs(rez(id,2));
        end
        %% polarne u dekartove
        %INPUTS
            % koordinate [ radius, ugao]
            % centarSlike [ row, column]
        %OUTPUTS
            % rez[ row, column ]
        function [ rez] = polarToDekartKoordinate(koordinate, centarSlike) 
            if numel(koordinate(1,:,:))==2
            rez=[0 0 ];
            for i =1:numel(koordinate(:,1,:))
                red       =  centarSlike(:,1) - koordinate(i,1)*sin(deg2rad(koordinate(i,2)));
                kolona    =  centarSlike(:,2) + koordinate(i,1)*cos(deg2rad(koordinate(i,2)));
                rez = [rez; red , kolona]; %red=V, kolona=H
            end
            rez = rez(2:end,:);
            else
                rez=[0 0 0];
            for i =1:numel(koordinate(:,1,:))
                red     =  centarSlike(:,1) - koordinate(i,1,:)*sin(deg2rad(koordinate(i,2,:)));
                kolona  =  centarSlike(:,2) + koordinate(i,1,:)*cos(deg2rad(koordinate(i,2,:)));
                rez     =  [rez; red ,kolona, 0];
            end
            rez = rez(2:end,:,:);
            end
        end
        %% racunanje gradienta slike
        function [xGrad, yGrad] = izracunajGradient(img,dt)
            imgSize = size(img)     ;
            %inicijalizacija matrica gde ce se cuvati slika
            xGrad   = zeros(imgSize);
            yGrad   = zeros(imgSize);
            if nargin == 1
                dt=1;%korak
            end
            %
            [xGrad,yGrad] = imgradientxy(img);
% % % % %             for row = (dt+1):imgSize(1)-dt
% % % % %                 for column = (dt+1):imgSize(2)-dt
% % % % % %                      img(row,column)
% % % % % %                      y1 = img(row-dt,column   )
% % % % % %                      y2 = img(row+dt,column   )
% % % % % %                      y2 - y1
% % % % %                     yGrad(row,column) = (img(row+dt,column   ) - img(row-dt,column   ))/2.;
% % % % %                     xGrad(row,column) = (img(row   ,column+dt) - img(row   ,column-dt))/2.;
% % % % % %                     x2 = img(row   ,column+dt)
% % % % % %                     x1 = img(row   ,column-dt)
% % % % % %                     x2 - x1
% % % % % %                 end
% % % % % %             end
        end
        %% racunanje gradienta slike
        function [xGrad, yGrad] = izracunajGradientOCT(img)
            imgSize = size(img)     ;
            %inicijalizacija matrica gde ce se cuvati slika
            xGrad   = zeros(imgSize);
            yGrad   = zeros(imgSize);
            dt      = 4             ;%korak
            %
            for row = (dt+1):imgSize(1)-dt
                for column = (dt+1):imgSize(2)-dt
                    xGrad(row,column) = img(row-dt,column   ) - img(row+dt,column   );
                    if(xGrad(row,column) < 0.1)
                        xGrad(row,column) = 0;
                    end
                    yGrad(row,column) = sum(img(row   ,column:column-dt)) - sum(img(row   ,column:column+dt));
                    if(yGrad(row,column) < 0.1)
                        yGrad(row,column) = 0;
                    end
                end
            end
            %prikazi rezultate \za debugg
            %f   = figure; ttl = ['original ']; subplot(2,3,2), subimage(img), title(ttl);ttl = ['xGrad ']; subplot(2,3,4), subimage(xGrad), title(ttl); ttl = [ 'yGrad ' ]; subplot(2,3,6), subimage(yGrad), title(ttl);            
            %close(f);
        end
%         %% racunanje gradienta slike
%         function [xGrad, yGrad] = izracunajGradientIVUS(img)
%             imgSize = size(img)     ;
%             %inicijalizacija matrica gde ce se cuvati slika
%             xGrad   = zeros(imgSize);
%             yGrad   = zeros(imgSize);
%             dKolona  = 2             ;%korak
%             dRow   = 6;
%             %
%             for row = (dRow+1):imgSize(1)-dRow
%                 for column = (dKolona+1):imgSize(2)-dKolona
%                     xGrad(row,column) = img(row-dRow,column   ) - img(row+dRow,column   );
%                     if(xGrad(row,column) < 0.1)
%                         xGrad(row,column) = 0;
%                     end
%                     yGrad(row,column) = mean(mean(img(row-dRow:row,  column-dKolona:column) -...
%                                                  img(row:row+dRow,  column:column+dKolona)));
%                     if(yGrad(row,column) < 0.1)
%                         yGrad(row,column) = 0;
%                     end
%                 end
%             end
%             %obrisi sum
%             imshow(yGrad);
%             vGrad = arsIMG.removeNoiseIVUS(xGrad);
%             imshow(vGrad);
%             %prikazi rezultate \za debugg
%             %f   = figure; ttl = ['original ']; subplot(2,3,2), subimage(img), title(ttl);ttl = ['xGrad ']; subplot(2,3,4), subimage(xGrad), title(ttl); ttl = [ 'yGrad ' ]; subplot(2,3,6), subimage(yGrad), title(ttl);            
%             %close(f);
%         end
%         %%
%         function [rez] = izracunajGradientIVUS2(img)
%             %idemo redom i racunamo koliko je slicno lumenu
%         end
        %%%%%%%%%%
        function plotIntenzitete(img)
            dim = size(img);
            visina =dim(1); duzina = dim(2);
            for i =1:visina
                x = [1:duzina];
                y = img(i,:);
                z = x; z(z>0)=i;
                hold on;
                line(x,y,z);
            end
        end
        %% racunanje gradienta slike, trebalo bi da radi brze
        function [xGrad, yGrad] = izracunajGradient2(img)
            imgSize = size(img)     ;
            %inicijalizacija matrica gde ce se cuvati slika
            dt      = 2             ;%korak
            xGrad   = zeros(imgSize);%gradijenti po horizontali
            xGrad(:,1+dt:end-dt) = img(:,1:end-2*dt) - img(:,1+2*dt:end);
            xGrad   = abs(xGrad);
            %
            yGrad   = img'; 
            yGrad(:,1+dt:end-dt) = yGrad(:,1:end-2*dt) - yGrad(:,1+2*dt:end);
            yGrad   = yGrad';
            yGrad   = abs(yGrad);
            %prikazi rezultate \za debugg
            %f   = figure; ttl = ['original ']; subplot(2,3,2), subimage(img), title(ttl);ttl = ['xGrad ']; subplot(2,3,4), subimage(xGrad), title(ttl); ttl = [ 'yGrad ' ]; subplot(2,3,6), subimage(yGrad), title(ttl);            
            %close(f);
        end
        %% racunanje gradienta slike, pod uglom Fi-za IVUS
        %INPUTS
            %img-slika
            %Fi -ugao
            %dim-dimenzije filtera
        function [xGrad, yGrad] = izracunajGradientPodUglomFi(img,fi,dim)
            debugMode = 0             ;
            imgSize   = size(img)     ;
            %   I----------I  ^y=row
            %   I 1 1 1 1 1I  |
            %   I----------I--|----->x=column
            %   I-1-1-1-1-1I
            %   I----------I
            dRow      = dim(1)        ;%korak
            dCol      = dim(2);       
            %%%%%%%%%%%%%%%%%%
            filtrY =  ones(dim);
            filtrX =  ones(dim);
            filtrY(round(dRow/2)+1:end,:) = -1;
            filtrX(:,round(dCol/2)+1:end) = -1;
            %zarotiraj filter za ugao
            %ref:http://www.mathworks.com/help/toolbox/images/ref/imrotate.html
            filtrY = imrotate(filtrY,fi,'bicubic');
            filtrX = imrotate(filtrX,fi,'bicubic');
            yGrad  = imfilter(img, filtrY); 
            xGrad  = imfilter(img, filtrX);
            if debugMode
                imshow(yGrad); hold on; imshow(filtrY); title('gradient Y');
                figure;
                imshow(xGrad); hold on; imshow(filtrX); title('gradient X');
                close all;
            end
        end
        %% odklanjanje suma, tj. pixela koji su usamljeni
        function [rez] = removeNoise(img)           
            %ako je avg suma oko izabranog pixela manja od njega onda je taj pixel crn 
            dt = 4; % blok je 5x5 piksela
            imgSize = size(img);
            rez = img;
            for row = (dt+1):imgSize(1)-dt
                for column = (dt+1):imgSize(2)-dt
                    sampleBlock = img( (row-dt):(row+dt),(column-dt):(column+dt));
                    %imshow(sampleBlock);
                    avgSum = mean(mean(sampleBlock));
                    %if  img(row,row) < avgSum
                    if  avgSum<0.07 %bilo 0.1
                        rez(row,column) = 0;
                    else
                        %imshow(sampleBlock);
                        
                    end%if                    
                end
            end           
            %imshow(rez); 
            %figure; imshow(img)
        end
        %% odklanjanje suma, tj. pixela koji su usamljeni
        function [rez] = removeNoiseOCT(img)           
            %ako je avg suma oko izabranog pixela manja od njega onda je taj pixel crn 
            dt        = 4            ; % blok je 5x5 piksela
            imgSize   = size(img)    ;
            rez       = img          ;
            kandidati = find(img>0.1);
            obrisani  = 1            ;
            for i = 1:numel(kandidati)                
                row    = mod(kandidati(i),imgSize(1))   ;
                column = (kandidati(i)-row)/imgSize(1)+1; 
                if(column>dt && column<(imgSize(2)-dt)) && (row>dt && row<(imgSize(1)-dt))
                    sampleBlock = img( (row-dt):(row+dt),(column-dt):(column+dt));
                    avgSum      = mean(mean(sampleBlock))                        ;
                    if  avgSum<0.07 %bilo 0.1
                        rez(row,column) = 0;                        
                    end%if    
                end                
            end               
            rez(rez<0.1)=0;
        end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [rez] = removeNoiseIVUS(img, min, drow,dcolumn)   
            if nargin==2
                drow    = 5;
                dcolumn = 5;
            end
            fltr = ones(drow,dcolumn);
            fltr = fltr / numel(fltr);
            rez  = imfilter(img, fltr);
%             imshow(rez);
            rez(rez<min)=0;
%             imshow(rez);
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function napraviKlipZaArThreat(aviPath, pozadina)
            aviInfo              = mmreader(aviPath)     ;
            ukupnoFrejmova       = aviInfo.NumberOfFrames;
            moviename = 'a11a.avi';
            aviobj = avifile(moviename,'compression','None','fps','2');
            pozadina = imread(pozadina);
            for brFrejma =1:ukupnoFrejmova
                img     = read(aviInfo,brFrejma);
                F       = pozadina              ;
                F(1:696,1:900,:)=img;
                aviobj  = addframe(aviobj,F)    ;                
            end
             aviobj = close(aviobj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% vrsi se segmentacija kroz frameove
        function [rez, rez2] = getFramesFromAviOCT(aviPath)
            %selektovanj BSplajna iz 2 angio
%             biPlane = arsBiPlaneAngio;
%             biPlane = biPlane.test(biPlane);
            biPlane = arsCubicBSplineOpenNonuniform         ;
            points  = [ 0 0 0; 400 400 400; 100 200 1500;  ];
            biPlane = biPlane.set(biPlane, points,3)        ;
            rez     = [ 0 0 0]; rez2 =rez;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            brPocetnogFrejma     = 1                     ;
            brZadnjegFrejma      = 260                   ;
            aviInfo              = mmreader(aviPath)     ;
            ukupnoFrejmova       = aviInfo.NumberOfFrames;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            snake = arsOpenSnake;
            snake = snake.setSettingsDefault(snake);
            %brPocetnogFrejma
%             snake = snake.setImgDataZaOCT(snake,read(aviInfo,brPocetnogFrejma));
%             snake = snake.setContourPointsOCT(snake)  ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %pravljenje avi fajla
               aviobj = avifile('simone_2.avi','compression','None','fps','1');
               aviobj.KeyFramePerSec = 1;
            %ucitavanje *.dat fajl 
            fid = fopen('exp1.txt', 'w');
            daLiTrebaPravitiFajl = 1;
            korak = 3;
            %upisivanje konture u fajla za pak
            if daLiTrebaPravitiFajl
                fprintf(fid, 'Number of branches\n');
                fprintf(fid, '2\n');
                fprintf(fid, 'Branches\n');
                fprintf(fid, 'ID    ParentID\n');
                fprintf(fid, '2     0\n');
                fprintf(fid, '1     2\n');
                fprintf(fid, 'Number of slices:\n');
                fprintf(fid, '%d\n', round(brZadnjegFrejma/korak));   
            end
            
            brSlajsa = 0;
            for brFrejma=brPocetnogFrejma:korak:brZadnjegFrejma%ukupnoFrejmova
                brFrejma                
                fig   = figure;
                imshow(read(aviInfo,brFrejma));
                %postavi podatke u snake-u
                snake = snake.setImgDataZaOCT(snake,read(aviInfo,brFrejma));
                snake = snake.setContourPointsOCT(snake)  ;   
                %iscrtaj konturu
                konturaSlajsa = arsOCT.iscrtajKonturuUDekartovomKooSistemu(snake);
                konturaSlajsa(:,3,:) = brFrejma*10;
                rez = [ rez; konturaSlajsa];
                x = konturaSlajsa(:,1,:); y = konturaSlajsa(:,2,:); z = konturaSlajsa(:,3,:);
                x = [ x;x(1)]; y = [ y;y(2)]; z = [ z;z(3)];
                hold on; line(x,y,z)
                %%%%%%%%%%%%% upisivanje u fajl
                if daLiTrebaPravitiFajl%da li se generise dat fajl
                    fprintf(fid, 'Slice:\n');
                    fprintf(fid, '%d\n', brSlajsa);
                    brSlajsa = brSlajsa+1;
                    fprintf(fid, 'Branch count\n');
                        fprintf(fid, '%d\n', 2);
                    for brBranch =1:2                    
                        fprintf(fid, 'Branch \n');
                        fprintf(fid, '%d\n', brBranch);
                        fprintf(fid, 'Number of points:\n');
                        fprintf(fid, '%d\n', numel(konturaSlajsa(:,1,:)));
                        konturaNaSplajnu = biPlane.getTackePresekaAuTrihedronSplajnaAtT(biPlane,konturaSlajsa,brFrejma/brZadnjegFrejma);
                        rez2 = [ rez2; konturaNaSplajnu];
                        for brTackeKonture = 1:numel( konturaNaSplajnu(:,1,:))
                            tacka = konturaNaSplajnu(brTackeKonture,:,:);
                            fprintf(fid, '%.6f %.6f %.6f\n', tacka(1),tacka(2),tacka(3));
                        end                    
                    end 
                end             
                %upisi u avi fajl
                   F       = getframe(fig)         ;
                   F.cdata = F.cdata(1:760,1:845,:);                
                   aviobj  = addframe(aviobj,F)    ;
                   close(fig);
                %
            end%for
               aviobj = close(aviobj);
               fclose(fid);%dat fajl
        end%getFrameVromAvi
        %% napravi AVI od niza jpg slika-frameova iz dicoma
        function napraviAVIizNizaJPG(filename, moviename)
            %filename osnova naziva jpg slika, na nju se dodaju samo 12345
            %moviename naziv klipa
            aviobj = avifile(moviename,'compression','None','fps','1');
               aviobj.KeyFramePerSec = 1;
               for i = 1:500
                 nazivSlike    = num2str(i,'%04i')               ;
                 nazivSlike    = [filename '_' nazivSlike '.jpg'];
                 trenutnaSlika = imread(nazivSlike)              ;
                 aviobj        = addframe(aviobj,trenutnaSlika)  ;
               end
               aviobj = close(aviobj);
        end
        %%
        function napraviAVIizNizaSlika(NizSlika, moviename)
            outputVideo = VideoWriter(moviename);
            outputVideo.FrameRate = 3;
            open(outputVideo);
            [pom pom N] = size(NizSlika);
            for i = 1:N
               img = mat2gray(NizSlika(:,:,i));
               writeVideo(outputVideo,img)
            end
            close(outputVideo);
        end
        %% 
        function [rez] =  MeanSquaredDifference(originalBlok, kandidatBlok)
            rez =  originalBlok - kandidatBlok;
            rez = rez .* rez                         ;
            rez = sum(sum( rez))                     ;
            MN = size(originalBlok);
            if isnan(rez)%akoje NaN stavi da bude nula
               rez = 0;
            end
            rez = rez / (MN(1)*MN(2));
        end
        %% The Mean Absolute Difference (MAD)
         function [rez] =  MeanAbsoluteDifference(originalBlok, kandidatBlok)
            rez = originalBlok - kandidatBlok        ;
            rez = abs(rez)                           ;
            rez = sum(sum( rez))                     ;
            MN  = size(originalBlok)                 ;
            if isnan(rez)%akoje NaN stavi da bude nula
               rez = 0;
            end
            rez = rez / (MN(1)*MN(2));
         end
         %% funkcija sinuta sa
         %% http://www.oneder.de/2006/12/05/studium/vorlesungen/medical-image-processing-defect-pixel-interpolation-in-matlab/
         function f = interpdefectimage(g, w, maxit, pad, debug)
        % INTERPDEFECTIMAGE Defect pixel interpolation.
        %
        %   f = INTERPDEFECTIMAGE(g, w, maxit, pad, debug) calculates
        %   the interpolated image f from the observed image g and
        %   the binary defect pixel map w. In w a zero-entry corresponds
        %   to a defect pixel and a 1 to a non-defect pixel. f,g and
        %   w are 2d-images of the same dimensions.
        %   The interpolation is stopped after maxit iterations.
        %   To increase the resolution in the frequency domain, 
        %   the image g and mask w are padded at each side with pad zeros.
        %   The debug parameter controls whether
        %   to show a visualization of the iterative process.

            g = padarray(g,[pad pad]);
            w = padarray(w,[pad pad],1);
            % precalculate useful things
            dim = size(g);            % image dimensions
            halfDim = ceil(dim / 2);  % half dimensions
            G = fft2(g);    % fourier transform of g
            W = fft2(w);    % fourier transform of w
            maxDeltaE_G_Ratio = Inf;
            maxDeltaE_G_Ratio_Tres = 1.0e-6;

            % initialization
            Fhat = zeros(dim);
            FhatNext = zeros(dim);

            % debug
            if (debug~=0)
                figure();
            end

            % iteration
            for i=1:maxit

                % check convergence criterion: stop
                % if the ratio of energy reduction
                % fell below the theshold.
                if (maxDeltaE_G_Ratio <= maxDeltaE_G_Ratio_Tres)
                    % disabled, so just unccoment the break
                    % for enabling the stopping criterion
                    %break;
                end

                % in the i-th iteration select the line pair s1,t1
                % which maximizes the energy reduction
                deltaE_G = abs(G(1:dim(1), 1:halfDim(2)+1));
                [maxDeltaE_G idx] = max(deltaE_G(:));
                idx = find(deltaE_G==maxDeltaE_G);
                r = randperm(numel(idx));
                [s1 t1] = ind2sub(dim, idx(r(1)));

                % calculate the ratio of energy reduction
                % in comparison to the last iteration
                if (i > 1) 
                    maxDeltaE_G_Ratio = abs( (maxDeltaE_G - energyHist(i-1)) / maxDeltaE_G );
                end

                energyHist(i) = maxDeltaE_G; % save the last energy reduction value

                s1f = s1 - 1; % Original row position in frequency domain
                t1f = t1 - 1; % Original column position in frequency domain

                % compute the corresponding linepair s2,t2:
                % mirror the positions at halfDim
                s2 = s1; t2 = t1;
                if s1f > 0
                    s2 = -1*(s1-halfDim(1)) + halfDim(1)+2;
                end
                if t1f > 0
                    t2 = -1*(t1-halfDim(2)) + halfDim(2)+2;
                end
                if s1f==0 && t1f==0 % special case
                    s2 = -1; t2 = -1;
                end

                % compute twice_s1 = mod(2*s1,dim(1)) and tiwce_t1 =
                % mod(2*t1,dim(2)), as we require them in the next step
                twice_s1 = mod(2*s1f,dim(1))+1;
                twice_t1 = mod(2*t1f,dim(2))+1;

                % estimate the new Fhat (FhatNext)
                specialCases = [0 0; 0 halfDim(2); halfDim(1) 0; halfDim(1) halfDim(2)];
                if any(specialCases == repmat([s1f t1f], 4, 1))
                    % handle any of the special cases: (0,0), (0,M/2), (N/2,0),
                    % (N/2,M/2)
                    FhatNext(s1, t1) = FhatNext(s1, t1) + (dim(1)*dim(2))*G(s1,t1)/W(1,1);
                else 
                    % handle the general case
                    tval = dim(1)*dim(2)*( G(s1,t1)*W(1,1) - conj(G(s1,t1))*W(twice_s1,twice_t1));
                    tval = tval / (abs(W(1,1))^2 - abs(W(twice_s1,twice_t1))^2);
                    % Accumulation: Update pair (s1,t1),(s2,t2)
                    FhatNext(s1, t1) = FhatNext(s1, t1) + tval;
                    FhatNext(s2, t2) = FhatNext(s2, t2) + conj(tval);
                end

                % End iteration step by forming the new error spectrum
                G = G - arsIMG.convhelper(FhatNext-Fhat, W, s1, t1);
                % make sure we don't get any rounding errors,
                % G(s1,t1) and G(s1,t2) should be zero
                G(s1,t1) = 0;
                if all([s2,t2]~=-1) G(s2,t2) = 0; end
                Fhat = FhatNext;

                % debug
                if (debug~= 0 && mod(i,20)==0)            
                    fhat = real(ifft2(Fhat)); % obtain the estimation of the real image
                    idx = find(w==0);         % find the important mask entries
                    f = g;
                    f(idx) = fhat(idx);
                    f = f(33:dim(1)-32,33:dim(2)-32); 
                    subplot(2,2,1);
                    imshow(f,[]);
                    title('Current interpolation result');

                    subplot(2,2,2);
                    imshow(fhat,[]);
                    title('Current approximation of f');

                    subplot(2,2,3);
                    plot([1:i], log(energyHist(1:i)/dim(1) +(abs(energyHist(1:i))==0)));
                    title('Error energy reduction');

                    subplot(2,2,4);
                    imshow(abs(G),[]);
                    title('Error spectrum G');

                    drawnow;
                end
            end

            fhat = real(ifft2(Fhat)); % obtain the estimation of the real image
            idx = find(w==0);         % find the important mask entries
            f = g;                    % return the result
            f(idx) = fhat(idx);
            f = f(pad+1:dim(1)-pad,pad+1:dim(2)-pad); 
        end


        function G = convhelper(F, W, s, t)
        % do the convolution of the m-times-n matrix F and W
        % s,t is the position of the selected
        % line pair, the convolution is simplified
        % in the following way:
        % G(k1,k2) = F(k1,k2) 'conv' W(k1,k2) 
        %          = (F(s,t)W(k1-s,k2-t) + F*(s,t)W(k1+s,k2+t)) / (MN)
        % where F* is the conjugate complex.

            sz = size(F);
            G = zeros(sz);
            sf = s-1;
            tf = t-1;
            F_st = F(s,t);
            F_stc = conj(F_st);

            [X Y] = meshgrid(0:sz(1)-1, 0:sz(2)-1);
            I = sub2ind(sz, X(:)+1, Y(:)+1);
            I_neg = sub2ind(sz, mod(X(:)-sf,sz(1))+1, mod(Y(:)-tf,sz(2))+1);
            I_pos = sub2ind(sz, mod(X(:)+sf,sz(1))+1, mod(Y(:)+tf,sz(2))+1);

            if (all(I_neg==I_pos))
                % special case
                G(I) = F_st.*W(I_neg) ./ prod(sz);
            else
                G(I) = (F_st.*W(I_neg) + F_stc.*W(I_pos)) ./ prod(sz);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% funkcija brise "prave" - horizontalne ili vertikalne edgeove duze od dt
        % INPUTS
            % - img - binarna slika
        % OUTPTUS
            % - rez - binarna slika bez nepotrebnih edgeova
        function [rez] = izbrisiPredugackePraveEdgeoveZaIVUS(img,ugaoNaKomeJeShum)
             debugMode = 1;
             bb  = edge(img,'canny');
             %%%IZBRISI SHUM
             bb = arsIVUS.izbrisiDeoSlikeNaUglu(bb, ugaoNaKomeJeShum, 10);
             %%%
             %%%
             BW3 = bwmorph(bb,'skel',inf)        ;
             CC  = bwconncomp(BW3,4)  ;
             rez = zeros(CC.ImageSize);
             %IZBRSI PRAVE LINIJE
             for i = 1:CC.NumObjects
%                  i
%                      numel(CC.PixelIdxList{i})
                 if(numel(CC.PixelIdxList{i})<15)                     
                     rez(CC.PixelIdxList{i}) = 1;
                 end
             end
             %IZBRISI KRATKE LINIJE
             rez = arsIMG.izbrisiKratkeDeloveBinarneSlike(rez,14);
             if debugMode
                imshow(rez);
             end
             se  = strel('diamond',3);
             bw2 = imdilate(rez,se)     ;
             if debugMode
                imshow(bw2);
             end
             BW3 = bwmorph(bw2,'skel',inf)        ;
             CC  = bwconncomp(BW3,8)   ;
             rez = BW3;
             if debugMode
                 figure; imshow(rez);
             end
             rez = arsIMG.izbrisiKratkeDeloveBinarneSlike(rez,50);
             CC  = bwconncomp(BW3,4)   ;
             pomRez = zeros(CC.ImageSize);
             pp = pomRez;
             %IZBRISI POZADINSKE LINIJE - ostavi unutrasnje
%              for ugao = 1:360
%                  pom =0;
%                  for precnik = 1:CC.ImageSize(2)/2
%                      if(pom)%ako je nasao konturu na uglu-predji na sledeci
%                          continue;
%                      end
%                      koordinate = arsIMG.polarToDekartKoordinate([precnik,ugao],round(CC.ImageSize/2) );
%                      koordinate = round(koordinate);
%                      i = koordinate(1);
%                      j = koordinate(2);
%                      pp(i,j)=1;
% %                      f = figure;
% %                      imshow(pp);
% %                      close(f);
%                      %nadji prvi pixel razlicit od nule
%                      if (rez(i,j)==0)
%                          pP = arsIMG.daLiSeNalziUOkoliniTacke([i,j],rez);
%                          i = pP.koo(1);
%                          j = pP.koo(2);
%                      end                         
%                      if(rez(i,j)==1)
%                          %nadji konturu kojoj taj pixel pripada
%                          idPixela = arsIMG.pixelsKoordinate2IDS(rez,[i,j]);
%                          for brKonture = 1:CC.NumObjects
%                              daLiPripadaKonturi = find(CC.PixelIdxList{brKonture}==idPixela);
%                              %ako je nasao konturu kojoj pripada
%                              if ~isempty(daLiPripadaKonturi)
%                                  %oboji tu tacku
%                                  pomRez(CC.PixelIdxList{brKonture})=1;
% %                                  imshow(pomRez);
%                                  %i izadji=predji na novi ugao
%                                  pom = 1;%
%                                  break;
%                              end
%                          end
%                      end
%                  end
%              end
               if debugMode
                    figure; imshow(rez);
               end
%              se  = strel('diamond',4);
%              bw2 = imdilate(rez,se)  ;
%              bw2 = im2uint8(bw2); 
%              bw2 = arsIMG.gaussianSmoothing(bw2, [ 10 10], 2);
%              imshow(bw2);
%              bw2 = arsIMG.anizotropnaDifuzija(bw2);
%              imshow(bw2);
        end
        % da li se nalazi u okolini tacke
        function [rez] = daLiSeNalziUOkoliniTacke(tacka,img)
            rez.is  = 0;
            rez.koo = tacka;
            dt  = 2;
            for i = tacka(1)-dt:tacka(1)+dt
                for j = tacka(2)-dt:tacka(2)+dt
                    if img(i,j)==1
                        rez.is = 1;
                        rez.koo= [ i, j];
                        return;
                    end
                end
            end
        end
        %% 
        function [rez] = anizotropnaDifuzija(img)            
            % anizotropna difuzija
            % function y = nldifc( u, lambda, sigma, m, stepsize, nosteps, varargin)
          brKorakaAnizotropneDifuzije = 1;%4
          rez = nldifc(img,...
                linspace(1.5,4.0,brKorakaAnizotropneDifuzije),...lambda  1.5,4
                linspace(3.0,6.3,brKorakaAnizotropneDifuzije),...sigma 3,0.3
                10,...m 10
                1000,...stepsize 1000
                brKorakaAnizotropneDifuzije,...nosteps
                1,2,'aos', 'grad', 'dfstep', 2, 'alt1','imscale', 'norm', 1);%varargin
        end
        %%  konvertuje jpgove iz DICOMA u mat fajlove koji se brze ucitavaju
        function konvertujJpg2Mat(osnovaSlike)
            dcm = dicominfo(osnovaSlike);
            for i=1:dcm.NumberOfFrames
                img1Path = [osnovaSlike '_'  num2str(i,'%.4i') '.jpg']    ;
                img      = imread(img1Path)                               ;
                savefile = [osnovaSlike '_'  num2str(i,'%.4i') '.mat']    ;
                save(savefile, 'img');
            end
        end
        %%
        function [rez] = izbrisiKratkeDeloveBinarneSlike(img, minDuzina)
             BW3 = img;
             CC  = bwconncomp(BW3,8)   ;
             rez = BW3;
             for i = 1:CC.NumObjects
                 numel(CC.PixelIdxList{i})
                 if(numel(CC.PixelIdxList{i})<minDuzina)       
                     koordinate = arsIMG.pxelIDS2Koordinate(rez,CC.PixelIdxList{i});
                     for j =1:numel(koordinate(:,1))
                         rez(koordinate(j,1),koordinate(j,2))=0;
                     end
                     
                 end
             end
        end
        %% 
        function [gx,gy]=gaussgradient(IM,sigma)
        %GAUSSGRADIENT Gradient using first order derivative of Gaussian.
        %  [gx,gy]=gaussgradient(IM,sigma) outputs the gradient image gx and gy of
        %  image IM using a 2-D Gaussian kernel. Sigma is the standard deviation of
        %  this kernel along both directions.
        %
        %  Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn)
        %  at Tsinghua University, Beijing, China.

        %determine the appropriate size of kernel. The smaller epsilon, the larger
        %size.
            epsilon=1e-2;
            halfsize=ceil(sigma*sqrt(-2*log(sqrt(2*pi)*sigma*epsilon)));
            size=2*halfsize+1;
            %generate a 2-D Gaussian kernel along x direction
            for i=1:size
                for j=1:size
                    u=[i-halfsize-1 j-halfsize-1];
                    hx(i,j)=arsIMG.gauss(u(1),sigma)*arsIMG.dgauss(u(2),sigma);
                end
            end
            hx=hx/sqrt(sum(sum(abs(hx).*abs(hx))));
            %generate a 2-D Gaussian kernel along y direction
            hy=hx';
            %2-D filtering
            gx=imfilter(IM,hx,'replicate','conv');
            gy=imfilter(IM,hy,'replicate','conv');
        end        
        %pom. funkcija za gaussgradient
        function y = dgauss(x,sigma)
            %first order derivative of Gaussian
            y = -x * arsIMG.gauss(x,sigma) / sigma^2;
        end
        %%%%%%%%%%%%%%%%%%%%%
        %% selektuje piksele unutar poligona definisang sa XY
        %INPUTS
            %img - slika
            %XY  - koordinate poligona
        %OTUPUTS
            %idPixela   - id piksela unutar region
            %rezPikseli - vrednosti piksela unutar regiona
        function [rezPikseli idPiksela rezMaska] = selektujPikseleUnutarKonture(img,XY)
            [m n]      = size(img);
            rezMaska   = poly2mask(XY(:,1), XY(:,2), m, n);
            idPiksela  = find(rezMaska>0);
            rezPikseli = img(idPiksela);
        end
        %%
        function rezIMG = cropIMG(img, crop)
            rezIMG = img(crop(1,1):crop(2,1), crop(1,2):crop(2,2)); 
        end
        %% Izracunaj shape context
        %INPUT
            %c            - kontura
        %OUTPUTS
            %shapeContext - shape context matrica
        function shapeContext = IzracunajShapeContextZaKonturuC(c)
            for iTacke  = 1 : numel(c(:,1))
                    trenutnaTacka = c(iTacke,:);
                    %postavi lokalni koo sistem u i-toj tacki
                    pomC          = [c(:,1)-trenutnaTacka(:,1), c(:,2)-trenutnaTacka(:,2)];
                    trenutnaTacka = round(trenutnaTacka);
                    %prebaci iz cartesian u log polar
                    [teta r]     = cart2pol(pomC(:,1), pomC(:,2));
                    teta         = radtodeg(teta)                ;% prebaci u stepene
                    teta(teta<0) = 360 + teta(teta<0)            ;% od 0 do 360
                    logr         = log(r)                   ;% r-osu prebaci u log
                %     logr         = log10(r);
                    % Sledecom linijom se izbegava loopovanje pri racunanju shape contexta.
                    % log-polar koordinate se prebacuju u koordinate shape-context regija 
                    pomLogRTeta= [int8(fix(logr(:)))+int8(mod(logr(:),1)>0), int8(teta(:)/nBinsTeta)+int8(mod(teta(:),nBinsTeta)>0)+1]; 
                    %Sada ga samo napakuj u matricu-histogram
                    shapeContext(:,:,iTacke) = zeros(bBinsR,nBinsTeta); % bBinsR x nBinsTeta
                    for i = 1 : numel(pomLogRTeta(:,1))
                        try
                        shapeContext(pomLogRTeta(i,1), pomLogRTeta(i,2), iTacke) = shapeContext(pomLogRTeta(i,1), pomLogRTeta(i,2), iTacke) + 1;
                        catch
                            bb=1;
                        end
                    end
            end
        end
    end%methods
end%bezier class deff