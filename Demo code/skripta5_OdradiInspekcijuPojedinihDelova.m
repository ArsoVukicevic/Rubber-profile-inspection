%%
%INPUTS
    %PROFIL KOJI SE ISPITUJE
    %koordinateBifurkacija - bifurcation points na profilu koji se ispituje
    %konturaProfila        - detektovana kontura profila koji se ispituje
    %rezGTnoveKoordinate   - REFERENTNI PROFIL u koordinatama novog profila
    %img                   - slika trenutnog profila
function  rez = skripta5_OdradiInspekcijuPojedinihDelova(tackeBifurkacija, centralneLinije, kontura,.... TRENUTNI PROFIL
                                                   rezGTnoveKoordinate,........................... GT profil
                                                   img,...
                                                   idPozicije,....
                                                   axesHandle,....
                                                   debugMode,...
                                                   trenutniPolozajMisa)
                                           
   if isstr(img)
       img = imread(img);
       img = img(:,:,1);
   end
   if ~exist('debugMode')
       debugMode = 1;
   end
   img = mat2gray(img);
   % Nadji corresponding tacke
   shapeContextBifurkacijaTrenutnogProfila = nadjiShapeContextZaTacke(tackeBifurkacija, centralneLinije, kontura);
   %% Tacke koje se ispituju   
   TackaKojaSeTrazi = [{[rezGTnoveKoordinate.TackeZaMatching.LevoA  ; rezGTnoveKoordinate.TackeZaMatching.LevoB ]},... Levo spoljasnje perce
                       {[rezGTnoveKoordinate.TackeZaMatching.LevoB  ; rezGTnoveKoordinate.TackeZaMatching.LevoC ]},... Levo horizontalno krilce
                       {[rezGTnoveKoordinate.TackeZaMatching.LevoD  ; rezGTnoveKoordinate.TackeZaMatching.LevoE ]},... Levo unutrasnje perce
                       {[rezGTnoveKoordinate.TackeZaMatching.LevoD  ; rezGTnoveKoordinate.TackeZaMatching.LevoF ]},... Levi vrat
                       {[rezGTnoveKoordinate.TackeZaMatching.LevoG  ; rezGTnoveKoordinate.TackeZaMatching.DesnoG]},... Baza
                       {[rezGTnoveKoordinate.TackeZaMatching.DesnoD ; rezGTnoveKoordinate.TackeZaMatching.DesnoF]},... Desni vrat
                       {[rezGTnoveKoordinate.TackeZaMatching.DesnoD ; rezGTnoveKoordinate.TackeZaMatching.DesnoE]},... Desno unutrasnje perce
                       {[rezGTnoveKoordinate.TackeZaMatching.DesnoB ; rezGTnoveKoordinate.TackeZaMatching.DesnoC]},... Desno horizontalno krilce
                       {[rezGTnoveKoordinate.TackeZaMatching.DesnoA ; rezGTnoveKoordinate.TackeZaMatching.DesnoB]}.... Desno spoljasnje perce
       ];
   if exist('trenutniPolozajMisa')
       for i = 1:numel(TackaKojaSeTrazi)
           rastojanje(i) = norm(mean(TackaKojaSeTrazi{i})-trenutniPolozajMisa);
       end
       [bzv idPozicije] = min(rastojanje);
   end
   if exist('idPozicije')
       TackaKojaSeTrazi={TackaKojaSeTrazi{idPozicije}};
   end
   for i = 1 : numel(TackaKojaSeTrazi)  
    idMatchingBifurkacije                   = nadjiMatchingTacku(shapeContextBifurkacijaTrenutnogProfila,......................................................................................... TRENUTNI PROFIL
                                                                 TackaKojaSeTrazi{i}, rezGTnoveKoordinate.koordinateCentralnihLinija, rezGTnoveKoordinate.konturaProfila);%... GT
%     debugMode=1
    if debugMode
        if ~exist('axesHandle')
            axesHandle=figure;
        end
        axes(axesHandle);
%         imshow(mat2gray(img));
        plotLine(kontura); plotLine(rezGTnoveKoordinate.konturaProfila);
        scatter(TackaKojaSeTrazi{i}(:,1), TackaKojaSeTrazi{i}(:,2));
        scatter(tackeBifurkacija(idMatchingBifurkacije,1), tackeBifurkacija(idMatchingBifurkacije,2));
    end
       %% Poravnavanje geometrija
   end
     
    brOkolnihTacki = 166; 
    rez=fitujGeometrijuDelaKojiSePosmatra(tackeBifurkacija(idMatchingBifurkacije,:), centralneLinije                               , kontura,.............................. trenutna kontura 
                                         TackaKojaSeTrazi{i}                      , rezGTnoveKoordinate.koordinateCentralnihLinija, rezGTnoveKoordinate.konturaProfila,... gt
                                         rezGTnoveKoordinate.tolerancije,...
                                         img,....
                                         axesHandle,....
                                         debugMode, brOkolnihTacki);
if nargin > 8  
%     debugMode=1;
    rez=fitujGeometrijuDelaKojiSePosmatra(trenutniPolozajMisa, centralneLinije                               , kontura,.... trenutna kontura 
                                          trenutniPolozajMisa, rezGTnoveKoordinate.koordinateCentralnihLinija, rez.Dicp,... gt
                                          rez.tolerancija,...
                                          img,....
                                          axesHandle,....
                                          debugMode, brOkolnihTacki+33);       
end
    rez.TackaKojaSeTrazi=TackaKojaSeTrazi;  
end



%vrsi se fitovanje 100 tacki oko tackaBifurkacijaGT na 100 tacki oko tackaBifurkacija
%INPUTS
    %
%OUTPUTS
    %rez - linije za iscrtavanje
function rez = fitujGeometrijuDelaKojiSePosmatra(tackaBifurkacija  ,   centralneLinije,   kontura,... trenutna kontura
                                                 tackaBifurkacijaGT, centralneLinijeGT, konturaGT,... gt       kontura
                                                 tolerancija,...
                                                 img,...
                                                 handlesID, debugMode, brOkolnihTacki)
pointCloud     = [kontura  ]; %; centralneLinije   ] ; 
pointCloudGT   = [konturaGT]; %; centralneLinijeGT ] ; 
for i = 1 : numel(tackaBifurkacija(:,1))
    pointCloud     = selektujNajblizihNTacki(pointCloud  , tackaBifurkacija  , brOkolnihTacki);
    pointCloudGT   = selektujNajblizihNTacki(pointCloudGT, tackaBifurkacijaGT, brOkolnihTacki);
end
    %% Run ICP (standard settings)
    M                = pointCloud'   ;
    D                = pointCloudGT' ;
    [Ricp Ticp ER t] = icp(M, D, 30) ;

    % Transform data-matrix using ICP result
    Dicp = Ricp * konturaGT' + repmat(Ticp, 1, numel(konturaGT(:,1)));
    
    for i = 1 : numel(tolerancija) %br toleracija
        tolerancija{i} = Ricp * tolerancija{i}' + repmat(Ticp, 1, numel(tolerancija{i}(:,1)));
        tolerancija{i} = tolerancija{i}';
    end
    debugMode=1
    if debugMode
        if ~exist('handlesID')
            handlesID=figure;
        end
        axes(handlesID);
        cla(handlesID,'reset');
        imshow(mat2gray(img));    
        plotLine(kontura);
        plotLine(Dicp', 2);
        scatter(tackaBifurkacija(:,1), tackaBifurkacija(:,2));
%         scatter(pointCloud(:,1)  , pointCloud(:,2)  );
%         scatter(pointCloudGT(:,1), pointCloudGT(:,2));
        for i = 1 : numel(tolerancija)
            plotLine(tolerancija{i});
        end
    end
    rez.tolerancija      = tolerancija     ;
    rez.kontura          = kontura         ; 
    rez.Dicp             = Dicp'           ;
    rez.tackaBifurkacija = tackaBifurkacija;
    rez.pointCloud       = pointCloud      ;
end

function rez = selektujNajblizihNTacki(pointCloud, tackaBifurkacije, brOkolnihTacki)
    for i = 1 : numel(tackaBifurkacije(:,1))
        d          = AngioIvusMath.arsNorm(AngioIvusMath.arsMinus(pointCloud,tackaBifurkacije(i,:)));
        [bzv id]   = sort(d);
        if i == 1
            rez        = pointCloud(id(1:brOkolnihTacki),:); %uzmi 100 najblizih tacki
        else
            rez = [rez;pointCloud(id(1:brOkolnihTacki),:)];
        end                
    end
end


%INPUTS
    %GT PROFIL
        % shapeContextBifurkacijaTrenutnogProfila - izracunat prethodno na osnovu: gtTackaBifurkacije, gtCentralneLinije, gtKontura
    %trenutni profil
        %tackeBifurkacija
        %centralneLinije
        %kontura
function [idMinimuma rez] = nadjiMatchingTacku(shapeContextBifurkacijaTrenutnogProfila,....
                                         gtTackaBifurkacije, gtCentralneLinije, gtKontura)   
    for iTacke = 1 : numel(gtTackaBifurkacije(:,1))                                 
        gtShapeContext = nadjiShapeContextZaTacke(gtTackaBifurkacije(iTacke,:), gtCentralneLinije, gtKontura);                       
        %nadji najblizu tacku
        for i = 1:numel(shapeContextBifurkacijaTrenutnogProfila(1,1,:))
    %         pom = shapeContextBifurkacijaTrenutnogProfila(:,:,i) - gtShapeContext;
    %         missmatch(1,i) = var(pom(:));%std(abs(pom(:)));
    %           missmatch(1,i) = hist_cost_2(gtShapeContext, shapeContextBifurkacijaTrenutnogProfila(:,:,i));
    %           missmatch(1,i) = trace(gtShapeContext)-trace(shapeContextBifurkacijaTrenutnogProfila(:,:,i));
              pom = shapeContextBifurkacijaTrenutnogProfila(:,:,i);
              missmatch(1,i) = regression(gtShapeContext(:)', pom(:)');
        end      
        [minimum idMinimuma(iTacke)] = max(missmatch(:));
    end
end

%fja racuna shape context za niz tacki (TackaBifurkacije) 
%INPUTS
    %TackaBifurkacije - tacke za koje se racuna SC
    %CentralneLinije, Kontura - oblak tacki
function shapeContext = nadjiShapeContextZaTacke(TackaBifurkacije, CentralneLinije, Kontura, nBinsTeta,nBinsR)
%
    debugMode = 0;
    rMax = 666;%pixela
    if ~exist('nBinsTeta')
        nBinsTeta = 120;
        nBinsR    = 100;
    end
    c = [Kontura; CentralneLinije]; % Oblak tacki
    for iTacke  = 1 : numel(TackaBifurkacije(:,1))
        trenutnaTacka = TackaBifurkacije(iTacke,:);
        if debugMode
            figure; % Prikaz cele konture
                scatter(TackaBifurkacije(:,1), TackaBifurkacije(:,2));
                scatter(CentralneLinije(:,1) , CentralneLinije(:,2) );
                plotLine(Kontura);
            figure; % Prikaz jedne tacke
            plotLine([0 0 0; trenutnaTacka]);
            plotLine(c);
        end
        %postavi lokalni koo sistem u i-toj tacki
        pomC          = [c(:,1)-trenutnaTacka(:,1), c(:,2)-trenutnaTacka(:,2)];
        trenutnaTacka = round(trenutnaTacka);
        %prebaci iz cartesian u log polar
        [teta r]      = cart2pol(pomC(:,1), pomC(:,2));
        id = find(r<rMax);
        r             = r(id);  
        teta          = teta(id);
        teta          = radtodeg(teta)                ;% prebaci u stepene
        teta(teta<0)  = 360 + teta(teta<0)            ;% od 0 do 360
        logr          = r;% u radu je bilo log(r)                        ;% r-osu prebaci u log
        logr(logr<0)  = 0;
    %     logr         = log10(r);
        % Sledecom linijom se izbegava loopovanje pri racunanju shape contexta.
        % log-polar koordinate se prebacuju u koordinate shape-context regija 
        % u radu je logr bilo floor(logr(:)/(log(rMax)/nBinsR))+1
        pomLogRTeta= [floor(logr(:)/((rMax+0.0001)/nBinsR))+1, floor(teta(:)/(360/nBinsTeta))+1]; 
        %Sada ga samo napakuj u matricu-histogram
        shapeContextPom           = zeros(nBinsR,nBinsTeta); % bBinsR x nBinsTeta
        [pomLogRTetaUnique,ia,ic] = unique(pomLogRTeta,'rows');
        pomLogRTetaUnique         = pomLogRTetaUnique(:,1) + (pomLogRTetaUnique(:,2)-1)*nBinsR;
        a_counts                  = accumarray(ic,1);
        shapeContextPom(pomLogRTetaUnique) = a_counts;
%         for i = 2 : nBinsR
% %             shapeContextPom(i,:)=shapeContextPom(i,:) / (6^(i-1));
%         end
        shapeContext(:,:,iTacke) = shapeContextPom;
    end
end

function iscrtajShapeContextUtacki(centar, nBinsR, nBinsTeta, rMax)
centar=centar([2 1]);
for r = log(1):log(rMax)/nBinsR:log(rMax+1)
    viscircles(centar,exp(1)^r);
end    
for th = 0 : (2*pi/nBinsTeta) : 2*pi 
    x = rMax * cos(th) + centar(1);
    y = rMax * sin(th) + centar(2);
    plotLine([x y 0; centar 0]);
end
end

function HC=hist_cost_2(BH1,BH2)
    % HC=hist_cost_2(BH1,BH2);
    %
    % same as hist_cost.m but BH1 and BH2 can be of different lengths

    [nsamp1,nbins]=size(BH1);
    [nsamp2,nbins]=size(BH2);

    BH1n = BH1./repmat(sum(BH1,2)+eps,[1 nbins])       ;
    BH2n = BH2./repmat(sum(BH2,2)+eps,[1 nbins])       ;
    tmp1 = repmat(permute(BH1n,[1 3 2]),[1 nsamp2 1])  ;
    tmp2 = repmat(permute(BH2n',[3 2 1]),[nsamp1 1 1]) ;
    HC   = 0.5*sum(((tmp1-tmp2).^2)./(tmp1+tmp2+eps),3);
end