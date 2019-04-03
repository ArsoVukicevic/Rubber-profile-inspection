%Skripta vrsi preklapanje tehnickog crteza sa trenutnog slikom
%INPTUS
    %koordinateCentar, xPravac,yPravac, scaleFaktor - kalibracija
    %imgPath - Slika sa fotoaparata
    %kontura - za preklapanje
%OUTPUTS - nema, radi se samo vizualizacija    
    %rezKontura - kontura referentnog profila u lokalnim koordinatama
    %rez        - rezGT u lokalnim koordinatama
function [rezKontura rez]= skripta3_PrecklopiDetektovanuSaoriginalnomKonturom(koordinateCentar, xPravac, yPravac, scaleFaktor, imgPath, kontura, debugMode)
 
rezGT = load('include\rezGT.mat')      ;
rezGT = rezGT.rezGT;

rez = konvertujIzReferentnogMMuTrenutniPixelKooSistem(rezGT, koordinateCentar, xPravac, yPravac, scaleFaktor);
rezKontura = rez.konturaProfila;
if debugMode
    imshow(imgPath);
    plotLine(rezKontura);
end
end

function rez = konvertujIzReferentnogMMuTrenutniPixelKooSistem(rezGT, koordinateCentar, xPravac, yPravac, scaleFaktor)
rez.konturaProfila             = pom(rezGT.konturaProfila            , koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.koordinateCentralnihLinija = pom(rezGT.koordinateCentralnihLinija, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
for i = 1 : numel (rezGT.tolerancije )
    rez.tolerancije{i}                = pom(rezGT.tolerancije{i}     , koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
end
%
rez.TackeZaMatching.DesnoA = pom(rezGT.TackeZaMatching.DesnoA, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoB = pom(rezGT.TackeZaMatching.DesnoB, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoC = pom(rezGT.TackeZaMatching.DesnoC, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoD = pom(rezGT.TackeZaMatching.DesnoD, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoE = pom(rezGT.TackeZaMatching.DesnoE, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoF = pom(rezGT.TackeZaMatching.DesnoF, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.DesnoG = pom(rezGT.TackeZaMatching.DesnoG, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
%--
rez.TackeZaMatching.LevoA = pom(rezGT.TackeZaMatching.LevoA, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoB = pom(rezGT.TackeZaMatching.LevoB, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoC = pom(rezGT.TackeZaMatching.LevoC, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoD = pom(rezGT.TackeZaMatching.LevoD, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoE = pom(rezGT.TackeZaMatching.LevoE, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoF = pom(rezGT.TackeZaMatching.LevoF, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
rez.TackeZaMatching.LevoG = pom(rezGT.TackeZaMatching.LevoG, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT);
%prebaci ostale tacke
end

function rezTacka = pom(TackaGT, koordinateCentar, xPravac, yPravac, scaleFaktor, rezGT)
vTranslacije = [-14.8, -33.7, 0]                             ;
rezTacka     = TackaGT * rezGT.scaleFaktor                   ;
rezTacka     = AngioIvusMath.arsPlus(rezTacka, vTranslacije) ;

%prebaci konturu iz mm u pixele
rezTacka = rezTacka/scaleFaktor;
rezTacka = AngioIvusMath.arsPlus(AngioIvusMath.arsPuta(xPravac, rezTacka(:,1)),...
                                 AngioIvusMath.arsPuta(yPravac, rezTacka(:,2))); 
rezTacka =  AngioIvusMath.arsPlus(rezTacka,koordinateCentar);
%prebaci ostale tacke
end
