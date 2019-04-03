classdef AngioIvusMath
    properties
        a=1;
    end
	methods(Static=true)
        %% interpolira 
        %% interpoliranje fejsa/qudrilateriala
        %      D._________________.C
        %       | |y              /
        %       | |              /
        %       | |------>x     /
        %      A.____________./B
        % 
        % INPUTS
            %cvorovi[brCvorova x 3]
            %fejso[brFesojsva x4]
            %X,Y br podela po pravcima x,y
        % OUTPUTS
            %rezTacke 
            %fejsovi
        function [rezTacke fejsovi]= interpolirajQuadrilaterial(cvorovi,fejs, X, Y)
            plotLine(cvorovi);
            A = cvorovi(fejs(:,1),:);
            B = cvorovi(fejs(:,2),:);
            C = cvorovi(fejs(:,3),:);
            D = cvorovi(fejs(:,4),:);
            if(X>1)%radi se podele, x=br podela
                AD = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(A,D, Y);
                BC = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(B,C, Y);
                plotLine(AD);
                plotLine(BC);
                Cvorovi = [ 0 0 0];
                for i =1:numel(AD(:,1))
                    pom = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(AD(i,:),BC(i,:), X);
%                     plotLine(pom);
                    Cvorovi = [ Cvorovi;pom ];
                end
                Cvorovi = Cvorovi(2:end,:);
%                 plotLine(Cvorovi);
                fejs= [ 0 1  X+1 X]
                prviRedFejsova = fejs;
                for i =1:X-1
                    prviRedFejsova(i,:) = fejs + i;
%                     plotLine(Cvorovi(prviRedFejsova(i,[1 2 3 4 1]),:));
                end
                fejsovi=prviRedFejsova;
                for i =1:Y-2
                    fejsovi=[fejsovi;prviRedFejsova+X*i];
                end
                for i =1:numel(fejsovi(:,1))
                    plotLine(Cvorovi(fejsovi(i,[1 2 3 4 1]),:));
                end            
                rezTacke = Cvorovi;
                fejsovi  = fejsovi;
%                 close all;
            else% x=0.2(vrati tacku na 0.2)
                AD       = AngioIvusMath.interpolirajRastojanjeIzmedju2TackeVratiTackuNaT(A,D, Y);
                BC       = AngioIvusMath.interpolirajRastojanjeIzmedju2TackeVratiTackuNaT(B,C, Y);
                rezTacke = AngioIvusMath.interpolirajRastojanjeIzmedju2TackeVratiTackuNaT(B,C, X);
            end            
        end
        %% hexa=[ A B C D E F G H]
        %       H ._______.G
        %        /|      /|
        %     E./_|___F./ |
        %      | D|____|_ /C
        %      | /     | /
        %     A._______./B
        %     
        % modRada = 1 (usitni na 27 elemenata)
        % modRada = 2 (ostavi jednu stranu neusitnjenu)
        function [rezCvorovi Elementi]=remeshHex(Cvorovi,Hexe, modRada)
            rezCvorovi=[0 0 0];
            Elementi = [1:8];
            for i=1:numel(Hexe(:,1))
                Hexa=Hexe(i,:);
                A=Cvorovi(Hexa(:,1),:);
                B=Cvorovi(Hexa(:,2),:);
                C=Cvorovi(Hexa(:,3),:);
                D=Cvorovi(Hexa(:,4),:);
                E=Cvorovi(Hexa(:,5),:);
                F=Cvorovi(Hexa(:,6),:);
                G=Cvorovi(Hexa(:,7),:);
                H=Cvorovi(Hexa(:,8),:);
                %------------*-----------*-------------
                AE = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(A,E,4);%BILO 4 , deli se na 3 dela
                BF = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(B,F,4);
                CG = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(C,G,4);
                DH = AngioIvusMath.interpolirajRastojanjeIzmedju2Tacke(D,H,4);

                [Cvorovi1 Fejsovi1] = AngioIvusMath.interpolirajQuadrilaterial([AE(1,:);BF(1,:);CG(1,:);DH(1,:)],[1 2 3 4], 4, 4);
                [Cvorovi2 Fejsovi2] = AngioIvusMath.interpolirajQuadrilaterial([AE(2,:);BF(2,:);CG(2,:);DH(2,:)],[1 2 3 4], 4, 4);
                [Cvorovi3 Fejsovi3] = AngioIvusMath.interpolirajQuadrilaterial([AE(3,:);BF(3,:);CG(3,:);DH(3,:)],[1 2 3 4], 4, 4);
                [Cvorovi4 Fejsovi4] = AngioIvusMath.interpolirajQuadrilaterial([AE(4,:);BF(4,:);CG(4,:);DH(4,:)],[1 2 3 4], 4, 4);

                if modRada==2
                    Cvorivi34 = (Cvorovi2+Cvorovi3)/2;%sloj izmedju 3 i 4
                    rezCvorovi_i = [Cvorovi1;....prvi sloj
                                  Cvorovi2;....drugi sloj
                                  Cvorovi4(1,:);  Cvorovi3([2 3],:)   ;  Cvorovi4(4 ,:);...
                                  Cvorovi3(5,:);  Cvorivi34([6   7],:);  Cvorovi3(8 ,:);...
                                  Cvorovi3(9,:);  Cvorivi34([10 11],:);  Cvorovi3(12,:);...
                                  Cvorovi4(13,:);  Cvorovi3([14 15],:);  Cvorovi4(16,:);
                                  ];
                    prviElement       = [17 18 2 1   21 22 6 5];%1 2 6 5 17 18 22 21
%                     prviElement       = [17 21 22 18  1 5 6 2];%17 18 22 21   1 2 6 5
                    prviRedElemenata  = [prviElement; prviElement+1; prviElement+2];
                    prviSlojElemenata = [prviRedElemenata; prviRedElemenata+4;prviRedElemenata+8];
                    Elementi_i        = [prviSlojElemenata;...prvi sloj elemenata-samo kocke
                                         prviSlojElemenata+16;....drugi sloj 
                                         %---treci sloj, ostalo jos 5 elementa
                                         33 34 38 37   45 46 42 41;.....sa desne | 33 45 46 34     37 41 42 38  ____ 34 46  45 33   38 42 41 37
                                         35 36 40 39   47 48 44 43;.....sa leve  | 36 48  47 35    40 44 43 39
                                         34 35 39 38   46 47 43 42;.....piramida ispod drugog sloja a iznad dna 34 35 47 46    38 39 43 42
                                         33 36 35 34   45 48 47 46;.....na dnu   | 33 36  48 45   34 35 47 46
                                         ];
                else
                    rezCvorovi_i = [Cvorovi1;....prvi sloj
                                  Cvorovi2;....drugi sloj
                                  Cvorovi3;...
                                  Cvorovi4;
                                  ];
                    prviElement       = [17 18 2 1   21 22 6 5];
                    prviRedElemenata  = [prviElement; prviElement+1; prviElement+2];
                    prviSlojElemenata = [prviRedElemenata; prviRedElemenata+4;prviRedElemenata+8];
                    Elementi_i        = [prviSlojElemenata;...prvi sloj elemenata-samo kocke
                                         prviSlojElemenata+16;....drugi sloj .
                                         prviSlojElemenata+32];
                end
                Elementi_i=Elementi_i+(i-1)*numel(rezCvorovi_i(:,1));
                rezCvorovi=[rezCvorovi;rezCvorovi_i];
                Elementi=[Elementi  ;Elementi_i];
            end
            rezCvorovi = rezCvorovi(2:end,:);
            Elementi   = Elementi(2:end,:);
        end
        %%
        function [rez]=zapreminaTetraedarskogModela(Cvorovi,Tetre)
            a = Cvorovi(Tetre(:,1),:);
            b = Cvorovi(Tetre(:,2),:);
            c = Cvorovi(Tetre(:,3),:);
            d = Cvorovi(Tetre(:,4),:);
            ZapreminaTetriABCD  = AngioIvusMath.arsDot((a-b),AngioIvusMath.arsCross(b-d,c-d))/6;
            rez = sum(abs(ZapreminaTetriABCD));
        end
        %% racuna zapreminu FEM modela / mreze sastavljene od Heksaedara
        %INPUTS
            %Hexe[ A[Ax Ay Az], B[Bx By Bz], C[Cx Cy Cz]]
        function[rez]=zapreminaHexaedronskogModela(Cvorovi,Hexe)
            if nargin==1%prosledjen je PAKC/PAKF objekat
                Hexe    = Cvorovi.settings.C8.Elementi;
                Cvorovi = Cvorovi.settings.C7.Cvorovi;
            end
            %svaka heksa se sastoji iz 4 tetre
            tetreSelektor = [1 2 3 6; 1 3 4 6; 1 4 5 6; 5 4 8 6; 4 8 6 3; 8 6 7 3; ];
            %obrazac za zapreminu tetre 
            %ref-http://stackoverflow.com/questions/9866452/calculate-volume-of-any-tetrahedron-given-4-points
            %%%%-http://en.wikipedia.org/wiki/Tetrahedron
            for i = 1:numel(tetreSelektor(:,1))
                trenutnaHexa = Hexe(:,tetreSelektor(i,:));
                if i==1
                    CvoroviABCD  = trenutnaHexa;
                else
                    CvoroviABCD  = [CvoroviABCD; trenutnaHexa];
                end
%                 plotLine(Cvorovi(trenutnaHexa,:));
%                 plotLine(Cvorovi([trenutnaHexa(1) trenutnaHexa(3)],:));
%                 plotLine(Cvorovi([trenutnaHexa(1) trenutnaHexa(4)],:));
%                 plotLine(Cvorovi([trenutnaHexa(4) trenutnaHexa(2)],:));
            end
            a = Cvorovi(CvoroviABCD(:,1),:);
            b = Cvorovi(CvoroviABCD(:,2),:);
            c = Cvorovi(CvoroviABCD(:,3),:);
            d = Cvorovi(CvoroviABCD(:,4),:);
            ZapreminaTetriABCD  = AngioIvusMath.arsDot((a-b),AngioIvusMath.arsCross(b-d,c-d))/6;
            rez = sum(abs(ZapreminaTetriABCD));
        end
        %% tacka se postavlja iz XY ravni 2D u XY u 3D(XYZ). Ravan je
        %  definisana sa vektorima XY, Z nije potreban
        function [rez] = postavi2DtackuU3D(Tacka, X,Y, O)
            X1 = AngioIvusMath.arsPuta([Tacka(:,1),Tacka(:,1),Tacka(:,1)],X);           
            Y1 = AngioIvusMath.arsPuta([Tacka(:,2),Tacka(:,2),Tacka(:,2)],Y); 
            rez=X1+AngioIvusMath.arsPlus(X1+Y1, O);
        end
        %% Postavi 2D konturu na 3D bifukracija patch (koji se sastvoji iz 2 polovine 2 razlicita patcha)
        %           Y1 ^ Y2
        %              |
        %   X1 <-------+-------->  X2
        %              |
        % X1 ==  Normala, Y1 == Binormala globalnog koo sistema patcha(BSplajna)
        %INPUTS
            % Tacka - 2D konture tacki koje se postavljaju. Tacke su poredjane u   y|_x koordinatnom sistemu. Prva tacka je na uglu 0 (1,)
            % X1,Y1 - orijentacija prve polovine patcha  (X1Y1)
        function rezTacka3D = postavi2DtackuU3D_BifurkacijaPatch(Tacka2D, X1, Y1, O, X2, Y2)
            for i = 1:numel(Tacka2D(:,1))
                if Tacka2D(i,1)>0 %odaberi X osu
                    X = X1;
                else
                    X = X2;
                end
                if Tacka2D(i,2)>0 %odaberi Y osu
                    Y = Y1;
                else
                    Y = Y2;
                end
                rezTacka3D(i,:) = O + X*Tacka2D(i,1) + Y*Tacka2D(i,2);
            end
        end
        %% INPUTS
            % T - tacka kojoj prva tacka treba biti najbliza
            % kontura - kontura sa tackama
            % Centar - centar kontura
            % Noormala na ravan
        % OUTPUTS
        function [rez]=sortiraj3DkonturuUodnosuNaTackuT(kontura,T,Centar,NormalaNaRavan)
            NormalaNaRavan=NormalaNaRavan/norm(NormalaNaRavan);
            %nadji tacku iz konture najblizu tacki T
           dmin=6666666;
           for i =1:numel(kontura(:,1,:))
               if norm(T-kontura(i,:,:))<dmin
                   dmin=norm(T-kontura(i,:,:));
                   prvaTacka = kontura(i,:,:);
               end
           end
           % osnovni ugao je prvaTacka,Centar. u odnosu na njega se mere uglovi
           uglovi = kontura(:,1,:)      ;
           Vo     = prvaTacka-Centar    ;
           Vo     = Vo/norm(Vo)         ;
           for i = 1:numel(uglovi)
               %vektor trenutna tacka-centar
               Vi       = kontura(i,:,:)-Centar               ;
               Vi       = Vi/norm(Vi)                         ;
               %nadji ugao izmedju Vo i Vi
               ugao     = acos(dot(Vo,Vi)/(norm(Vo)*norm(Vi)));
               %provera da li je pozitivan ili negativan ugao
               pomV=cross(Vo,Vi);
               pomV = pomV/norm(pomV);
               if ugao ~= 0%ne poklapaju se tacke
                   dist=norm(pomV-NormalaNaRavan);
                   if dist>1
                       ugao=2*pi-ugao;
                   end
                   uglovi(i) = ugao;
               else
                   uglovi(i) = 0;
               end%if
           end
           sort(uglovi);
           p = sortrows(unique([uglovi,kontura],'rows'));
           rez=p(:,2:end,:);
%            x=rez(:,1,:);y=rez(:,2,:);z=rez(:,3,:);
%            line(x,y,z);
        end
        %% rotiraj 2D tacku/ke oko tacke za ugao alfa
        %ref:http://answers.yahoo.com/question/index?qid=20120411030225AAZ7BZw
        function [rez] = rotiraj2DTackuAOkoTackeTZaUgaoAlfa(A, T, alfa)
            debugMode=0;
            brTacaka=numel(A(:,1));
            if debugMode
                plotLine(A); hold on; plotLine([0 0], T);
            end
            Z=A(:,3,:);
            A =[A(:,1,:) A(:,2,:)];
            pomT = A;
            pomT(:,1)=T(1);
            pomT(:,2)=T(2);
            A = A - pomT;            
            matricaRotacije = [cos(deg2rad(alfa)), -sin(deg2rad(alfa));...
                               sin(deg2rad(alfa)),  cos(deg2rad(alfa))];
            for i = 1:brTacaka                
                A(i,:)=matricaRotacije*A(i,:)';
            end
            rez= [A + pomT, Z];
            if debugMode
                hold on; plotLine(rez);
            end
        end
        %% rotiraj 3D tacku oko linije 
        %P1[a,b,c] P2[d,e,f]
        function [rez]=rotiraj3DTackuTOkoLinijeLPoP1(T,P1,P2,teta)
            quat = arsQuanternion;
            quat = quat.setRealniImaginarniDeo(quat,teta,P2-P1);
%             if numel(T(:,1))>1
%                 pomP1=T;
%                 pomP1(:,1)=P1(1);pomP1(:,2)=P1(2);pomP1(:,3)=P1(3);
%             else
%                 pomP1 = P1;
%             end
%             V    = T-pomP1;
            V    = T-P1;
            rez = P1 + quat.rotirajVektorOkoQuanterniona(V,quat);
            %
%             plotLine(T,P1);plotLine(rez,P2);
%             a=P1(1);b=P1(2);c=P1(3);
%             d=P2(1);e=P2(2);f=P2(3);            
%             pravacLinije=P2-P1;
%             u=pravacLinije(1);v=pravacLinije(2);w=pravacLinije(3); 
%             L = sum(pravacLinije .* pravacLinije);
%             cosTeta=cos(rad2deg(teta));
%             sinTeta=sin(rad2deg(teta));
%             sqrtL=sqrt(L);
%             matricaRotacije=...
%             [ (u*u+(v*v+w*w)*cosTeta),...
%              (u*v*(1-cosTeta)-w*sqrtL*sinTeta),...
%              (u*w*(1-cosTeta)+v*sqrtL*sinTeta),...
%              ((a*(v*v+w+w)-u*(b*v+c*w))*(1-cosTeta)+(b*w-c*w)*sqrtL*sinTeta);...
%              ...DRUI RED...
%              (u*w*(1-cosTeta)+w*sqrtL*sinTeta),...
%              (v*v+(u*u+w*w)*cosTeta),...
%              (u*w*(1-cosTeta)-u*sqrtL*sinTeta),...
%              ((b*(v*v+w*w)-v*(a*u+c*w))*(1-cosTeta)+(c*u-a*w)*sqrtL*sinTeta);...
%              ...TRECI RED...
%              (u*w(1-cosTeta)-v*sqrtL*sinTeta),...
%              (u*w*(1-cosTeta)+u*sqrtL*sinTeta),...
%              (w*w+(u*u+v*v)*cosTeta),...
%              ((c*(u*u+v*v)-w*(a*u+b*v))*(1-cosTeta)+(a*v-b*u)*sqrtL*sinTeta);...
%              0,0,0,1];
%          rez=matricaRotacije * [T(:);1];            
        end
        %% rotira tacku(e) oko druge tacke, vektor rotacije je
        % INPUTS
            %TackeZaRotaciju
            %TackaOkoKojeSeRotira
            %alfa - rotacija oko Xose  --- STEPENi - konvertuje se u radijane!!!---
            %beta - rotacija oko Yose
            %gama - rotacija oko Zose
        function [rez] = rotiraj3DTackuTOkoTackeA(TackeZaRotaciju, TackaOkoKojeSeRotira, alfa, beta, gama)
            %1 od TackeZaRotaciju oduzmi TackaOkoKojeSeRotira
            pomTackeZaRotaciju = AngioIvusMath.arsMinus(TackeZaRotaciju,TackaOkoKojeSeRotira);
            M = makehgtform('xrotate',deg2rad(alfa),'yrotate',deg2rad(beta),'zrotate',deg2rad(gama));
%             M1 = rotx(alfa)*roty(beta)*rotz(gama);
            for i = 1:numel(pomTackeZaRotaciju(:,1))
                pom              = M * [ pomTackeZaRotaciju(i,:)' ; 1]                    ;
                pomRotirano      = pom(1:3)'                                              ;            
                rez(i,:)         = AngioIvusMath.arsPlus(pomRotirano,TackaOkoKojeSeRotira);            
            end       
            %
% %             f = figure; hold on;
% %             plotLine(pomRotirano);
% %             plotLine([TackaOkoKojeSeRotira; rez]);
% %             plotLine([TackaOkoKojeSeRotira; TackeZaRotaciju]);
% %             close(f);
        end
        %% funkcija nalazi srednju tacku najblizeg rastojanja dve linije
        %INPUTS - linija L1 i L2
            %linija L1 je definisana sa Po,P1
            %linija L2 je definisana sa Qo,Q1
        %OUTPUTS - najkraca duz koja spaja L1 i L2
            %Psc,Qtc - definise najkracu duz
            %M       - tacka na sredini duzi Psc,Qtc
        function [Psc, Qtc, M, tp, tq] = najblizaTackaIzmedjuLinija(Po, P1, Qo, Q1)   
% % % % % %             u = P1-Po;
% % % % % %             wo = Po-Qo;
% % % % % %             v = Q1-Qo;
% % % % % %             a = dot(u,u);
% % % % % %             b = dot(u,v);
% % % % % %             c = dot(v,v);
% % % % % %             d = dot(u,wo);
% % % % % %             e = dot(v,wo);
% % % % % %             if((a*c-b*b) ~= 0)
% % % % % %               Sc =(b*e-c*d)/(a*c-b*b);
% % % % % %               Tc =(a*e-b*d)/(a*c-b*b);
% % % % % %               Psc = Po+Sc*u;
% % % % % %               Qtc = Qo+Tc*v;   
% % % % % %              else %end if
% % % % % %             %linije su paralelne
% % % % % %             Tc = e/c;
% % % % % %             Qtc = Qo+Tc*(Q1-Qo);  
% % % % % %             Psc = Qtc;
% % % % % %             end
% % % % % %             M = (Qtc + Psc)/2;            
            u = P1-Po;
            wo = Po-Qo;
            v = Q1-Qo;
            a = AngioIvusMath.arsDot(u,u);
            b = AngioIvusMath.arsDot(u,v);
            c = AngioIvusMath.arsDot(v,v);
            d = AngioIvusMath.arsDot(u,wo);
            e = AngioIvusMath.arsDot(v,wo);
            if((a .* c  - b .* b) ~= 0)
              Sc =(b.*e-c.*d)./(a.*c-b.*b);
              Tc =(a.*e-b.*d)./(a.*c-b.*b);
              Psc = Po+[Sc Sc Sc].*u;
              Qtc = Qo+[Tc Tc Tc].*v;   
            else %end if
            %linije su paralelne
            Tc = e./c;
            Qtc = Qo+Tc.*(Q1-Qo);  
            Psc = Qtc;
            end
            M = (Qtc + Psc)/2;
            %parametri(za optimizaciju)
            tp =  AngioIvusMath.arsNorm(Po-M) ./ AngioIvusMath.arsNorm(Po-P1);
            tq =  AngioIvusMath.arsNorm(Qo-M) ./ AngioIvusMath.arsNorm(Qo-Q1);
        end
        %%
         function rez = getTackuPresekaDveLinije(Po, P1, Qo, Q1)
           [Psc, Qtc, M] = najblizaTackaIzmedjuLinija(Po, P1, Qo, Q1);
           if((norm(Psc-Qtc) + norm(Qtc-M))<0.0001)%linije se seku
               rez = M;
           else%linije se ne seku vrati :
                %Psc-Qtc najkraca duz izmedju L1 i L2
                %M  - srednja tacka dui Psc-Qtc
               rez = [Psc, Qtc, M];
           end
         end
         %% nadji presek linije i BSplajna. Za odredjivanje patcheva pri 3d rekonstrukciji
         function rez = getTackuPresekaLinijeIBsplajna(BSplajn1 ,BSplajn2, Po, N)
            debugMode = 0;
            tackeSplajna1 = BSplajn1.interpolateFromAtoB(BSplajn1, 0,0.01,1);
            tackeSplajna2 = BSplajn1.interpolateFromAtoB(BSplajn2, 0,0.01,1);
            if debugMode
                plotLine([Po;Po+N*10]); plotLine([Po;Po-N*10]);
            end
           for i=1:numel(tackeSplajna1(:,1))-1
               if debugMode
                   f=figure;
                   plotLine([Po;Po+N*10]);plotLine([tackeSplajna1(i,:);tackeSplajna1(i+1,:)]);plotLine([tackeSplajna2(i,:);tackeSplajna2(i+1,:)]);
                   close(f);
               end
               T1 = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(Po, Po+N*16, tackeSplajna1(i,:), tackeSplajna1(i+1,:));
%                if isstr(T1)
%                    T1 = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(Po, Po+N-10, tackeSplajna1(i,:), tackeSplajna1(i+1,:));
%                end
               T2 = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(Po, Po-N*16, tackeSplajna2(i,:), tackeSplajna2(i+1,:));
%                if isstr(T1)
%                    T2 = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(Po, Po-N*10, tackeSplajna2(i,:), tackeSplajna2(i+1,:));
%                end
               if ~isstr(T1)
                rez.T1       = T1(1,:);
                rez.duzinaT1 = AngioIvusMath.arsNorm(T1-Po);
               end
               if ~isstr(T2)
                rez.T2       = T2(1,:);
                rez.duzinaT2 = -AngioIvusMath.arsNorm(T2-Po);
               end
           end
            if exist('rez')%promeni orijentaciju i odradi opet
                if rez.duzinaT1 > AngioIvusMath.arsNorm(rez.T1-Po+N)                
                    rez.duzinaT1 = -rez.duzinaT1;
                    rez.duzinaT2 = -rez.duzinaT2;
                end
            else
                rez = AngioIvusMath.getTackuPresekaLinijeIBsplajna(BSplajn1,BSplajn2, Po, -N)
            end
%          if numel(rez(:,1))==1
%              rez='error';%nema preseka
%          else
%              rez=rez(2,:);
% %          end
         end
         %%
         function [rez Psc Qtc tPsc tQtc] = getTackuIzmeldjuDveLinije(Po, P1, Qo, Q1)
           [Psc, Qtc, M] = AngioIvusMath.najblizaTackaIzmedjuLinija(Po, P1, Qo, Q1);
               rez =  Psc + (Qtc - Psc)/2;
         end
         %% rastojanje imedju dve linije
         function rez = getRazstojanjeIzmedjuDveLinije(Po, P1, Qo, Q1)
           [Psc, Qtc, M] = AngioIvusMath.najblizaTackaIzmedjuLinija(Po, P1, Qo, Q1);
               rez = norm(Psc-Qtc);
         end
         %% rastojanje imedju dve linije
         function rez = getRazstojanjeIzmedjuTackeILinije(pt, v1, v2)
            a = v1 - v2;
            b = pt - v2;
            d = norm(cross(a,b)) / norm(a);
         end
         %% rastojanje izmedju 2 tacke
         function [r] = rastojanjeIzmedjuDveTacke(A,B)
             if numel(B)<numel(A)
                 bckpB = B;
                 B= A;
                 B(:,1)=bckpB(1); 
                 B(:,2)=bckpB(2);
                 B(:,3)=bckpB(3);
             end
             r = A - B;
             r = r .* r;
             if numel(r(1,:)>2)
                 r=r(:,1)+r(:,2);
             else
                r=r(:,1)+r(:,2)+r(:,3);
             end
             r=sqrt(r);
         end
         %% arsDot ako treba da se skalarno pomnozi niz vektora
         %skalarni proizvod arsDot([Ax,Ay,Az],[Bx,By,Bz])=Ax*Bx+Ay*By+Az*Bz
         function [rez]=arsDot(A,B)
             if numel(A(1,:))==3%3D
                 if  numel(A(:,1))==1 && numel(B(:,1))>1
                     pom = B; pom(:,1) = A(1); pom(:,2) = A(2); pom(:,3) = A(3); A = pom;
                 elseif numel(B(:,1))==1 && numel(A(:,1))>1
                     pom = A; pom(:,1) = B(1); pom(:,2) = B(2); pom(:,3) = B(3); B = pom;
                 end
                rez = A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3);
             else%2D
                rez = A(:,1).*B(:,1) + A(:,2).*B(:,2);
             end
         end         
         function [rez]=arsNorm(A)%podrazumeva se da su vektori zapisani [ x,y,z; x,y,z...]
             if numel(A(1,:))==1
                 A=A';
             end
             if numel(A(1,:))==3%3D
                 rez = A(:,1).*A(:,1) + A(:,2).*A(:,2) + A(:,3).*A(:,3);
             else%2D
                 rez = A(:,1).*A(:,1) + A(:,2).*A(:,2);
             end
             rez = sqrt(rez);
         end 
         function [rez]=arsUnit(A)
              if numel(A(1,:))==1
                 A=A';
             end
             intenziteti = AngioIvusMath.arsNorm(A);
             rez = A ./[intenziteti intenziteti intenziteti];
         end 
         function [rez]=arsVektToJedinicniVektor(A)
             rez = A(:,1).*A(:,1) + A(:,2).*A(:,2) + A(:,3).*A(:,3);
             rez = A ./ sqrt(rez);
         end 
         %%
         function [rez]=ars2D_to_3D(tacka2D)
             rez      = [tacka2D tacka2D(:,1)];
             rez(:,3) = 0;
         end
         %%
         function [rez]=arsCross(A,B)
             if numel(A(1,:))==3%3D
                 if  numel(A(:,1))==1 && numel(B(:,1))>1
                     pom = B; pom(:,1) = A(1); pom(:,2) = A(2); pom(:,3) = A(3); A = pom;
                 elseif numel(B(:,1))==1 && numel(A(:,1))>1
                     pom = A; pom(:,1) = B(1); pom(:,2) = B(2); pom(:,3) = B(3); B = pom;
                 end
             end
             x =   A(:,2).*B(:,3)-A(:,3).*B(:,2) ;
             y = -(A(:,1).*B(:,3)-A(:,3).*B(:,1));
             z =   A(:,1).*B(:,2)-A(:,2).*B(:,1) ;             
             rez = [x y z];
         end 
         %A[brTacaka,3] B[Bx,By,Bz]
         function [rez]=arsPlus(A,B)
             if numel(A(1,:))==3%3D
                 if  numel(A(:,1))==1 && numel(B(:,1))>1
                     pom = B; pom(:,1) = A(1); pom(:,2) = A(2); pom(:,3) = A(3); A = pom;
                 elseif numel(B(:,1))==1 && numel(A(:,1))>1
                     pom = A; pom(:,1) = B(1); pom(:,2) = B(2); pom(:,3) = B(3); B = pom;
                 end
             elseif numel(A(1,:))==2%2D
                 if  numel(A(:,1))==1 && numel(B(:,1))>1
                     pom = B; pom(:,1) = A(1); pom(:,2) = A(2); A = pom;
                 elseif numel(B(:,1))==1 && numel(A(:,1))>1
                     pom = A; pom(:,1) = B(1); pom(:,2) = B(2); B = pom;
                 end
             end
             rez = A+B;
         end
         
         function [rez]=arsMinus(A,B)
             if numel(A(1,:))==3%3D
                 if  numel(A(:,1))==1 && numel(B(:,1))>1
                     pom = B; pom(:,1) = A(1); pom(:,2) = A(2); pom(:,3) = A(3); A = pom;
                 elseif numel(B(:,1))==1 && numel(A(:,1))>1
                     pom = A; pom(:,1) = B(1); pom(:,2) = B(2); pom(:,3) = B(3); B = pom;
                 end
             end
             rez = A-B;
         end
         
         function [rez]=arsPuta(A,B)
             if numel(B)==3 %drugi argument je [x,y,z]
                 pom = A;
                 pom(:,1)=B(1); pom(:,2)=B(2); pom(:,3)=B(3);
                 rez = A .* pom;
             elseif numel(B)==1
                 rez = A * B;
             elseif (numel(A)==3) & numel(B(1,:))==1 %prvi je jedan vektor N x [x y z] drugi je skalar N x a
                 brRedovaA = numel(B(:,1));
                 A = AngioIvusMath.vektorizuj(A,brRedovaA);
                 B = [B(:) B(:) B(:)];
                 rez = A .* B;
             else
                 rez = A .* B;
             end
         end
         %% mnozi dve matrice, A i B mogu biti kombinacija vektora i skalara
         % A i B su CELLS
         %ref, http://en.wikipedia.org/wiki/Matrix_multiplication Matrix product (two matrices)
         function [rez] = arsMnozenjeMatrica(A,B)
             [n  m] = size(A);
             [m2 p] = size(B);
             if m==m2
             for i=1:n
                 for j=1:p
                     for k=1:m
                         if k==1
                            rez{i,j} = A{i,k}*B{k,j};
                         else
                            rez{i,j} = rez{i,j} + A{i,k}*B{k,j}; 
                         end
                     end
                 end
             end
             else
                 rez= 'error';
             end
         end
         %% 
         function [rez] = arsMat2Cell(A)
             [I J] = size(A);
             for i=1:I
                for j=1:J
                    rez{i,j}=A(i,j);
                end
             end
         end
         %% 
         function [rez] = vektorizuj(Tacka, brRedova)
             if numel(brRedova)>1
                brRedova = numel(brRedova(:,1)); %kao brRedova moze i da se prosledi niz tacaka na kji treba da se vektorizuje
             end
             rez = zeros(brRedova,3);
             rez(:,1) = Tacka(1); rez(:,2) = Tacka(2); rez(:,3) = Tacka(3);
         end
         %% tacka preseka Linije PoP1 i ravni defisnisane sa Vo i noralom.
         % Ako je daLiSeTrazePresekSamoIzmedjUPoP1==1 ondah vraca rezultat
         % samo ako je tacka preseka izmedju PoP1. Radi i za proizvoljan
         % br tacaka/Linija. napravljeno za seckanje stl/quadrliateriala
            %INPUTS
                %linija definisana sa 2 tacke Po i P1 
                %ravan definisana sa tackom Vo i normalom n
                %opciono daLiSeTrazePresekSamoIzmedjUPoP1==1 vrati tacku preseka samo ako se nalazi izmedju PoP1 segmenta
            %OUTPUTS
                %rez            - tacka/cke 
                %idTackipreseka - id rez tacaka iz linija PoP1 
         function [rez idLinijaPreseka idLinijaKojeNeSeku]= getTackuPresekaLinijeIRavniParalelizovano(Po, P1, Vo, n, daLiSeTrazePresekSamoIzmedjUPoP1)
%              for i=1:numel(Po(:,1))
%                  plotLine(Po(i,:),P1(i,:));
%              end
             nbckp  = n ;
             Vobckp = Vo;
             idLinijaPreseka=0;
         %reff: Line-Plane Intersection
         %%http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
         %%
         if(numel(Vo(:,1))<2)%ako nije prosledjen vec paralelizovan podatak
           pomVo = Po; pomVo(:,1)=Vo(1); pomVo(:,2)=Vo(2); pomVo(:,3)=Vo(3); Vo = pomVo; clear('pomVo');     
           pomNo = Po; pomNo(:,1)=n(1);  pomNo(:,2)=n(2);  pomNo(:,3)=n(3);  n  = pomNo; clear('pomNo'); 
         end
           u = P1 - Po;
           w = Po - Vo;
           %ako su svi rezultati=0, linija je paralelna, nema tacke preseka
           pom = AngioIvusMath.arsDot(n,u);
           brTacakaPreseka = find(pom~=0);
           %%%%%%%%%%%
%            for i=1:numel(brTacakaPreseka)
%                j=brTacakaPreseka(i);
%              plotLine(Po(j,:),P1(j,:));
%            end
           %%%%%%%%%%
           if(isempty(brTacakaPreseka))
               rez = 'error';
               %linija je paralelna sa ravni, nema tacke preseka
           else%nadji tacku preseka
               Po_bckp   = Po;
               P1_bckp   = P1;
               u        = u(brTacakaPreseka,:);
               w        = w(brTacakaPreseka,:);
               Po       = Po(brTacakaPreseka,:);
               P1       = P1(brTacakaPreseka,:);
               n        = n(brTacakaPreseka,:);
               pom      = w; 
               if numel(nbckp)==3%poseban slucaj
                   pom(:,1) = nbckp(1); 
                   pom(:,2) = nbckp(2);
                   pom(:,3) = nbckp(3);
                   n        = pom;
               end
               s        = -(AngioIvusMath.arsDot(n,w)) ./ (AngioIvusMath.arsDot(n,u));%skalar
               rez      = Po + [s s s] .* u            ;
               if nargin>4 %vrati samo tacke koje su izmedju PoP1
                   idLinijaPreseka = find(s<=1 & s>=0);
                   rez             = Po(idLinijaPreseka,:) + [s(idLinijaPreseka,:) s(idLinijaPreseka,:) s(idLinijaPreseka,:)] .* u(idLinijaPreseka,:)            ;
                   idLinijaPreseka = brTacakaPreseka(idLinijaPreseka);
                   idLinijaKojeNeSeku = 1:numel(P1(:,1));
                   idLinijaKojeNeSeku = setdiff(brTacakaPreseka,idLinijaKojeNeSeku);% ???? zasta ovo treba?
%                    for i =1:numel(idLinijaPreseka)
%                        id=idLinijaPreseka(i);
%                        plotLine([Po_bckp(id,:);P1_bckp(id,:)]);
%                    end
%                    plotLine(rez);
               end
           end%if
%            p = figure;
%            plotLine(rez); hold on ;
%            plotLine([Po;P1]);
%            close(p);
         end
         %% paralelizovano ima vise tacki/linija a jedna ravan!
    function [rez idLinijaPreseka]= getTackuPresekaLinijeIRavniParalelizovano2(Po, P1, Vo, n, daLiSeTrazePresekSamoIzmedjUPoP1)
%              for i=1:numel(Po(:,1))
%                  plotLine(Po(i,:),P1(i,:));
%              end
             nbckp  = n ;
             Vobckp = Vo;
             idLinijaPreseka=0;
         %reff:http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
           
           pomPo = Vo; pomPo(:,1)=Po(1); pomPo(:,2)=Po(2); pomPo(:,3)=Po(3); Po = pomPo;  clear('pomPo');     
           pomP1 = Po; pomP1(:,1)=P1(1); pomP1(:,2)=P1(2); pomP1(:,3)=P1(3); P1  = pomP1; clear('pomP1'); 
           
           u = P1 - Po;
           w = Po - Vo;
           %ako su svi rezultati=0, linija je paralelna, nema tacke preseka
           pom = AngioIvusMath.arsDot(n,u);
           brTacakaPreseka = find(pom~=0);
           %%%%%%%%%%%
%            for i=1:numel(brTacakaPreseka)
%                j=brTacakaPreseka(i);
%              plotLine(Po(j,:),P1(j,:));
%            end
           %%%%%%%%%%
           if(isempty(brTacakaPreseka))
               rez = 'error';
               %linija je paralelna sa ravni, nema tacke preseka
           else%nadji tacku preseka
               Po_bckp   = Po;
               P1_bckp   = P1;
               u        = u(brTacakaPreseka,:);
               w        = w(brTacakaPreseka,:);
               Po       = Po(brTacakaPreseka,:);
               P1       = P1(brTacakaPreseka,:);
               pom      = w; 
               if numel(nbckp)==3%poseban slucaj
                   pom(:,1) = nbckp(1); 
                   pom(:,2) = nbckp(2);
                   pom(:,3) = nbckp(3);
                   n        = pom;
               end
               s        = -(AngioIvusMath.arsDot(n,w)) ./ (AngioIvusMath.arsDot(n,u));%skalar
               rez      = Po + [s s s] .* u            ;
               if nargin>4 
                   idLinijaPreseka = find(s<=1 & s>=0);
                   rez             = Po(idLinijaPreseka,:) + [s(idLinijaPreseka,:) s(idLinijaPreseka,:) s(idLinijaPreseka,:)] .* u(idLinijaPreseka,:)            ;
                   idLinijaPreseka = brTacakaPreseka(idLinijaPreseka);
%                    for i =1:numel(idLinijaPreseka)
%                        id=idLinijaPreseka(i);
%                        plotLine([Po_bckp(id,:);P1_bckp(id,:)]);
%                    end
%                    plotLine(rez);
               end
           end%if
           p = figure;
           plotLine(rez); hold on ;
           plotLine([Po;P1]);
           close(p);
         end
         %% upravi tacke dobijene kao presek quadrilateriala(omotacha arsPAKC-F) i ravni
            %INPUTS
                %Cvorovi                - x,y,z
                %LinijeCvorova          - node1 node2 idElementa
                %IdQuadrialetrialLinija - idQuadrialteriala iz koga je linija napravljena
                %CvoroviPreseka         - cvorovi preska linija i ravni prethodno pozvane funkcijq quadrilaterialSlicer
         function [rezTacke]=upraviTackeSlajsaQuadrilateriala(Cvorovi, LinijeCvorova, IdQuadrialetrialLinija,CvoroviPreseka)
%              for i =1:numel(LinijeCvorova(:,1))
%                  plotLine(Cvorovi(LinijeCvorova(i,[1 2])',:));
%              end
             %             node1  node2                           idFejsa                idElementa 
             pom                                  = [ sort(LinijeCvorova(:,[1 2]), 2, 'ascend') IdQuadrialetrialLinija LinijeCvorova(:,3) ];
             [bckpFejsovi idBckpFejsoviOrgLinije] = sortrows(pom);clear('pom');
             [a2 b2 c2]                           = unique(bckpFejsovi(:,3), 'first');
             [a3 b3 c3]                           = unique(bckpFejsovi(:,3), 'last');
             %A-pokazivac na susedni fejs iTe linije. treca kolona u
             %bckpFejsovi su fejsovi, prve dve su sortirani cvorivi-linije.
             %Pokazivac A omogucava da se za dve merdzovane linije npr.
             %[node1 node2 fejs1] [1 2 fejs6] sa reda fejs1 skoci na fejs6
             %i tako nastavi konektovanje
             A=bckpFejsovi(:,3);A(b2)=b3;A(b3)=b2;
             clear('a2');clear('b2');clear('c2');
             clear('a3');clear('b3');clear('c3');
             %     node1 node2          veza          fejs element
             tabela=[bckpFejsovi(:,[1 2])  A    bckpFejsovi(:,[3 4])];
             %pocni od prve tacke
             prviId         = 1     ;
             trenutniId     = prviId;
             rezTackeId     = prviId;
             rezTacke       = CvoroviPreseka(idBckpFejsoviOrgLinije(trenutniId),:);
             rez            = 666 ;
             %imporvizacija do while naredbe
             brRedovaTabele = numel(tabela(:,1));
             %broj ponadjenih kontura
             brKontura = 0;
             doWhileExitFlag =true;
             while doWhileExitFlag
                 %ako je zatvoren krug END
                 if trenutniId==0%%% ZATVORIO JE KRUG, predji na novu konturu ili ako ih nema vise izadji
                      [trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag] = ....
                            AngioIvusMath.pomFjaZaAngioIvusMath_upraviTackeSlajsaQuadrilateriala(....
                                trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag);%%kraj reda
                 end
                 trenutnaLinija = tabela(trenutniId,[1 2])                                       ;
%                  plotLine(Cvorovi(trenutnaLinija(:),:))                                          ;
                 rezTacke       = [rezTacke;CvoroviPreseka(idBckpFejsoviOrgLinije(trenutniId),:)];
%                  plotLine(rezTacke(end-1:end,:))                                                 ;
                 %ako su tacke identicne predji na susedni fejs
                 if  doWhileExitFlag && trenutniId<brRedovaTabele && mean(tabela(trenutniId,[1 2])==tabela(trenutniId+1,[1 2]))==1 
                     %predji na susedni FEJS 
                     pom               = trenutniId            ;
                     trenutniId        = tabela(trenutniId+1,3); 
                     tabela(pom,  3)   = 0                     ;
                     tabela(pom+1,  3) = 0                     ;
                     %flaguj trenutni da se ne vraca na njega                      
                 elseif trenutniId>1 &&  mean(tabela(trenutniId,[1 2])==tabela(trenutniId-1,[1 2]))==1   
                     pom               = trenutniId            ;
                     trenutniId        = tabela(trenutniId-1,3);  
                     tabela(pom-1,3)   = 0                     ;
                     tabela(pom,  3)   = 0                     ;                     
                 end
                 rezTackeId = [rezTackeId;trenutniId];
                 %ako je zatvoren krug sacuvaj konture i proveri da li ih ima josh
                 if trenutniId==0%%% ZATVORIO JE KRUG, predji na novu konturu ili ako ih nema vise izadji
                      [trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag] = ....
                            AngioIvusMath.pomFjaZaAngioIvusMath_upraviTackeSlajsaQuadrilateriala(....
                                trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag);%%kraj reda
                 end
             end
         end
       %% pomocna funkcija za AngioIvusMath.upraviTackeSlajsaQuadrilateriala
       function [trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag] = pomFjaZaAngioIvusMath_upraviTackeSlajsaQuadrilateriala(trenutniId,prviId,rezTacke,rezTackeId,tabela,brKontura,rez,doWhileExitFlag)
           %ako je trenutni Id jednak pocetnom znaci da je zatvorio krug
%            if trenutniId==prviId%%% ZATVORIO KRUG 
%                 plotLine(rezTacke)          ;
%                 tabela(rezTackeId,3)=0                   ;
                next               = find(tabela(:,3)>0) ;
                if isempty(next)%%% KRAJ-nema vise kontura-----------------
                   if brKontura ==0%ako je nasao samo 1
                     brKontura     = brKontura+1      ;                              
                     slajs.kontura = rezTacke         ;
                     rez           = slajs            ;
                   else%ako ih ima vishe                              
                     slajs.kontura = rezTacke(2:end,:);
                     rez           = [rez;slajs]      ;
                   end
                   %%%%%%%%%%%%%KRAJ SECKANjA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   doWhileExitFlag = false;%%%NEMA VISE KONTURA IZADjI IZ ROOT FJE
                   %%% profiltriraj konture(izbrisi iste tacke)
                   brPronadjenihKontura=numel(rez);
                   %tacke su u stekovima vec sortirane, ali ima duplih
                   %susednih cvorova. nadji unique tacke, i onda iz tih
                   %Idova izdvoj jedinstvene-sortirane tacke
                   for i = 1:brPronadjenihKontura                     
                       [a b c]              = unique(rez(i).kontura, 'rows','first');
                       idJedisntvenihTacaka = sort(b)                               ;
                       rez(i).kontura       = rez(i).kontura(idJedisntvenihTacaka,:);
%                        plotLine(rez(i).kontura );
                   end
                   %Merdzuj konture. u posebnim slucajevima moze da nadje 2 iste konture 
                   clear('rezTacke');%reset
                   rezTacke(1).kontura = rez(1).kontura;
                   brJedinstvenihKontura=1;
                   for i=2:brPronadjenihKontura
                       %vidi da li vec postoji takva kontura u rez steku
                       daLiVecPostojiItaKonturaUSteku=false;
                       for j=1:numel(rezTacke)
                           [a b c] = unique([rezTacke(j).kontura;rez(i).kontura], 'rows','first');
                           if numel(b) == (numel(c)/2)
                               daLiVecPostojiItaKonturaUSteku = true;
                           end
                       end
                       %ako ne postoji ista kontura dodaj je
                       if ~daLiVecPostojiItaKonturaUSteku
                           brJedinstvenihKontura                   = brJedinstvenihKontura+1;
                           rezTacke(brJedinstvenihKontura).kontura = rez(i).kontura         ;
                       end
                   end
%                    for j=1:numel(rezTacke)
%                       plotLine(rezTacke(j).kontura );
%                    end                           
                   return;%------------------------------------------------
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else%%%%% IMA JOSH KONTURA, resetuj indikatore i nastavi
%                    figure;
                   trenutniId = next(1)   ;%uzmi prvi i nastavi
                   prviId     = trenutniId;
                   %sacuvaj konturu
                   if brKontura ==0
                     brKontura     = brKontura+1      ;                              
                     slajs.kontura = rezTacke         ;
                     rez           = slajs            ;
                   else                              
                     slajs.kontura = rezTacke(2:end,:);
                     rez           = [rez;slajs]      ;
                   end
                   rezTacke=[0 0 0];
                end
%            end
       end
       %% QUADRILATERIAL SLICER (SECKA OMOTAC DATA)
        %vraca tacke preseka ravni(definisane sa normala,binormala,tangenta iz
        %CenterLine) i Quadrilateriala koji opisuju omotach modela
        %INPUTS 
            %Cvorovi         - jedinstvena lista cvorova[brCvorova x3 ] {x y z}
            %Quadrilateriali - fejsovi na obodu [brFejsova x 5] { node1 node2 node3 node4 IdElemantaKomeFejsPripada}
         	%CentarLine      - linja koja definise B-Spline
        %OUTPUTS
            %rezTacke        = rez tacke
        function  [rez] = quadrilaterialCenterlineSlicer(Cvorovi,Quadrilateriali,CentarLine)
            %%  NAPRAVI OD QUADRILATERIALA LINIJE(DVA CVRORA)
                % 1. napravi Linije
                %selektor{idNode1 idNode2 Quadrilateriali(5)=idElementa}
                selektor= [1 2 5; 2 3 5; 3 4 5; 4 1 5];
                Linije  = Quadrilateriali(:,selektor(1,:)); % Linija=[idNode1 idNode2 idElementa]
                pom = 1:numel(Quadrilateriali(:,1));
                pom =pom(:);
                idLinijeQuadrilaterial=pom;                
                for i=2:4
                    Linije =[Linije; Quadrilateriali(:,selektor(i,:))];%br linija je 4*brQuadrialteriala
                    idLinijeQuadrilaterial=[idLinijeQuadrilaterial;pom];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %provera za Quad
%                 for i =1:numel(Quadrilateriali(:,1))
%                    plotLine(Cvorovi(Quadrilateriali(i,[1 2])',:));
%                    plotLine(Cvorovi(Quadrilateriali(i,[2 3])',:));
%                    plotLine(Cvorovi(Quadrilateriali(i,[3 4])',:));
%                    plotLine(Cvorovi(Quadrilateriali(i,[4 1])',:));
%                 end
                  %provera za Linije
                  brQuat = numel(Quadrilateriali(:,1));
%                 for i =1:brQuat
%                    plotLine(Cvorovi(Linije(i,[1 2]),:));
%                    plotLine(Cvorovi(Linije(i+brQuat,[1 2]),:));
%                    plotLine(Cvorovi(Linije(i+2*brQuat,[1 2]),:));
%                    plotLine(Cvorovi(Linije(i+3*brQuat,[1 2]),:));
%                 end
                clear('pom');
                % 2. napraiv normale binormale ravni
                if isstr(CentarLine)
                    %ucitaj iz txta
                else
        %             splajn                  = arsCubicBSplineOpenNonuniform;
        %             splajn                  = splajn.set(splajn, CentarLine,3);
        %             normaleBinormaleSplajna = splajn.testQuaternionOrjentacijuPreseka(splajn, CentarLine, 3, [0.01:0.1:1]);
        %             splajn                  = normaleBinormaleSplajna.sveTackeSplajna;
                    splajn = CentarLine;
        %               splajn = load('splajn_arsVTK_215red.mat');
        %               splajn = splajn.splajn;
                end
                normale = CentarLine;
                %za svaku od linija nadji presek sa ravni na segmentu node1 node2
                brpodelaSplajna=numel(normale(:,1));
                for i = 1:brpodelaSplajna            
                    Vo = splajn(i,:);
                    n  = normale(i,:);
                    [tackeItogPreseka idLinijaIzKojihSuDobijeneTacke] = AngioIvusMath.getTackuPresekaLinijeIRavniParalelizovano(Cvorovi(Linije(:,1),:), Cvorovi(Linije(:,2),:), Vo, n,666);
                    %iscrtaj quadrialteriale koji su preseceni
% %                for j =1:numel(idLinijaIzKojihSuDobijeneTacke(:,1))
% %                    i=idLinijaIzKojihSuDobijeneTacke(j);
% %                    plotLine(Cvorovi(Linije(i,[1 2])',:));
% %                    plotLine(Cvorovi(Linije(i,[1 2])',:));
% %                    plotLine(Cvorovi(Linije(i,[1 2])',:));
% %                    plotLine(Cvorovi(Linije(i,[1 2])',:));
% %                end
% % %                close all;
% %                figure;
%                for j =1:numel(idLinijaIzKojihSuDobijeneTacke(:,1))
%                    i=idLinijaIzKojihSuDobijeneTacke(j);
%                    i=idLinijeQuadrilaterial(i);
%                    quad = Quadrilateriali(i,1:4);
%                    plotLine(Cvorovi(quad([1 2])',:));
%                    plotLine(Cvorovi(quad([2 3])',:));
%                    plotLine(Cvorovi(quad([3 4])',:));
%                    plotLine(Cvorovi(quad([4 1])',:));
%                end
%                plotLine(tackeItogPreseka);                
                    %dobijene tacke organizuj u zatvorene konture
                    %-ako ne postoje tacke preseka zadaj povrsinu nula a
                    %tacke konture stavi error
                    if isstr(tackeItogPreseka)%==error
                        rezPovrsineKontura  =  0      ;
                        kontureItogSlajsa   = 'error' ;
                    else
                        kontureItogSlajsa = AngioIvusMath.upraviTackeSlajsaQuadrilateriala(...
                            Cvorovi,...
                            Linije(idLinijaIzKojihSuDobijeneTacke,:),...
                            idLinijeQuadrilaterial(idLinijaIzKojihSuDobijeneTacke),...
                            tackeItogPreseka);
                        %povrsine po konturama(po slajsovima moze imati vise kontura
                        for brKonture=1:numel(kontureItogSlajsa)
                            %----------------------------------------                                                  3D tacke,        normala  
                            rezPovrsineKontura(brKonture) = AngioIvusMath.getPovrsinuPoligonaURavni(kontureItogSlajsa(brKonture).kontura,   n   );
                        end
                    end
                    %napakuj rezultate
                    slajs.povrsineKontura= rezPovrsineKontura;
                    slajs.konture        = kontureItogSlajsa ;
                    rez(i)               = slajs             ;
                    end
                
        end
         %% presek linije i ravni
          %INPUTS 
            %Linija - definisana tackama Po i P1
            %Ravan  - definisana sa tackom Vo i vektor normale n
          %OUTPUTS
            %rezTackePreseka - rez tacke
            %id              - idTacki koje su izmedju PoP1 
            %tParametar      - ako je vece od nule onda je tacka preseka ispred n
         function [rezTackePreseka id tParametar] = getTackuPresekaLinijeIRavni(Po, P1, Vo, n,modVratiSveTacke)
             ars =0;
             if numel(Po(:,1))>1 && numel(Po)>3%vektorizacija
                 pom = Po;
                 pom(:,1)=Vo(1); pom(:,2)=Vo(2); pom(:,3)=Vo(3); Vo = pom;
                 n = AngioIvusMath.arsUnit(n);
                 pom(:,1)=n(1) ; pom(:,2)=n(2) ; pom(:,3)=n(3) ; n  = pom;
%                  pom(:,1)=P1(1); pom(:,2)=P1(2); pom(:,3)=P1(3); P1 = pom;
                if numel(P1(:,1))==1
                    pom(:,1)=P1(1); pom(:,2)=P1(2); pom(:,3)=P1(3); P1 = pom;
                end
             elseif numel(Vo(:,1))>1 && numel(Vo)>3%vektorizacija
                  pom = Vo; pom(:,1)=Po(1); pom(:,2)=Po(2); pom(:,3)=Po(3); Po = pom;
                  pom = Vo; pom(:,1)=P1(1); pom(:,2)=P1(2); pom(:,3)=P1(3); P1 = pom;
             else
                 P1=P1(:)'; Po=Po(:)'; Vo=Vo(:)'; n=n(:)'; 
             end
         %reff: %http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
           u = P1 - Po;
           w = Po - Vo;
           if(sum(AngioIvusMath.arsDot(n,u)) == 0) 
               rez = 'error';
               tParametar = rez;
               id = rez;
               %linija je paralelna sa ravni, nema tacke preseka
           else%nadji tacku preseka
               s   = -(AngioIvusMath.arsDot(n,w))./(AngioIvusMath.arsDot(n,u)); 
               rez = Po + [s(:) s(:) s(:)] .* u;
               tParametar = s(:);
           end%if
           if ars
               rez = rez';
           end
           rezTackePreseka = rez;
           id = find( s>=0 & s<=1);
%            p = figure;
%            plotLine(rez,P1); hold on ;
%            set(p,'Color','red');
%            plotLine(Po,P1);
         end
         %% ged presek linije i sfere
         %http://www.ahinson.com/algorithms_general/Sections/Geometry/IntersectionOfParametricLineAndSphere.pdf
            %Linija (Po, d) 
                %Po - ta?ka iz koje kre?e linija (PoX, PoY, Po)
                %d  - jedini?ni vektor pravca    (dX ,  dy, dz)
            %Sfera  (C,  r)
                %C  - centar sfere (Cx, Cy, Cz)
                %d  - radijus (skalar)
         function [idSferaKojeSeceLinija, tackaPreseka1, tackaPreseka2] = getPresekLinijeSfere(Po, d, C, r)
             a = sum( [d .* d]');
             b = 2 * sum ( [ ( d .* (Po - C)) ]');
             c = sum( [(Po-C) .* (Po-C)]' ) - r.*r;
             %
             t = b.*b - 4*(a.*c);
             idGdeSece = find(t>=0);
             t1 = (-b(idGdeSece) + sqrt(t(idGdeSece))) ./ (2*a(idGdeSece));
             t2 = (-b(idGdeSece) - sqrt(t(idGdeSece))) ./ (2*a(idGdeSece));
             % rezultati
             idSferaKojeSeceLinija = idGdeSece;
             tackaPreseka1 = Po(idGdeSece,:) + d(idGdeSece,:) .* [ t1(:) t1(:) t1(:)];
             tackaPreseka2 = Po(idGdeSece,:) + d(idGdeSece,:) .* [ t2(:) t2(:) t2(:)];
         end
         %% nalazi presek trougla/trouglova{T1,T2,T3} sa ravni{Po,P1}
         %INPUTS
            %Trougao/Trouglovi{T1,T2,T3} moze i vise troulova iz STLa
            %Ravan{Po,P1} - SAMO dve tacke
         %OUTPUTS
            %tacke preseka linija trougla i ravni
         function [rez idTrougla]=getPresekTrouglaIRavni(T1,T2,T3,Po,P1)
             pom = T1;
             pom(:,1)=P1(1);pom(:,2)=P1(2);pom(:,3)=P1(3);
             P1 = pom;
             pom(:,1)=Po(1);pom(:,2)=Po(2);pom(:,3)=Po(3);
             Po = pom;
             normalaRavni       = P1-Po;%AngioIvusMath.arsUnit(P1-Po); 
             [tackaPresekaT1T2 idTrouglaT1T2]  = AngioIvusMath.getTackuPresekaLinijeIRavniParalelizovano(T1, T2, Po, normalaRavni,'uzmiSamoTackeIzmedju');%getTackuPresekaLinijeIRavniParalelizovano
             [tackaPresekaT1T3 idTrouglaT1T3]  = AngioIvusMath.getTackuPresekaLinijeIRavniParalelizovano(T1, T3, Po, normalaRavni,'uzmiSamoTackeIzmedju');%getTackuPresekaLinijeIRavni
             [tackaPresekaT3T2 idTrouglaT3T2]  = AngioIvusMath.getTackuPresekaLinijeIRavniParalelizovano(T3, T2, Po, normalaRavni,'uzmiSamoTackeIzmedju');
             rez = [ 0 0 0];
             rez = unique([tackaPresekaT1T2;tackaPresekaT1T3;tackaPresekaT3T2], 'rows');
             sveTacke     = [tackaPresekaT1T2;tackaPresekaT1T3;tackaPresekaT3T2];
             sviTrouglovi = [idTrouglaT1T2;idTrouglaT1T3;idTrouglaT3T2]; 
             [a b] = unique(sviTrouglovi);
             sviTrouglovi = sviTrouglovi(b);
             idTrougla = sviTrouglovi;
             for i=1:numel(sviTrouglovi)
                 plotLine( [T1(sviTrouglovi(i),:); T2(sviTrouglovi(i),:); T3(sviTrouglovi(i),:);T1(sviTrouglovi(i),:)]);
             end
% %              if ~isstr(tackaPresekaT1T2)
% %                  pom = abs((AngioIvusMath.arsNorm(T1-T2)-(AngioIvusMath.arsNorm(tackaPresekaT1T2-T1)+AngioIvusMath.arsNorm(tackaPresekaT1T2-T2))))<0.00001;
% %                  pom = find(pom>0);
% %                  if ~isempty(pom)
% %                     rez=[rez;tackaPresekaT1T2(pom,:)];
% %                  end
% %              end
% %              if ~isstr(tackaPresekaT1T3)
% %                  pom = abs((AngioIvusMath.arsNorm(T1-T3)-(AngioIvusMath.arsNorm(tackaPresekaT1T3-T1)+AngioIvusMath.arsNorm(tackaPresekaT1T3-T3))))<0.00001;
% %                  pom = find(pom>0);
% %                  if ~isempty(pom)
% %                     rez=[rez;tackaPresekaT1T3(pom,:)];
% %                  end
% %              end
% %              if ~isstr(tackaPresekaT3T2)
% %                  pom = abs((AngioIvusMath.arsNorm(T3-T2)-(AngioIvusMath.arsNorm(tackaPresekaT3T2-T3)+AngioIvusMath.arsNorm(tackaPresekaT3T2-T2))))<0.00001;
% %                  pom = find(pom>0);
% %                  if ~isempty(pom)
% %                     rez=[rez;tackaPresekaT3T2(pom,:)];
% %                  end
% %              end
             if numel(rez)==3
                 rez='error';
             else
                 rez=rez(2:end,:,:);
             end             
         end
         %% presek dve linije P i Q na segmentima PoP1 i QoQ1 
         %INPUTS
           %Linija P: tacke Po i P1
           %Linija Q: tacke Qo i Q1
         %OUTPUTS
           %rez:tacka preseka dve linije koja se nalazi na segmentima PoP1 i
           %QoQ1 linija P i Q. Ukoliko se ne seku na datim segmentmia vraca
           %informaciju o gresci
         %reff:http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm
         function [rez] = getTackuPresekaLinijaNaSegmentima(Po,P1,Qo,Q1)
            [Psc, Qtc, M] = AngioIvusMath.najblizaTackaIzmedjuLinija(Po, P1, Qo, Q1);
            rez = 'error';
           if((norm(Psc-Qtc) + norm(Qtc-M))<0.0001)%linije se seku
               %da li se tacka M nalazi na segmentu PoP1 i QoQ1???
               %racunaju se duzine duzi pa se uporedjuje da li je zbir
               %duzina MPo i MP1 jednak PoP1(da li su tacke izmedju)
               MPo  = norm(Po-M) ;
               MP1  = norm(P1-M) ;
               PoP1 = norm(P1-Po);
               %ako se tacka preseka nalazi na segmentu PoP1 racunaj za QoQ1
               if(abs((MPo+MP1)-PoP1)<0.001)
                   MQo  = norm(Qo-M) ;%duzina MPo duzi
                   MQ1  = norm(Q1-M) ;
                   QoQ1 = norm(Q1-Qo);                  
                  if(abs((MQo+MQ1)-QoQ1)<0.001)
                      rez = M;
                  end
               else% linije se seku ali tacka preseka nije na oba segmentima
                   rez = 'error';
               end 
           else%linije se ne seku vrati error
               rez = 'error';
           end    
          end%getTackuPresekaLinijaNaSegmentima
         %% tacka preseka linije L(PoP1) i cetvorougla(ABCD)
         %podrazumeva se da su tacke ABCD i linija PoP1 u ravni
         %INPUTS:
            %Linija: P1,P2      = P1....P2
            %Cetvorougao:ABCD   =  A....B
            %                      :    :
            %                      C....D
         %OUTPUTS:
            %rez: T1, T2 tacke preseka ili 'error' ako ih nema
         function [rez] = getTackePresekaLinijeP1P2iCetvorouglaABCD(Po,P1,A,B,C,D)
             %da li se linija(P1,P2) nalazi u ravni (ABCD)
             n = cross(A-B, C-B); %vektor normale na ravan             
             if(dot(P1-Po, n) == 0)%ako je vektor P2P1 normalan na vekor normale ravni 
                 ABCD = [A;B;C;D;A];
                 TackePreseka = [666,666,666];%pomocna inicijaliacija
                 for i=1:4
                     T = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(Po,P1,ABCD(i,:,:),ABCD(i+1,:,:));
                     if(~strcmp(T,'error'))
                         TackePreseka = [ TackePreseka; T ];
                     end
                 end
                 TackePreseka = TackePreseka(2:3,:,:);
                 if (numel(TackePreseka(:,1,:))==2)%ako postoje dve presecne tacke
                    rez = TackePreseka;
                 else %ne sece ABCD u dve tacke
                    rez = 'error'; 
                 end
             else %Linija nije u istoj ravni sa ABCD                 
                 rez = 'error'; %nema tacko preseka
             end%if
         end%getTackePresekaLinijeP1P2iCetvorouglaABCD
         %% tacka preseka 3D linije sa trouglom u prostoru
         % linija i trougao NE SMEJU biti u istoj ravni
         %INPUTS
            %Linija: Po,P1
            %Trougao:A,B,C
         %Outputs: tacka preseka ili 'error' ako je nema
         function [rez id tParametar] = getTackuPresekaLinijeP1P2iTrouglaABC(Po,P1,A,B,C)
             n = AngioIvusMath.arsCross(A-B, C-A); %vektor normale na ravan 
             %izbaci trouglove koji su apralelni
             kandidati = find(AngioIvusMath.arsDot(P1-Po, n) ~= 0);
             n  = n(kandidati,:);
             A  = A(kandidati,:);
             B  = B(kandidati,:);
             C  = C(kandidati,:);
             if(AngioIvusMath.arsDot(P1-Po, n) ~= 0)%ako vektor PoP1 NIJE normalan na vekor normale trougla(ravni) 
                  [TackaPreseka  id tParametar] = AngioIvusMath.getTackuPresekaLinijeIRavni(Po,P1, A,n);
                  daLijeUTrouglu = AngioIvusMath.daLiSeTackaNalaziUnutarTrougla(A,B,C,TackaPreseka);
                  if(~daLijeUTrouglu)%ukoliko nije u trouglu 
                     TackaPreseka = 'error';
                  end
             end
             rez = TackaPreseka;
         end%getProjekcijuTackeNaPravu         
        %% da li se tacke P1 i P2 nalaze na istoj strani duzi AB
        %ref:http://www.blackpawn.com/texts/pointinpoly/default.html
        %INPUTS 
            %Linija - definisana tackama A i B
            %Tacke  - P1 i P2
        function [rez id] = daLiSeTackeNalazeNaistojStraniDuziAB(A,B,P1,P2)
            cp1 = AngioIvusMath.arsCross(AngioIvusMath.arsMinus(B,A), AngioIvusMath.arsMinus(P1,A));
            cp2 = AngioIvusMath.arsCross(AngioIvusMath.arsMinus(B,A), AngioIvusMath.arsMinus(P2,A));
            rez = zeros(numel(cp2(:,1)),1);
            id = find(AngioIvusMath.arsDot(cp1,cp2) >=-0.00001);
            rez(id)=1;
        end%daLiSeTackeNalazeNaistojStraniDuziAB
        %% da li se tacka T nalazi na liniji AB izmedju te dve tacke?
        %INPUTS
            %A,B - tacke koje definisu liniju
            %T   - jedna ili vise tacki za koje se proverava
        %OUTPUTS
            %rezBool                         - vraca [0 1 0 0 1 1....numel(T)]
            %rezIdTacakaNaSemgmentu          - id tacaka na segmentu = find(rezBool==true)
            %rezIdTacakaNaSemgmentuSortirano - sortirano od A ka B
        function [rezBool rezIdTacakaNaSemgmentu rezIdTacakaNaSemgmentuSortirano] = daLiSeTackaNalaziIzmedjuDuziAB(A,B,T)
            %ako se poziva vektorizovano
            pomA=T; pomA(:,1)=A(1); pomA(:,2)=A(2);  pomA(:,3)=A(3); A=pomA;
            pomB=T; pomB(:,1)=B(1); pomB(:,2)=B(2);  pomB(:,3)=B(3); B=pomB;
            %ako se tacka nalazi na liniji onda je vektorski proizvod AT AB jednak nuli
            AB = B-A;
            AT = T-A;
            daLiJeNaLiniji   = AngioIvusMath.arsCross(AB,AT)        ;%daLiJeNaLiniji=0, znaci da se nalazi na liniji
            daLiJeNaLiniji   = AngioIvusMath.arsNorm(daLiJeNaLiniji);
            idTacakaNaLiniji = find(abs(daLiJeNaLiniji)<0.000001)   ;
            if ~isempty(daLiJeNaLiniji)
                %nadji tacke koje se nalaze na segmentu 
                jedinicniAB = AngioIvusMath.arsUnit(AB);
                t           = (T - A) ./ AB            ;
                t           = t(:,1)                   ;
                tNaSegmentu = find(t>=0 & t<=1)        ;
            end
            rezIdTacakaNaSemgmentu          = intersect(idTacakaNaLiniji,tNaSegmentu);
            rezBool                         = zeros(numel(t),1)                      ;
            rezBool(rezIdTacakaNaSemgmentu) = 1                                      ;
            [ a idSortirano]                = sort(t(rezIdTacakaNaSemgmentu))        ;
            rezIdTacakaNaSemgmentuSortirano = rezIdTacakaNaSemgmentu(idSortirano)     ;
        end
        %% da li se tacka T nalazi unutar trougla ABC
        %ref:http://www.blackpawn.com/texts/pointinpoly/default.html
        %INPUTS 
            %Trougao - definisana tackama A, B i C
            %Tacka   - T
        function [rez] = daLiSeTackaNalaziUnutarTrougla(A,B,C, T)
            rez= 0;
            if(AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(A,B,C,T) &&...
                    AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(B,C,A,T) &&...
                    AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(C,A,B,T))
                rez=1;return;                
                else
                    rez=0;return;
            end
            %-pokusaj paralelizacije
            for i = 1:numel(T(:,1))
                iT = T; iT(:,1)=T(i,1); iT(:,2)=T(i,2); iT(:,3)=T(i,3);
                
                
%                   daLije =  AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(A,B,C,T(i,:)) .*...
%                             AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(B,C,A,T(i,:)) .*...
%                             AngioIvusMath.daLiSeTackeNalazeNaistojStraniDuziAB(C,A,B,T(i,:));
%                 rez(i) = sum(daLije);%moze da bude samo u jednom trouglu
            end
        end%
        %%
        function [rez] = daLiJeTackaUnutarRavniABCD(A,B,C,D, tacka)
        %reff:http://www.blackpawn.com/texts/pointinpoly/default.html
        %proverava da li je tacka unutar trouglova
            %da li je tacka u trouglu ABC            
            %  A....B
            %  :    :
            %  C....D
            if( AngioIvusMath.daLiSeTackaNalaziUnutarTrougla(A,B,C,tacka) ||...
                AngioIvusMath.daLiSeTackaNalaziUnutarTrougla(B,C,D,tacka))
                rez = 1;
            else
                rez = 0;
            end%
        end%daLiJeTackaUnutarRavniABCD
        %% isto kao daLiSuTackeUnutarTetre samo sto je prosledjen niz Tetri
        function [rez] = daLiSuTackeUnutarTetri(Cvorovi, Tetra, Tacke)  
            rez = daLiSuTackeUnutarTetre(Cvorovi, Tetra, Tacke);%prekucao sam funkciju u C++ tako da moze da radi za jednu/vise
        end
        %% CircleFitByPratt(XY)
        %INPUTS
            %XY-2D koordinate tacki
            %brPodelaZaRezTacke-krug treba da se fituje i da se vrati brPodelaZaRezTacki na krugu
        %OUTPUTS
            %Par[a b R] - a,b centar nafitovanog kruga, R precnik
            %rezTacke   - tacke kruga
        function [Par rezTacke] = fitujKrugNaTacke2D(XY, brPodelaZaRezTacke)
        debugMode=0;
        %http://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit--pratt-method-
        %--------------------------------------------------------------------------
        %  
        %     Circle fit by Pratt
        %      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
        %      Computer Graphics, Vol. 21, pages 145-152 (1987)
        %
        %     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
        %
        %     Output: Par = [a b R] is the fitting circle:
        %                           center (a,b) and radius R
        %
        %     Note: this fit does not use built-in matrix functions (except "mean"),
        %           so it can be easily programmed in any programming language
        %
        %--------------------------------------------------------------------------

        n = size(XY,1);      % number of data points
        centroid = mean(XY);   % the centroid of the data set
        %     computing moments (note: all moments will be normed, i.e. divided by n)
        Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;
        for i=1:n
            Xi = XY(i,1) - centroid(1);  %  centering data
            Yi = XY(i,2) - centroid(2);  %  centering data
            Zi = Xi*Xi + Yi*Yi;
            Mxy = Mxy + Xi*Yi;
            Mxx = Mxx + Xi*Xi;
            Myy = Myy + Yi*Yi;
            Mxz = Mxz + Xi*Zi;
            Myz = Myz + Yi*Zi;
            Mzz = Mzz + Zi*Zi;
        end
        Mxx = Mxx/n;
        Myy = Myy/n;
        Mxy = Mxy/n;
        Mxz = Mxz/n;
        Myz = Myz/n;
        Mzz = Mzz/n;
        %    computing the coefficients of the characteristic polynomial
        Mz = Mxx + Myy;
        Cov_xy = Mxx*Myy - Mxy*Mxy;
        Mxz2 = Mxz*Mxz;
        Myz2 = Myz*Myz;

        A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
        A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
        A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
        A22 = A2 + A2;

        epsilon=1e-12; 
        ynew=1e+20;
        IterMax=20;
        xnew = 0;
        %    Newton's method starting at x=0
        for iter=1:IterMax
            yold = ynew;
            ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
            if (abs(ynew)>abs(yold))
                disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
                xnew = 0;
                break;
            end
            Dy = A1 + xnew*(A22 + 16*xnew*xnew);
            xold = xnew;
            xnew = xold - ynew/Dy;
            if (abs((xnew-xold)/xnew) < epsilon), break, end
            if (iter >= IterMax)
                disp('Newton-Pratt will not converge');
                xnew = 0;
            end
            if (xnew<0.)
                fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
                xnew = 0;
            end
        end
        %    computing the circle parameters
        DET = xnew*xnew - xnew*Mz + Cov_xy;
        Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;
        Par    = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];
        if ~exist('brPodelaZaRezTacke')
            brPodelaZaRezTacke=16;
        end
        pomXZ=XY; XY(:,3)=0;
        srednjiPrecnikUlaza=mean(AngioIvusMath.arsNorm(pomXZ));
        if Par(3)/srednjiPrecnikUlaza > 1.217 
            Par(3)=srednjiPrecnikUlaza;
        end
        if norm(Par(1:2))/srednjiPrecnikUlaza>0.16
            Par(1)=0;
            Par(2)=0;
        end
        alfa=0:2*pi/brPodelaZaRezTacke:2*pi; alfa=alfa(1:end-1);
        x=Par(3)*cos(alfa)+Par(1);
        y=Par(3)*sin(alfa)+Par(2);
        rezTacke=[x(:) y(:)]; rezTacke(:,3)=0;
        if debugMode
            figure; axis equal; hold on;            
            plot(XY(:,1), XY(:,2), 'rS');
            plotLine(rezTacke);
        end
        end%CircleFitByPratt
        %% FITOVANjE 3D kruga (za potrebe generisanja patchve)
        %% fitovanje elipse
        function [ellipse_t rezTacke] = fitujElipsuNaTacke2D(XY, brPodelaZaRezTacke )
        % fit_ellipse - finds the best fit to an ellipse for the given set of points.
        %
        % Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
        %
        % Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
        %           axis_handle - optional. a handle to an axis, at which the estimated ellipse 
        %                         will be drawn along with it's axes
        %
        % Output:   ellipse_t - structure that defines the best fit to an ellipse
        %                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
        %                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
        %                       phi         - orientation in radians of the ellipse (tilt)
        %                       X0          - center at the X axis of the non-tilt ellipse
        %                       Y0          - center at the Y axis of the non-tilt ellipse
        %                       X0_in       - center at the X axis of the tilted ellipse
        %                       Y0_in       - center at the Y axis of the tilted ellipse
        %                       long_axis   - size of the long axis of the ellipse
        %                       short_axis  - size of the short axis of the ellipse
        %                       status      - status of detection of an ellipse
        %
        % Note:     if an ellipse was not detected (but a parabola or hyperbola), then
        %           an empty structure is returned

        % =====================================================================================
        %                  Ellipse Fit using Least Squares criterion
        % =====================================================================================
        % We will try to fit the best ellipse to the given measurements. the mathematical
        % representation of use will be the CONIC Equation of the Ellipse which is:
        % 
        %    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
        %   
        % The fit-estimation method of use is the Least Squares method (without any weights)
        % The estimator is extracted from the following equations:
        %
        %    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
        %
        %    where:
        %       A   - is the vector of parameters to be estimated (a,b,c,d,e)
        %       x,y - is a single measurement
        %
        % We will define the cost function to be:
        %
        %   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
        %            = (X*A+f_c)'*(X*A+f_c) 
        %            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
        %
        %   where:
        %       g_c(x_c,y_c;A) - vector function of ALL the measurements
        %                        each element of g_c() is g(x,y;A)
        %       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
        %       f_c            - is actually defined as ones(length(f),1)*f
        %
        % Derivation of the Cost function with respect to the vector of parameters "A" yields:
        %
        %   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
        %
        % Which yields the estimator:
        %
        %       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
        %       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %
        % (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
        %  
        % NOW, all that is left to do is to extract the parameters from the Conic Equation.
        % We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
        %
        %    Recall the conic representation of an ellipse:
        % 
        %       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
        % 
        % We will check if the ellipse has a tilt (=orientation). The orientation is present
        % if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
        % tilt of the ellipse.
        %
        % If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
        % we will remove the tilt of the ellipse so as to remain with a conic representation of an 
        % ellipse without a tilt, for which the math is more simple:
        %
        % Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
        %
        % We will remove the orientation using the following substitution:
        %   
        %   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
        %   
        %   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
        %
        %   where:      c = cos(phi)    ,   s = sin(phi)
        %
        %   and simplify...
        %
        %       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
        %           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
        %
        %   The orientation is easily found by the condition of (B_new=0) which results in:
        % 
        %   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
        %   
        %   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
        %   all the other constants A`,C`,D`,E` can be found.
        %
        %   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
        %   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c 
        %   C` = A*s^2 + B*c*s + C*c^2
        %
        % Next, we want the representation of the non-tilted ellipse to be as:
        %
        %       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
        %
        %       where:  (X0,Y0) is the center of the ellipse
        %               a,b     are the ellipse "radiuses" (or sub-axis)
        %
        % Using a square completion method we will define:
        %       
        %       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
        %
        %       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
        %                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
        %
        %       which yields the transformations:
        %       
        %           X0  =   -D`/(2*A`)
        %           Y0  =   -E`/(2*C`)
        %           a   =   sqrt( abs( F``/A` ) )
        %           b   =   sqrt( abs( F``/C` ) )
        %
        % And finally we can define the remaining parameters:
        %
        %   long_axis   = 2 * max( a,b )
        %   short_axis  = 2 * min( a,b )
        %   Orientation = phi
        %
        %

        % initialize
        orientation_tolerance = 1e-3;

        % empty warning stack
        warning( 'AA' );

        % prepare vectors, must be column vectors
        x = XY(:,1);
        y = XY(:,2);

        % remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
        mean_x = mean(x);
        mean_y = mean(y);
        x = x-mean_x;
        y = y-mean_y;

        % the estimation for the conic equation of the ellipse
        X = [x.^2, x.*y, y.^2, x, y ];
        a = sum(X)/(X'*X);

        % check for warnings
        if ~isempty( lastwarn )
        %     disp( 'stopped because of a warning regarding matrix inversion' );
        %     ellipse_t = [];
        %     return
        end

        % extract parameters from the conic equation
        [a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

        % remove the orientation from the ellipse
        if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )

            orientation_rad = 1/2 * atan( b/(c-a) );
            cos_phi = cos( orientation_rad );
            sin_phi = sin( orientation_rad );
            [a,b,c,d,e] = deal(...
                a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
                0,...
                a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
                d*cos_phi - e*sin_phi,...
                d*sin_phi + e*cos_phi );
            [mean_x,mean_y] = deal( ...
                cos_phi*mean_x - sin_phi*mean_y,...
                sin_phi*mean_x + cos_phi*mean_y );
        else
            orientation_rad = 0;
            cos_phi = cos( orientation_rad );
            sin_phi = sin( orientation_rad );
        end

        % check if conic equation represents an ellipse
        test = a*c;
        switch (1)
        case (test>0),  status = 'Sve ok';
        case (test==0), status = 'Parabola found';  warning( 'fit_ellipse: Did not locate an ellipse' );
        case (test<0),  status = 'Hyperbola found'; warning( 'fit_ellipse: Did not locate an ellipse' );
        end

        % if we found an ellipse return it's data
        if (test>0)

            % make sure coefficients are positive as required
            if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end

            % final ellipse parameters
            X0          = mean_x - d/2/a;
            Y0          = mean_y - e/2/c;
            F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
            [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );    
            long_axis   = 2*max(a,b);
            short_axis  = 2*min(a,b);

            % rotate the axes backwards to find the center point of the original TILTED ellipse
            R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
            P_in        = R * [X0;Y0];
            X0_in       = P_in(1);
            Y0_in       = P_in(2);

            % pack ellipse into a structure
            ellipse_t = struct( ...
                'a',a,...
                'b',b,...
                'phi',orientation_rad,...
                'X0',X0,...
                'Y0',Y0,...
                'X0_in',X0_in,...
                'Y0_in',Y0_in,...
                'long_axis',long_axis,...
                'short_axis',short_axis,...
                'status','' );
        else
            % report an empty structure
            ellipse_t = struct( ...
                'a',[],...
                'b',[],...
                'phi',[],...
                'X0',[],...
                'Y0',[],...
                'X0_in',[],...
                'Y0_in',[],...
                'long_axis',[],...
                'short_axis',[],...
                'status',status );
%             X0=0; Y0=0;
        end

        % check if we need to plot an ellipse with it's axes.
        if (test>0)%(nargin>2) & ~isempty( axis_handle ) & (test>0)

            % rotation matrix to rotate the axes with respect to an angle phi
            R = [ cos_phi sin_phi; -sin_phi cos_phi ];

            % the axes
        % %     ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
        % %     horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
        % %     new_ver_line    = R*ver_line;
        % %     new_horz_line   = R*horz_line;

            % the ellipse
            theta_r         = 0:2*pi/brPodelaZaRezTacke:2*pi;theta_r=theta_r(1:end-1);
            ellipse_x_r     = X0 + a*cos( theta_r );
            ellipse_y_r     = Y0 + b*sin( theta_r );
            rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
            rezTacke = rotated_ellipse; rezTacke(:,3)=0;
        %     rotated_ellipse=rotated_ellipse';
            % draw
        %     hold_state = get( axis_handle,'NextPlot' );
        %     set( axis_handle,'NextPlot','add' );
        %     plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
        %     plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
%             hold on; plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
        %     set( axis_handle,'NextPlot',hold_state );
        end
        end
        %%
        function [rez]=getElipsuCentroidHekse(Cvorovi,Hekse)
            rez = Cvorovi(Hekse(:,1),:);
            for i=2:8
                rez = rez+ Cvorovi(Hekse(:,i),:);
            end
            rez = rez/8;
        end
        %% da li je tacka unutar hekse(za refinement heksi kod arsStent)
        %inputs 
            %Cvorovi [ tacke 
            %Tetra - JEDNA ILI VISE!!! enumeracija cvorova krece od 1!
            %tacka - SAMO JEDNA!!!!!!!
        function [rez] = daLiSuTackeUnutarTetre(Cvorovi, Tetra, Tacke)  
%             %ref\http://steve.hollasch.net/cgindex/geometry/ptintet.html
%             for i = 1:numel(Tacke(:,1))
%                 if i == 34953
%                     ars0=111;
%                 end
%                 Tacka = Tacke(i,:);
%                 d0 = Cvorovi(Tetra,:)     ;
%                 a  = [1 1 1 1]            ;
%                 d0 = [d0 a(:)]            ;
%                 D0 = det(d0)              ;
%                 % determinanta D1
%                 d1 = d0; d1(1,1:3) = Tacka;
%                 D1 = det(d1)              ;
%                 % determinanta D2
%                 d2 = d0; d2(2,1:3) = Tacka;
%                 D2 = det(d2)              ;
%                 % determinanta D3
%                 d3 = d0; d3(3,1:3) = Tacka;
%                 D3 = det(d3)              ;
%                 % determinanta D4
%                 d4 = d0; d4(4,1:3) = Tacka;
%                 D4 = det(d4)              ;
%                 if abs(D0-(D1+D2+D3+D4))<0.000001%provera
%                     %If any other Di=0, then P lies on boundary i (boundary i being that boundary formed by the three points other than Vi). 
%                     if (D1*D2*D3*D4) == 0
% %                         rez(i)=i; continue; -pogresan uslov!boundary==na
% %                         nekoj od linija/pravaca koje cine tetru
%                     end
%                     if D0 >0%svi treba da budu 0 ili veci od nule
%                         if D1 >=0 & D2 >=0 & D3 >=0 & D4 >=0
%                             rez(i)=i; continue;
%                         else
%                             rez(i)=0;continue;
%                         end
%                     else
%                         if D1 <=0 & D2 <=0 & D3 <=0 & D4 <=0
%                             rez(i)=i;continue;
%                         else
%                             rez(i)=0;continue;
%                         end
%                     end
%                 else
%                     if (D1*D2*D3*D4) == 0
%                         rez(i)=i; continue;
%                     else
%                         rez(i) = -1;
%                     end
%                 end
%             end%for tetre
%             rez = rez(:);
           rez = daLiSuTackeUnutarTetre(Cvorovi,Tetra-1,Tacke);
           rez = rez>0;
        end
        %%
        %function returs stlPoints fromat and ABC format if its needed,if not - just delete it and adopt to your needs
        function [stlPoints Apoints Bpoints Cpoints ]= sphereTriangulation(numIterations, radius)
            if nargin <1    
                radius        = 1; % unit circle
                numIterations = 3; % enough for testing
            elseif nargin<2 %set default    
                radius        = 1;
            end%no more error handlers from me since I am not supposed to know for what You will use it :)
            %basic Octahedron reff:http://en.wikipedia.org/wiki/Octahedron
            %     ( ?1, 0, 0 )
            %     ( 0, ?1, 0 )
            %     ( 0, 0, ?1 )
            A=[1 0 0]*radius;
            B=[0 1 0]*radius;
            C=[0 0 1]*radius;
            %from +-ABC crete initial triangles which define oxahedron
            triangles=[ A;  B;  C;...
                        A;  B; -C;....
                      % -x, +y, +-Z quadrant
                       -A;  B;  C;...
                       -A;  B; -C;...
                      % -x, -y, +-Z quadrant
                       -A; -B;  C;...
                       -A; -B; -C;...
                      % +x, -y, +-Z quadrant
                        A; -B;  C;...
                        A; -B; -C];%-----STL-similar format 
             %for simplicity lets break into ABC points...
             selector         = 1:3:numel(triangles(:,1))-1 ;
             Apoints = triangles(selector  ,:)    ;
             Bpoints = triangles(selector+1,:)    ;
             Cpoints = triangles(selector+2,:)    ;
            %in every of numIterations 
            for iteration = 1:numIterations
                %devide every of triangle on three new
                %        ^ C
                %       / \
                % AC/2 /_4_\CB/2
                %     /\ 3 /\
                %    / 1\ /2 \
                % A /____V____\B           1st              2nd              3rd               4th        
                %        AB/2
                %new triangleSteck is [ A AB/2 AC/2;     AB/2 B CB/2;     AC/2 AB/2 CB/2    AC/2 CB/2 C]
                AB_2   =  (Apoints+Bpoints)/2;
                %do normalization of vector 
                AB_2   =  AngioIvusMath.arsUnit(AB_2)*radius;...same for next 2 lines
                AC_2   =  (Apoints+Cpoints)/2; AC_2   =  AngioIvusMath.arsUnit(AC_2)*radius;
                CB_2   =  (Cpoints+Bpoints)/2; CB_2   =  AngioIvusMath.arsUnit(CB_2)*radius;
                Apoints = [ Apoints;... A point from 1st triangle
                            AB_2;...    A point from 2nd triangle
                            AC_2;...    A point from 3rd triangle 
                            AC_2];...   A point from 4th triangle..same for B and C
                Bpoints = [AB_2; Bpoints; AB_2; CB_2   ];
                Cpoints = [AC_2; CB_2   ; CB_2; Cpoints];
            end
            %now tur points back to STL-like format....
            numPoints = numel(Apoints(:,1))                                 ; 
            selector  = 1:numPoints                                         ;
            selector  = selector(:)                                         ;    
            selector  = [selector, selector+numPoints, selector+2*numPoints];    
            selector  = selector'                                           ;
            selector  = selector(:)                                         ;%thats it
            stlPoints = [Apoints; Bpoints; Cpoints]                         ;
            stlPoints = stlPoints(selector,:)                               ;
        end
        %%
        %INPUTS
            %jedna! HEKSA
            %mnogo tacki
        function [idTackiKojeSuUHeksi] = daLiSuTackeUnutarHekse(Cvorovi,Tacke)    
            Tetre = DelaunayTri(Cvorovi(:,1),Cvorovi(:,2),Cvorovi(:,3));
            TetreId = Tetre.Triangulation;
            for i = 1:numel(TetreId(:,1))
                daLiJe(:,i) = AngioIvusMath.daLiSuTackeUnutarTetre(Cvorovi, TetreId(i,:), Tacke);
            end      
            try
                id                  = unique(daLiJe);
                idTackiKojeSuUHeksi = id(3:end);
            catch
                idTackiKojeSuUHeksi = 0;
            end
        end
        %% 
        %INPUTS
            %Cvorovi          - koordinate cvorova
            %Elementi         - hekse
            %RefinementRegion - definisan oblakom tacaka!(posle se vrsi mesh tog oblaka i tako se definise region)
        function [idHeksiUnutarRefinementRegiona] = getHekseUnutarOblakaTacaka(Cvorovi,Elementi,RefinementRegion)
            % #1 nadji tezista svih elemenata
            tezistaHeksi = (Cvorovi(Elementi(:,1),:)+.....
                            Cvorovi(Elementi(:,2),:)+.....
                            Cvorovi(Elementi(:,3),:)+.....
                            Cvorovi(Elementi(:,4),:)+.....
                            Cvorovi(Elementi(:,5),:)+.....
                            Cvorovi(Elementi(:,6),:)+.....
                            Cvorovi(Elementi(:,7),:)+.....
                            Cvorovi(Elementi(:,8),:))/8;
            % #2 nadji hekse koje su untar refinement regiona
            upisiVTK_HEXE(RefinementRegion,[1 2 3 4 5 6 7 8],'arsStent.RefinementBoxa.vtk'); 
            [idHeksiUnutarRefinementRegiona] = AngioIvusMath.daLiSuTackeUnutarHekse(RefinementRegion, tezistaHeksi);
            upisiVTK_HEXE(Cvorovi,Elementi(idHeksiUnutarRefinementRegiona,:),'arsStent.RefinementHexeUnutarRefinementBoxa.vtk');
        end
        %% presek Linije(Po,P1) kroz Piramidu(ABCD-osnova, i F-vrh)
        % INPUTS: 
            %Linija(Po,P1)
            %Piramida: osnova(ABCD), vrh(F)
        % Outputs: tacke A i B
        function [rez] = getTackePreseka3DPiramideI3DLinije(Po,P1, A,B,C,D,F)
            debugMode=0;
            %Piramida je napravljena od 4 trougla i cetvorougla=osnove
            % osnova se deli na dva trougla, dakle ima 6 trouglova  
            %  A.....B
            %   :   :
            %  C.....D
            ABCD = [ A; B; D; C; A];
            TackePreseka = [666,666,666];%pomocna inicijaliacija
                %trazi presecne tacke po omotacu
                 for i=1:4
                     T = AngioIvusMath.getTackuPresekaLinijeP1P2iTrouglaABC(Po,P1,F,ABCD(i,:,:),ABCD(i+1,:,:));
                        if(~strcmp(T,'error'))
                         TackePreseka = [ TackePreseka; T ];
                     end
                 end                 
                 TackePreseka = TackePreseka(2:end,:,:);
                 TackePreseka = unique(TackePreseka, 'rows');
                  %  A.....B
                  %   :   :
                  %  C.....D
                 CBD = [B;D;C];
                 for i=1:2
                     T = AngioIvusMath.getTackuPresekaLinijeP1P2iTrouglaABC(Po,P1,A,CBD(i,:,:),CBD(i+1,:,:));
                     if(~strcmp(T,'error'))
                         TackePreseka = [ TackePreseka; T ];
                     end
                 end
                 TackePreseka = unique(TackePreseka, 'rows');
                 %trazi presecne tacke po osnovi
                 if (numel(TackePreseka(:,1,:))==2)%ako postoje dve presecne tacke
                    rez = TackePreseka;
                 else %ne sece ABCD u dve tacke
                    rez = 'error'; 
                 end                   
                 if debugMode
                  AngioIvusMath.nacrtajPiramidu(A,B,C,D,F);
                  if(~strcmp(rez,'error'));
                    plotLine(TackePreseka(1,:,:),TackePreseka(2,:,:));
                  else                      
                    plotLine(Po,P1);
                  end
                 end%debugMode
        end%getTackePreseka3DPiramideI3DLinije
        %% Projektuje epipolarnu linijju(PoP1) na ravan slike ABCD iz
        %% fokalne tacke F
        function [rez] = getProjekcijuEpipolarneLinjePoP1naSlikuABCDizF(Po,P1,A,B,C,D,F)
            tackePresekaLijePoP1saPiramidomABCDF = AngioIvusMath.getTackePreseka3DPiramideI3DLinije(Po,P1,A,B,C,D,F);
            %da li linija sece piramidu???(da li se epipolarna linija vidi
            %na slici????)
            if(~strcmp(tackePresekaLijePoP1saPiramidomABCDF,'error'))
            T1 = tackePresekaLijePoP1saPiramidomABCDF(1,:,:);
            T2 = tackePresekaLijePoP1saPiramidomABCDF(2,:,:);
            %projektuj tacke T1 i T2 na sliku ABCD iz fokalnog izvora F
            vektNormaleRavniABCD = cross(B-A, C-A);%vektor normale na ravan slike ABCD
            projekcijaTackeT1 = AngioIvusMath.getTackuPresekaLinijeIRavni(T1,F,A,vektNormaleRavniABCD);
            projekcijaTackeT2 = AngioIvusMath.getTackuPresekaLinijeIRavni(T2,F,A,vektNormaleRavniABCD);
            rez = [projekcijaTackeT1; projekcijaTackeT2];
            else%linija PoP1 ne sece piramidu ABCDF
                rez = 'error';
            end%if
        end
        %% presecne tacke B-Splajna i prave linije PoP1
        % INPUTS:
            % BSplajn = otvoreni/zatvoreni splajn
            % Linija  = PoP1
        function [rez paramT] = getPresecneTackeBSplajnaILinijePoP1( BSplajn, Po,P1)
            %interpolira se Splajn od t=0 do t=1 sa odredjenim pomerajem 
            %i trazi se presek pravog segmenta sa linijom
            paramT=0;
            T1 = BSplajn.interpolate(BSplajn,0);
            rez = [666,666,666];
            for i=0.02:0.002:1
                T2 = BSplajn.interpolate(BSplajn,i); 
                tackaPreseka = AngioIvusMath.getTackuPresekaLinijaNaSegmentima(...
                                                             T1,T2, Po,P1);
                %ako postoji tacka preseka
                if(~strcmp(tackaPreseka,'error'))
                    rez = [rez;tackaPreseka];
                    paramT=[paramT;i];
                end%if
                T1 = T2;                
            end%for
            %ako ima vise od 1 rezultata znaci da ima i preseka
            if(numel(rez)>3)
                rez    = rez(2:end,:,:);
                paramT = paramT(2:end);
                %za slucaj da je tacka preseka na sastavu segmenta
                rez    = unique(rez, 'rows');
                paramT = unique(paramT);
            else
                rez    = 'error';
                paramT = 'error';
            end             
        end%
        %% projekcija vektora A na vektor B
        %    A / 
        %     / )teta
        %   O.---------> B
        % Vektori A i B krecu iz iste tacke! ako se pocetak A ne poklapa sa
        % pocetkom vekt. B treba projektovati pocetak(tacku) vektora A na B
        % i traslirati vektor A na tu tacku pa primeniti ovu metodu
        %INPUTS:
            % A(x,y,z), B(x,y,z)
        %OUTPUTS:
            % rez(x,y,z)
        function [rez] = getProjekcijuVektoraAnaVektorB(A,B)
        %reff: http://en.wikipedia.org/wiki/Vector_projection
        %Neka su A i B dva vektora. Projekcija vektora A na vektor B je
        %vektor C koji ima isti pravac kao B i intenzitet |a|cos(teta). Gde
        %je cos(teta)= dot(a,b) / ( |a|*|b|)
        duzinaVektA =  AngioIvusMath.arsNorm(A);
        kosinusTeta =  AngioIvusMath.arsDot(A,B) ./ (AngioIvusMath.arsNorm(A).*AngioIvusMath.arsNorm(B));
        B=AngioIvusMath.vektorizuj(B,numel(kosinusTeta));
        rez         =  B .* [duzinaVektA duzinaVektA duzinaVektA] .* [kosinusTeta kosinusTeta kosinusTeta];        
        end%getProjekcijuVektoraAnaVektorB
        %% projekcija Tacke A na Liniju PoP1
        %          A
        %   Po .__________. P1
        %INPUTS:
            %Tacka u 3D: A(x,y,z)
        %OUTPUTS:
            %Linija u 3D: Po(x,y,z),P1(x,y,z)
        function [rez] = getProjekcijuTackeAnaLinijuPoP1(A,Po,P1)
        %reff: http://cs.nyu.edu/~yap/classes/visual/03s/hw/h2/math.pdf
        %racuna se vektor PoA i zatim njegova projekcija na vektor PoP1
        if (numel(A(1,:))==1)
            A=A'; Po=Po'; P1=P1';
        end
        APo                     = AngioIvusMath.arsMinus(A,Po); 
        PoP1                    = AngioIvusMath.arsMinus(P1,Po);
        vektProjekcijeAPoNaPoP1 = AngioIvusMath.getProjekcijuVektoraAnaVektorB(APo,PoP1);
        rez                     = AngioIvusMath.arsPlus(vektProjekcijeAPoNaPoP1, Po);
        end%getProjekcijuTackeAnaLinijuPoP1
        %% projektuje tacku na osu i vraca polozaj tacke duz ose. 
        %INPUTS
            %Tacka - niz tacaka koje se konvertuju
            %O     - koordinatni pocetak(pocetak linije)
            %V     - vektor pravca
        %OUTPUTS
            %rezTacka - 
            %tRezTacka - parametar duzine(sa znakom!)
        function [rezTacka tRezTacka] = getProjekcijuTackeNaOsu(Tacka, O,V)
            V            = AngioIvusMath.arsUnit(V);
            TackaOV      = AngioIvusMath.getProjekcijuTackeAnaLinijuPoP1(Tacka, O, O+V);
            tRezTacka    = AngioIvusMath.arsMinus(TackaOV,O);
            tRezTacka    = AngioIvusMath.arsNorm(tRezTacka);
            znakRezTacka = AngioIvusMath.arsDot(AngioIvusMath.arsMinus(TackaOV,O),V);
            znakRezTacka = znakRezTacka./abs(znakRezTacka);            
            tRezTacka    = tRezTacka .*znakRezTacka;
            O=AngioIvusMath.vektorizuj(O,numel(tRezTacka));
            V=AngioIvusMath.vektorizuj(V,numel(tRezTacka));
            rezTacka= O + V .* [tRezTacka tRezTacka tRezTacka]; 
        end
        %% projekcija tacke na ravan
        %INPUTS:
            %T- tacka koja se projektuje
            %Ravan:
                %A - tacka u ravni
                %n - vektor normale
        %OUTPUTS: 
            %tacka u ravni
        function [rez] = getProjekcijuTackeAnaRavan(T,TackaRavni,normalaRavni)
            %svodi se na presek linije i ravni. U tacki T se ispaljuje zrak sa pravcem normale ravni
            brTacki=numel(T(:,1));
            TackaRavni   = AngioIvusMath.vektorizuj(TackaRavni,brTacki);
            normalaRavni = AngioIvusMath.vektorizuj(normalaRavni,brTacki);
            rez = AngioIvusMath.getTackuPresekaLinijeIRavniParalelizovano(T, T+normalaRavni, TackaRavni, normalaRavni);
        end%getProjekcijuTackeAnaLinijuPoP1
        %% funkcije za graficki prikaz
        function nacrtajTrougao(A,B,C)     
            for i = 1:numel(A(:,1))
                PlotLine(A,B); PlotLine(C,B); PlotLine(C,A);
            end
            axis equal;
        end%
        %%
        function nacrtajCetvorougao(A,B,C,D)
            if nargin == 1
                pom = A;
                A   = pom(1,:); B   = pom(2,:); C   = pom(3,:); D   = pom(4,:);
            end
            plotLine([A;B]); plotLine([B;D]); plotLine([C;D]); plotLine([C;A]);
            axis equal;
        end%
        %%
        function nacrtajQuadrilaterial(A,B,C,D)
            if nargin==2
                Cvorovi        = A;
                Quadrialateriali=B;
                A=Cvorovi(Quadrialateriali(:,1),:);
                B=Cvorovi(Quadrialateriali(:,2),:);
                C=Cvorovi(Quadrialateriali(:,3),:);
                D=Cvorovi(Quadrialateriali(:,4),:);
            end
         for i=1:numel(A(:,1))
             AngioIvusMath.nacrtajCetvorougao(A(i,:), B(i,:), D(i,:), C(i,:));
         end
        end
        %%
        function nacrtajHexaedar(A,B,C,D, E,F,G,H)
            if nargin ==1
                pom = A;
                A = pom(1,:,:);
                B = pom(2,:,:);
                C = pom(3,:,:);
                D = pom(4,:,:);
                E = pom(5,:,:);
                F = pom(6,:,:);
                G = pom(7,:,:);
                H = pom(8,:,:);
            end
            AngioIvusMath.nacrtajCetvorougao(A,B,D,C);
            AngioIvusMath.nacrtajCetvorougao(A,B,E,F);
            AngioIvusMath.nacrtajCetvorougao(C,B,G,F);%3
            AngioIvusMath.nacrtajCetvorougao(G,F,H,E);%4
            AngioIvusMath.nacrtajCetvorougao(A,D,E,H);%5
            AngioIvusMath.nacrtajCetvorougao(D,C,H,G);
            axis equal;
        end%
        %%
        function nacrtajPiramidu(A,B,C,D,F)
            F = F(:)';
            AngioIvusMath.nacrtajCetvorougao(A,B,C,D);
            plotLine([F;A]); 
            plotLine([F;B]); 
            plotLine([F;C]); 
            plotLine([F;D]);
            axis equal;
            xlabel('X osa');ylabel('Y osa');zlabel('Z osa');
        end%
        %% konvertuj ugao u radijane
        function [rad] = konvertujUgaoUUradijane(ugao)
            rad = (ugao*pi)/180;
        end
        %% konvertuj ugao u radijane
        function [ugao] = konvertujRadijaneUUgao(radijan)
            ugao = (radijan*180)/pi;
        end
        %% konvertuj binarni broj u doubl
        function [rez] = bin2Double(binarniBroj)
            rez = 0;
            brBinarnihBrojeva = numel(binarniBroj);
            for trenutnaPozicijaBinarnogBroja=1:brBinarnihBrojeva
               rez = rez + 2^(trenutnaPozicijaBinarnogBroja-1) * binarniBroj(brBinarnihBrojeva - trenutnaPozicijaBinarnogBroja +1.);
            end
            
        end
       %% nadji pikove 2D signala-nalazi SVE preveojne tacke
        % INPUTS - 2Dsignal
        %   . maxPik
        %  / \   /
        % /   \./minPik
        function [rez] = nadjiPrevojneTacke2DSignala(signal2D)
            signalA            = signal2D;
            signalB            = [0; signal2D(1:end-1)];
            razlika            = signalA - signalB;
            razlika(razlika>0) =  1; %pozitivan, prvi izvod>0, rastuca
            razlika(razlika<0) = -1;%opadajuca fja
            %nadji pikove, tamo gde fja menja znak
            signalA            = razlika;
            signalB            = [0; razlika(1:end-1)];
            razlika(razlika>0) =  1;%maxPik
            razlika(razlika<0) = -1;%minPik, gde je = 0, tu nema pika.
            razlika            = signalA - signalB;
            maxPik             = find(razlika ==-2);
            minPik             = find(razlika == 2);
            rez.maxPik         = maxPik-1;%korekcija za 1 polje jer selektuje prvi sledeci a ne mesto gde se desio pik
            rez.minPik         = minPik-1;
        end
        %% 
        function [rez]  = usrednjavanje1DSignala(signal1D, korak,debug)
            if nargin == 1
                korak = 3;
                debug = 0;
            elseif nargin==2
                debug = 0;
            end
            brElemenata = numel(signal1D);
            rez = signal1D;
            for i =1+korak:brElemenata-korak
                rez(i)= mean(signal1D(i-korak:i+korak));
            end
            if debug
                plot(signal1D);
                hold on;
                plot(rez,'color','r');
            end
        end
        %
        function [rez] = usrednjavanje1DSignalaIOdsecanjeVelikihCurvatura(signal1D, korak,debug)
            if(numel(signal1D(1,:))==2)
                xBrKolone = signal1D(:,2);
                signal1D  = signal1D(:,1);
            end
            singalUsrednjen = AngioIvusMath.usrednjavanje1DSignala(signal1D, korak,debug);
            pom = abs(signal1D - singalUsrednjen);
            pom(pom>5)=666;
            pom = AngioIvusMath.usrednjavanje1DSignala(pom, korak,debug)
            pom(pom<100)=1;
            pom(pom>100)=0;
%             pom = signal1D .* pom;
            rez = [0 0];
            for i = 1:numel(pom)
                if pom(i)
                    rez = [ rez; xBrKolone(i) signal1D(i)];
                end
            end
            rez = rez(2:end,:);
        end
        %%  Po - prva tacka, P1 - druga tacka, dt - br podela
        function [rez] = interpolirajRastojanjeIzmedju2Tacke(Po,P1,dt)
            if(dt==1)
                rez = [Po;P1];
            else
            rez = Po;
            rezX = linspace(Po(1), P1(1), dt+1);
            rezY = linspace(Po(2), P1(2), dt+1);
            rezZ = linspace(Po(3), P1(3), dt+1);
            rez  = [rezX', rezY', rezZ'];
            end
        end
        %%
        function [rez] = interpolirajSplajnom(Tacke, brPodela)
            %2D
            if numel(Tacke(1,:))==2
                Tacke1 = [Tacke Tacke(:,1)];
                Tacke1(:,3)=0;
            end
            splajn = arsBSplineOpenNonuniform;
            splajn = splajn.set(splajn, Tacke1, 3, 66) ;
            noveTacke =  splajn.interpolateFromAtoB(splajn, 0, 1/brPodela, 1);
            if numel(Tacke(1,:))==2
                rez = noveTacke(:,[1 2]);
            end
        end
        %% Funkcija vraca (lokalni) deo (globalnog) splajna. 
        %napraljeno za potrebe patient-specific savijanja stenta
        %INPUTS
            %tackeSplajna                      - Tacke koje opisuju splajnkoji se posmatra
            %duzinaPosmatranogDela             - Duzina dela splajna koji se posmatra
            %brTacakaPosmatranogDela           - Br tacaka sa koliko se interpolira posmatrani deo
            %tPocetakPosmatranogDelaDuzSplajna - Pozicija posmatranog dela duz splajna
        %REZ
            %rezTacke - posmatrani deo interpoliran sa brTacakaPosmatranogDela
        function [rezTacke] = interpolirajSplajnNaOdredjeniBrPodelaPocevsiOdT(tackeSplajna, duzinaPosmatranogDela, brTacakaPosmatranogDela, tPocetakPosmatranogDelaDuzSplajna)
        %Nadji tangentu na prvu tacku %%1 Napravi B-Spline od kontrolsnih tacaka putanje i napravi 
            %%1 Info o stentu
                stent.brTacakaStenta = brTacakaPosmatranogDela; % broj podela stenta po duzini
                stent.duzinaStenta   = duzinaPosmatranogDela  ; % Duzina stenta u mm
            %%2 Parametrizuj putanju (deo arterije) na kojoj je stent
                arterija.Tacke3Dt                     = tackeSplajna                                                                                    ;
                arterija.splajn                       = arsCubicBSplineOpenNonuniform                                                                   ;
                arterija.splajn                       = arterija.splajn.set(arterija.splajn, arterija.Tacke3Dt, 3)                                      ;
                arterija.duzina                       = AngioIvusMath.getDuzinuLinije(arterija.splajn.interpolateFromAtoB(arterija.splajn, 0, 0.005, 1)); % Ukupna duzina arterije
                arterija.tPocetakStentaUnutarArterija = tPocetakPosmatranogDelaDuzSplajna                                                               ;
                arterija.tKrajStentaUnutarArterije    = arterija.tPocetakStentaUnutarArterija + stent.duzinaStenta/arterija.duzina                      ; % Tacka duz arterije gde se stent zavrsava
                arterija.koraciInterpolacije          = [arterija.tPocetakStentaUnutarArterija:(arterija.tKrajStentaUnutarArterije-arterija.tPocetakStentaUnutarArterija)/(stent.brTacakaStenta-1):arterija.tKrajStentaUnutarArterije];
            %%3 Interpoliraj deo Arterije sa stentom na zadati broj tacki
                rezTacke                              = arterija.splajn.testQuaternionOrjentacijuPreseka(.....
                                                                           arterija.splajn,................
                                                                           arterija.Tacke3Dt,..............
                                                                           3,..............................
                                                                           arterija.koraciInterpolacije,...
                                                                           0.001);.........................
        end%interpolirajSplajnNaOdredjeniBrPodelaPocevsiOdT
        %%
        function [rez] = interpolirajRastojanjeIzmedju2TackeVratiTackuNaT(Po,P1,t)
            v = AngioIvusMath.arsMinus(P1,Po);
            d = AngioIvusMath.arsNorm(v);
            d = [d d d];
            v = v./d; 
            v = v .* d .* [t(:) t(:) t(:)];
            rez = Po + v;
        end
        %% duzina linije odredjene nizom tacaka
        function [rez]= izracunajDuzinuPolilajna(nizTacaka)
            brTacaka = numel(nizTacaka(:,1,:))-1;
            rez=0;
            for i =1:brTacaka
                rez=rez+ norm(nizTacaka(i+1,:,:)-nizTacaka(i,:,:));
            end
        end
        %% izbrisi IDove koji su blizu | za rekonstrukciju 3D ivusa
        function [rez] = izbrisiVisakIDova(id,dmin)
            if nargin ==1
                dmin=20;%minimalni razmak izmedju 2 frejma | rad srca 20~30
            end
            rez = id(1);
            for i=2:numel(id)-1
                if dmin<(id(i)-rez(end))
                    rez=[rez; id(i)];
                end
            end
            if dmin<(id(end)-rez(end))
            	rez=[rez; id(end)];
            else
                rez=[rez(1:end-1); id(end)];
            end
        end
        %% Prebaci tacke iz globalnog u lokalni koordinatni sistem
        %INPUTS 
            %Tacke-tacke koje se prebacuju
            %Koordinatni sistem 2: O2, X2, Y2, Z2
        function [rezTacke rezTackeX rezTackeY rezTackeZ]=konvertujTackeIzGlobalnogKoordinatnogSistemaULokalni(Tacke, O2,X2,Y2,Z2)
%             O2=AngioIvusMath.vektorizuj(O2,Tacke);
%             X2=AngioIvusMath.vektorizuj(X2,Tacke);
%             Y2=AngioIvusMath.vektorizuj(Y2,Tacke);
%             Z2=AngioIvusMath.vektorizuj(Z2,Tacke);
            [TackeXY]= AngioIvusMath.getProjekcijuTackeAnaRavan(Tacke, O2, Z2);
            [rezTackeX tRezTackeX] = AngioIvusMath.getProjekcijuTackeNaOsu(TackeXY, O2, X2);
            [rezTackeY tRezTackeY] = AngioIvusMath.getProjekcijuTackeNaOsu(TackeXY, O2, Y2);
            [TackeXZ]= AngioIvusMath.getProjekcijuTackeAnaRavan(Tacke, O2, Y2);
            [rezTackeZ tRezTackeZ] = AngioIvusMath.getProjekcijuTackeNaOsu(TackeXZ, O2, Z2);
            rezTacke=[tRezTackeX,tRezTackeY,tRezTackeZ];
            rezTacke(isnan(rezTacke))=0;
        end
        %%
        function [rezTacke]=konvertujTackeIzLokalnogKoordinatnogSistemaGlobalni(Tacke, O2,X2,Y2,Z2)
            brTacki = numel(Tacke(:,1));
            O2 = AngioIvusMath.vektorizuj(O2,brTacki);
            X2 = AngioIvusMath.vektorizuj(X2,brTacki);
            Y2 = AngioIvusMath.vektorizuj(Y2,brTacki);
            Z2 = AngioIvusMath.vektorizuj(Z2,brTacki);
            TackeX = Tacke(:,[1 1 1]);
            TackeY = Tacke(:,[2 2 2]);
            TackeZ = Tacke(:,[3 3 3]);
            rezTacke = O2 + TackeX .* X2 + TackeY .* Y2 + TackeZ .* Z2; 
        end
        
        %% iz oblaka tacaka izdvaja onekoji su na ravni
        %INPUTS
            %TACKA ABC (tacke koje definisu ravan)[A;B;C]
        %OUTPUTS
            %TackeKandidati - oblak tacaka
        function [CvoroviNaRavni, IdCvorovaNaRavni, RastojanjaCvorovaOdRavni] = getTackeNaRavniABC(TackeABC, TackeKandidati)
            Cvorovi             = TackeKandidati;
            A                   = TackeABC(1,:);
            B                   = TackeABC(2,:);
            C                   = TackeABC(3,:);
            BminusA             = B-A;
            CminusA             = C-A;
            normalaRavni        = cross(CminusA,BminusA);
            normalaRavni        = normalaRavni / norm(normalaRavni);
            normalaCvorovi      = Cvorovi;
            normalaCvorovi( :,1)= normalaRavni(1);
            normalaCvorovi(:,2) = normalaRavni(2);
            normalaCvorovi(:,3) = normalaRavni(3);
            tackaP1Cvorovi      = Cvorovi;
            tackaP1Cvorovi(:,1) = A(1);
            tackaP1Cvorovi(:,2) = A(2);
            tackaP1Cvorovi(:,3) = A(3);
            %%%%jna ravni
            f=(Cvorovi-tackaP1Cvorovi);
            jnaRavni          = normalaCvorovi .* f;
            rez               = jnaRavni(:,1)+jnaRavni(:,2)+jnaRavni(:,3);
            IdCvorovaNaRavni  = find(abs(rez)<0.003); %oni koji su mnogo blizu ravni
            CvoroviNaRavni    = Cvorovi(IdCvorovaNaRavni);
            RastojanjaCvorovaOdRavni  =rez;
        end
        %% get rasojanje izmedju ravni i tacke
        function rastojanja = getRastojanjeTackeOdRavni(Tacka, TackaRavni, VektNormaleRavni)
            %vektorizacija
            pom = Tacka;  pom(:,1) = TackaRavni(1)      ; pom(:,2) = TackaRavni(2)      ; pom(:,3) = TackaRavni(3)      ; TackaRavni       = pom;
            pom = Tacka;  pom(:,1) = VektNormaleRavni(1); pom(:,2) = VektNormaleRavni(2); pom(:,3) = VektNormaleRavni(3); VektNormaleRavni = pom;
            %reff:http://geomalgorithms.com/a04-_planes.html
            sn = -AngioIvusMath.arsDot(VektNormaleRavni,(Tacka - TackaRavni));
            sd =  AngioIvusMath.arsDot(VektNormaleRavni,VektNormaleRavni);
            sb = sn ./ sd;
            B  = Tacka + [sb sb sb] .* VektNormaleRavni;
            rastojanja = AngioIvusMath.arsNorm(Tacka-B);
        end
        %% 
        function [rez]=getPorvrsinuQuadrilateriala(Cvorovi, Quadrilateriali)
            [Cvorovi Trouglovi] = AngioIvusMath.konvertujQuadrilaterialeToTrouglove(Cvorovi, Quadrilateriali);
            upisiSTL(Cvorovi,Trouglovi,'QuadrilaterialiIzkojihSeRacunaPovrsina.stl');
            idCvoroviTrouglovi  = Trouglovi';             
            CvoroviTrouglovi    = Cvorovi(idCvoroviTrouglovi(:),:);
            rez                 = AngioIvusMath.getPovrsinuStl(CvoroviTrouglovi);
        end
        %%        .C
        %        / \  
        %     b /    \ a
        %      /   c   \
        %   A .----------.B
        function [rez]=getPovrsinuStl(stlPath)
            if isstr(stlPath)
                [cvoroviZaOgranicavanje] = arsVTK.ucitajTrougloveTackeizSTLa(stlPath);
            else
                cvoroviZaOgranicavanje   = stlPath;
            end
            %cvoroviZaOgranicavanje=[A1; B1; C1; A2; B2 C2;...itd]
            A=cvoroviZaOgranicavanje(1:3:end,:);%prva pa svaka treca
            B=cvoroviZaOgranicavanje(2:3:end,:);
            C=cvoroviZaOgranicavanje(3:3:end,:);
            %---------------------------------------------
            a = AngioIvusMath.arsNorm(C-B);
            b = AngioIvusMath.arsNorm(C-A);
            c = AngioIvusMath.arsNorm(A-B);
            %Area = sqrt(s*(s-a)*(s-b)*(s-c)) where s=(a+b+c)/2 is the semi-perimeter.
            s = (a+b+c)/ 2;
            rez = sqrt(s .* (s-a) .* (s-b) .* (s-c)); 
            rez=sum(rez);
        end
        %% preracunava silu u pritisak 
        function [pritisak_MPa] =  preracunajSiluUPritisakKojiTrebaZadatiNaFejsove(stlPath,silaNjutni)
                povrsinaSTLa        = AngioIvusMath.getPovrsinuStl(stlPath);% [mm2]
                povrsinaSTLa_metar2 = povrsinaSTLa * 0.000001              ;% [m2 ]
                pritisak            = silaNjutni/povrsinaSTLa_metar2;% [Pa ]
                pritisak_MPa        = pritisak * 0.000001                  ;% [MPa]
        end            
        %% dot(A,B)= |A||B|cos(A,B)
        function [rez]=rgetUgaoIzmedju2Vektora(V1,V2)
            cosTeta = AngioIvusMath.arsDot(V1,V2) ./ (AngioIvusMath.arsNorm(V1) .* AngioIvusMath.arsNorm(V2) );
            teta = acos(cosTeta);
            rez  = rad2deg(teta);
        end
        %%
        function rez=getDuzinuLinije(tacke)
            rez = tacke(2:end,:) - tacke(1:end-1,:);
            rez = rez .* rez;
            rez = sum(rez');
            rez = sqrt(rez);
            rez = sum(rez);
        end
        %% racuna povrsinu poligona odedjenog tackama u ravni ili ravni u 3D
        % tacke poligona ne moraju biti sortirane, ali moraju biti u ravni i
        % moraju opisivati konveksan poligon(npr obod vene)
        %tacke su SORTIRANE
        %-napravljeno za racunanje povrsine slajsa arterije-pogledaj
        %funckiju AngioIvusMath.getPresekTrouglaIRavni i arsVTK.stlCenterlineSlicer
        function [rezArea]=getPovrsinuPoligonaURavni(tackePoligona,normala)            
            if nargin<2%ako  nije prosledjena normala, racuna kao 2D 
                normala = [ 0 0 1];
            end
            %normalizuj normalu(unit)
            normala             = normala/(sum(sqrt(normala.*normala)));
            %ref http://paulbourke.net/geometry/area3d/
            % Area of a arbitrary planer polygon--^^^^
            %    ._________________.
            %    |                 |_____. P4
            %    |                       |
            %Pn-1.____.Pn      P2._______.P3
            %         |         /
            %     Po .|_______./ P1
            %
            % Po=Po  Pi-1 = 1:n-1 Pi=2:n, ako bi indexi isli on 0, dakle +1
            P_i_minus_1               = tackePoligona(2:end-1,:)         ;
            Pi                        = tackePoligona(3:end,:)           ;
            Po                        = Pi                               ;            
            Po(:,1)                   = tackePoligona(1,1)               ;
            Po(:,2)                   = tackePoligona(1,2)               ;
            Po(:,3)                   = tackePoligona(1,3)               ;
            %nadji povrsinice trouglica PoP1P2 PoP2P3...PoPn-1Pn
            Area                      = AngioIvusMath.arsCross((P_i_minus_1-Po),(Pi-Po))/2;
            bckpArea                  = Area                             ;            
            Area                      = Area .* Area                     ;
            Area                      = Area(:,1)+Area(:,2)+Area(:,3)    ;
            Area                      = sqrt(Area)                       ;
            bckpArea                  = bckpArea ./ [Area Area Area]     ;
            bckpNormalaRavni          = bckpArea                         ;
            bckpNormalaRavni(:,1)     = normala(1)                       ;
            bckpNormalaRavni(:,2)     = normala(2)                       ;
            bckpNormalaRavni(:,3)     = normala(3)                       ;
            normale                   = bckpNormalaRavni; clear('bckpNormalaRavni');
            bckpArea(isnan(bckpArea)) = 0                                ;
            znakPovrsine              = AngioIvusMath.arsDot(bckpArea,normale);
            rezArea                   = sum(Area.*znakPovrsine)          ;
        end
        %% nadji povrsinu fejsa
        function [rez]=getPovrsinuFejsa(tackeFejsa)
            %  A      a        B
            %   ._____________.
            %  d|             | b
            %   |_____________.C
            %   D      c
            a = norm(tackeFejsa(2,:)-tackeFejsa(1,:));
            b = norm(tackeFejsa(3,:)-tackeFejsa(2,:));
            c = norm(tackeFejsa(4,:)-tackeFejsa(4,:));
            d = norm(tackeFejsa(4,:)-tackeFejsa(1,:));
            D = norm(tackeFejsa(4,:)-tackeFejsa(2,:));
            D2= norm(tackeFejsa(3,:)-tackeFejsa(1,:));
            % zbir povrsina 2 trougla
            %Area = sqrt(s*(s-a)*(s-b)*(s-c)) where s=(a+b+c)/2 is the semi-perimeter.
            s1 = (a+b+D)/2;
            A1 = sqrt(s1*(s1-a)*(s1-b)*(s1-D));
            s2 = (c+d+D2)/2;
            A2 = sqrt(s2*(s1-c)*(s2-d)*(s2-D));
            rez = A1+A2;            
        end
        %% resavanje kubicne jendacine | pravljeno za resavanje kubicne
        %% jednacine za odredjivanje glavnih napona iz invarijanti|
        %% arsSTENT-SN analiza
        % funkcija oblika a*X^3 + b*X^2 + c*X^1 + d = 0;
        function [rez] = resiKubicnuJednacinu(a,b,c,d)
            a=a(:); b=b(:); c=c(:); d=d(:);
            %ref : http://www.1728.org/cubic2.htm
            a_2 = a.*a; a_3 = a_2.*a;
            b_2 = b.*b; b_3 = b_2.*b;           
            f   = ((3*c ./ a) - (b_2 ./ a_2))/3            ;
            g   = ((2*b_3./a_3)-(9*b.*c./a_2)+(27*d./a))/27;
            h   = (g.*g/4)+(f.*f.*f/27)                    ;
%             pom = find(h
%             if h(1)>0%If h>0, there is only 1 real root and is solved by another method
                indexElemenataVecihOdNule = find(h>0);
                if numel(indexElemenataVecihOdNule)>0
                    R = -(g(indexElemenataVecihOdNule)/2)+sqrt(h(indexElemenataVecihOdNule));
                    S = nthroot(R, 3);
                    T = -(g(indexElemenataVecihOdNule)/2)-sqrt(h(indexElemenataVecihOdNule));
                    U = nthroot(T, 3); 
                    X1= (S+U)-(b(indexElemenataVecihOdNule)./(3*a(indexElemenataVecihOdNule)));
                    x1(indexElemenataVecihOdNule) = X1;%ima samo jedan relan rooot
                    x2(indexElemenataVecihOdNule) = X1;
                    x3(indexElemenataVecihOdNule) = X1;
                end
%             elseif h(1)==0 && g(1)==0 && h(1)==0
%             elseif h(1)<0% When h <= 0, as is the case here, all 3 roots are real and we proceed as follows:
                indexElemenataManjihOdNule = find(h<=0);
                if numel(indexElemenataManjihOdNule)>0
                    i = sqrt((g(indexElemenataManjihOdNule).*g(indexElemenataManjihOdNule))/4-h(indexElemenataManjihOdNule))  ;
                    j = nthroot(i, 3)     ;
                    K = acos(-(g(indexElemenataManjihOdNule)./(2*i))) ;
                    L = -j                ;
                    M = cos(K/3)          ;
                    N = sqrt(3) * sin(K/3);
                    P = -(b(indexElemenataManjihOdNule)./(3*a(indexElemenataManjihOdNule)))       ;
                    %---------------------
                    x1(indexElemenataManjihOdNule) = 2*j .* cos(K/3) - (b(indexElemenataManjihOdNule)./(3*a(indexElemenataManjihOdNule)));
                    x2(indexElemenataManjihOdNule) = L.*(M+N)+P;
                    x3(indexElemenataManjihOdNule) = L.*(M-N)+P;
                end
%             end
            rez = [x1(:) x2(:) x3(:)];
        end
         %%
         function rez=kubicnaFunkcija(x,a,b,c)
             rez= x*x*x - x*x*a + x*b - c;
         end
         %% ugao izmedju dva vektora
         function [ThetaInDegrees] = getUgaoIzmedjuDvaVektora(A,B)
             %vektorizacija
             if numel(A(:,1))>1 && numel(B(:,1))==1
                 pom = A; pom(:,1)=B(1); pom(:,2)=B(2); pom(:,3)=B(3); B = pom;
             elseif numel(A(:,1))==1 && numel(B(:,1))>1
                 pom = B; pom(:,1)=A(1); pom(:,2)=A(2); pom(:,3)=A(3); A = pom;
             end
            CosTheta       = AngioIvusMath.arsDot(A,B) ./( AngioIvusMath.arsNorm(A) .* AngioIvusMath.arsNorm(B));
            ThetaInDegrees = acos(CosTheta)*180/pi;
         end   
         %% rotiraj tacku oko linije
         %INPUTS
            %Tacka   - za rotiranje 
            %A,B     - osa rotacije
            %newUgao - ugao za koji se rotira oko ose
         %OUTPUTS
            %rez     - rotirana tacka
         function [rez] = rotirajTackuOkoLinije(Tacka, A, B, newUgao)
             quaternion     = arsQuanternion;
             newOsaRotacije = B-A;
             quaternion     = quaternion.setUgaoOsuRotacije(quaternion, newUgao, newOsaRotacije);
             for i = 1:numel(Tacka(:,1))                 
                 newTacka   = Tacka(i,:)-A;
                 newTacka   = quaternion.rotirajVektorOkoQuanterniona(newTacka,quaternion);
                 rez(i,:) = newTacka + A;
%                  plotLine([rez(i,:);newTacka]);
             end
%              plotLine(rez);
%              plotLine(Tacka);
%              plotLine( [A;B]);
         end
         %% warp path
         function [kombinacija values] = pathWarping(D)
            [N M] = size(D);
            n=1;m=1;brKombinacija=2; kombinacija=[1 1]; values=0;
            while (n <= N )
                if n>=N & m>=M
                    kombinacija(end,:)=size(D);
                    return;
                end
                while ( m <= M )
                    if n>=N & m>=M
%                         imshow(D/max(max(D)));hold on;
%                         plot(kombinacija(:,2),kombinacija(:,1),'LineWidth',2);
                        kombinacija(end,:)=size(D);
                        return;
                    elseif n==N
                        dtwKandidati       = D(n,m+1); n = n; m=m+1;
                        dtwKandidatiIndexi = [n,m+1];
                    elseif m==M
                        dtwKandidati       = [D(n+1,m)]; n = n+1; m=M;
                        dtwKandidatiIndexi = [n+1,m   ];
                    else                        
                        dtwKandidati       = [D(n+1,m+1); D(n+1,m); D(n,m+1)];
                        dtwKandidatiIndexi = [n+1,m+1   ; n+1,m   ; n,m+1   ];
                    end
                    %------
                    [value,number]               = min(dtwKandidati)           ;
                    kombinacija(brKombinacija,:) = dtwKandidatiIndexi(number,:);
                    values(number)               = value;
                    n = dtwKandidatiIndexi(number,1); m = dtwKandidatiIndexi(number,2);
                    brKombinacija                = brKombinacija+1             ;
                    [n m]
                end
            end
         end
         %%
         function [Cvorovi, Quadrilateriali] = konvertujTrougloveToQuadrilateriale(Cvorovi, Trouglovi)
             A = Cvorovi(Trouglovi(:,1),:);
             B = Cvorovi(Trouglovi(:,2),:);
             C = Cvorovi(Trouglovi(:,3),:);  ...   C
             D = (A+B+C)/3;             ...      .....
             E = (A+B)/2;               ...    G...D...F
             F = (B+C)/2;               ...   /.........
             G = (A+C)/2;               ...  A.....E.....B
%              Cvorovi = [E; D; G; A;... A; G; D; E;...
%                         D; F; C; G;... G; C; F; D;... 
%                         E; B; F; D;....D; F; B; E];
             %ovako je lakse za enumeraciju
             Cvorovi         = [E;D;E; D;F;B; G;C;F; A;G;D];
             pom             = [1:numel(A(:,1))*3]';
             Quadrilateriali = [pom pom+pom(end) pom+2*pom(end)  pom+3*pom(end)];
         end
         %%
         function [Cvorovi Trouglovi] = konvertujQuadrilaterialeToTrouglove(Cvorovi, Quadrilateriali)
             Trouglovi = [Quadrilateriali(:,[1 2 3]); Quadrilateriali(:,[1 3 4])];
         end
         %% Izocentar se nalazi na preseku linija koje spajaju temena i polovinu naspramnih strana
         % ref: http://www.mathopenref.com/trianglecentroid.html
         function [rez] = izracunajCentroidTrouglaAB(A,B,C)
             rez = AngioIvusMath.getTackuIzmeldjuDveLinije(A, (B+C)/2, C, (A+B)/2);
             plotLine([A;B;C;A]);
             plotLine([rez; (B+C)/2]);
             plotLine([rez; (A+C)/2]);
             plotLine([rez; (A+B)/2]);
             axis equal;
         end
         %% Izocentar se nalazi na preseku linija koje prolaze kroz temena i polovinu ugla na tom temenu
         % ref: http://www.mathopenref.com/triangleincenter.html
         function [rez vAO vBO vCO] = izracunajIncentarTrouglaAB(A,B,C)
             vAB = AngioIvusMath.arsUnit(B-A);
             vAC = AngioIvusMath.arsUnit(C-A);
             pravacPolovinaUglaKodTemenaA = A + (vAB+vAC)/2;
             vCA = AngioIvusMath.arsUnit(A-C);
             vCB = AngioIvusMath.arsUnit(B-C);
             pravacPolovinaUglaKodTemenaC = C + (vCA+vCB)/2;
             rez = AngioIvusMath.getTackuIzmeldjuDveLinije(A, pravacPolovinaUglaKodTemenaA, C,pravacPolovinaUglaKodTemenaC);
                 plotLine([A;B;C;A]);
                 plotLine([rez; A]);
                 plotLine([rez; B]);
                 plotLine([rez; C]);
                 axis equal;
            vAO = AngioIvusMath.arsUnit(rez-A);
            vBO = AngioIvusMath.arsUnit(rez-B);%
            vCO = AngioIvusMath.arsUnit(rez-C);
         end
         %% Pravi tacke elipse (ili dela elipse) u zavisnosti od ulaza
         %INPUTS
            % O      - centar elipse u prostoru
            % vX, vY - vektor pravca x/y ose
            % a, b   - parametri elispe 
         function rezTacke = napraviElipsu(O, vX, vY, a, b, pocetniUgao,dPi, krajnjiUgao)
            vX = AngioIvusMath.arsUnit(vX); vY = AngioIvusMath.arsUnit(vY);
            alfa=pocetniUgao:(krajnjiUgao-pocetniUgao)/dPi:krajnjiUgao; alfa(end)=krajnjiUgao;             
            x = cos(alfa)*a; %polovinaKrugaX(end)=polovinaKrugaX(1);
            y = sin(alfa)*b; %polovinaKrugaY(end)=polovinaKrugaY(1);
            rezTacke = AngioIvusMath.arsPuta(vX, x(:)) + AngioIvusMath.arsPuta(vY, y(:));
            rezTacke = AngioIvusMath.arsPlus(rezTacke, O);
%             plotLine(rezTacke);
         end
         %% Na plus X a=a1, na -X a=a2
         %INPUTS
            %vX1 == +x, vX2 == -x
         function rezTacke = napraviElipsu_BifurkacijaPatch(O, vX1, vX2, vY, a1, a2, b, pocetniUgao,dPi, krajnjiUgao)
            vX1 = AngioIvusMath.arsUnit(vX1); vX2 = AngioIvusMath.arsUnit(vX2); vY = AngioIvusMath.arsUnit(vY);
            alfa=pocetniUgao:(krajnjiUgao-pocetniUgao)/dPi:krajnjiUgao; alfa(end)=krajnjiUgao;  
            idA1 = find(cos(alfa)>=0);% +X
            idA2 = find(cos(alfa)< 0);% -X 
            x(idA1) = cos(alfa(idA1))*a1;  y(idA1) = sin(alfa(idA1))*b; %polovinaKrugaX(end)=polovinaKrugaX(1);
            x(idA2) = cos(alfa(idA2))*a2;  y(idA2) = sin(alfa(idA2))*b; %polovinaKrugaY(end)=polovinaKrugaY(1);
            x = x(:); y = y(:);
            rezTacke(idA1,:) = AngioIvusMath.arsPuta(vX1, x(idA1)) + AngioIvusMath.arsPuta(vY, y(idA1));%plus  X
            rezTacke(idA2,:) = AngioIvusMath.arsPuta(vX2, abs(x(idA2))) + AngioIvusMath.arsPuta(vY, y(idA2));%minus X
            rezTacke = AngioIvusMath.arsPlus(rezTacke, O);
%             plotLine(rezTacke);
         end  
         %% Funkcija proverava simetricnost matrice
         %INPUTS
            %matrica-matrica koja se proverava
         %OUTPUTS
            %elementi iznad dijagonale minus ispod
         function rez=proveriSimetricnostMatrice(matrica)
%            pomMatrica = ones(size(matrica));  
           elmentiIznadDijagonale = triu(matrica);
           elmentiIspodDijagonale = tril(matrica');
           rez = elmentiIznadDijagonale-elmentiIspodDijagonale;
           rez = rez.*rez;
           rez = mean(abs(rez(:)));
         end
         %%
         %  A.________.D
         %   |        |
         %   |        |
         %   |        |
         % B.|________.C
         function rez=proveriSimetricnostPodMatrice(matrica, k, y0)
             y0 = round(y0)
             if nargin==1
                 AngioIvusMath.proveriSimetricnostMatrice(matrica);
             else
                 [xSize,ySize]=size(matrica);
                 A  = [0       0              0]; 
                 B  = [0       ySize          0];
                 C  = [xSize   0              0]; 
                 D  = [xSize   ySize          0];
                 Po = [-xSize  (k*-xSize+y0)  0]; 
                 P1 = [2*xSize (k*2*xSize+y0) 0];
                 plotLine([A;B;C;D]);
                 plotLine([Po;P1]);
                 [rez] = round(AngioIvusMath.getTackePresekaLinijeP1P2iCetvorouglaABCD(Po,P1,A,B,C,D));
                 rez   = sortrows(rez);
                 %koriguje ako je van granica matrice
                 rez(rez<1)=1; rez=rez(:,1:2);  
                 if rez(2,1)>xSize
                    rez(2,1)=xSize;
                 end
                 if rez(2,2)>ySize
                    rez(2,2)=ySize;
                 end
                 %obe dimenzije moraju biti iste
                 dimPodMatrice = min(rez(2,:)-rez(1,:));
                 %nadji asimetricnost matrice
%                  if (y0>0)%dijagonala podmatrice je iznad dijagonale glavne matrice
                 podMatrica = matrica(rez(1):rez(1)+dimPodMatrice,rez(3):rez(3)+dimPodMatrice);
                 rez        = AngioIvusMath.proveriSimetricnostMatrice(podMatrica);
%                  else
%                  end
                 %--                 
             end
         end         
         %% Evaluacija furijeovog reda (parametri i nterval su poznati)
         %INPUTS
            %x - x-osa, interval na kome se izvrasava periodicna fja
            %w -  w, Ao, A, B - koeficijenti furijevog reda, videti obrazac
         %OUTPUT
            %Fx - y-osa
         %ref - http://www.mathworks.com/help/curvefit/fourier.html
             function [Fx x_Fx]=FurijeovRedEval(x, w, Ao, A, B)
                 if nargin==3%Za potrebe optimizacije
                     N=Ao;%argin(3);
                     X=x;%argin(1);                     
                     x=w;%argin(2);
                     w=X(1);
                     Ao=X(2);
                     A =X(3:2+N);
                     B =X(3+N:end);
                 end
                 iwx = 1:numel(A);
                 for i=1:numel(x)
                     x_i=A;
                     x_i(1:end)=x(i);
                     Fx(i)= Ao + sum(A.*cos(x_i.*iwx*w)) + sum(B.*sin(x_i.*iwx*w));
                 end
                 x_Fx=[x(:),Fx(:)];
%                  plotLine(x_Fx);
             end
          %% Fitovanje furijeovog reda na podatke
          %ref-http://www.mathworks.com/help/curvefit/fourier.html
          %INOUTS
            %X,Y - podaci na koje se fituje periodicna kriva
            %N   - br. harmonika furijeovog reda
          %OUTPUTS
            %w, Ao, A, B - koeficijenti furijevog reda, videti obrazac
          function fja=FitujFurijeovRed(X,Y,N)
              if nargin<3
                  N=2;
              end
              fja = fit(X,Y,['fourier' num2str(N)]);              
          end
          %% Funckija za uzimanje koordinate-geometrije sa slike
          %INPUTS
            % Img              - Putanja ka slici ili matrica slike
            % KooPocetakPixeli - Piksel gde treba da se nalazi koo. pocetak
            % xScaleFactor     - Veli?ina pixela u mm u horizontalnom pravcu
            % yScaleFactor     - Veli?ina pixela u mm u vertikalnom   pravcu     
            % brKorakaInterpolacije - ako treba da se interoplira splajnom!
          %OUTPUTS
            % rez - koordinate
          function [rez rez1]= pokupiKoordinateSaSlike(Img, KooPocetakPixeli, xScaleFactor, yScaleFactor, brKorakaInterpolacije)
              %Parametri 
              if isstr(Img)
                  Img = imread(Img);
                  Img = Img(:,:,1);
              end
              if ~exist('KooPocetakPixeli')
                  KooPocetakPixeli = [1 1];
                  Xmm              = 1    ;
                  Ymm              = 1    ;
              elseif ~exist('Xmm') || ~exist('Xmm')                  
                  Xmm              = 1    ;
                  Ymm              = 1    ;
              end
              % Pokupi koordinate sa slike
                h = figure  ;
                hold on     ;
                x = 0       ;
                while(numel(x)<3)%mora da unese najmanje 3 koordinate
                 imshow(Img);
                 hold     on;          
                 axis    off;
                 text(16,36,'Left click to get points, right click to get end point','FontSize',[12],'Color', 'r');
                 [region x y] = roipoly ;      
                 close(h);
                 rez = [x(1:end-1) y(1:end-1)];
                end%while
                rez(:,1) =  (rez(:,1)-KooPocetakPixeli(1))/xScaleFactor;
                rez(:,2) = -(rez(:,2)-KooPocetakPixeli(2))/yScaleFactor;
                rez(:,3) = 0;
          rez1='Ne treba interpolacija';
          if exist('brKorakaInterpolacije')%treba interpolacija
              rez1 = rez;
              splajn = arsCubicBSplineOpenNonuniform;
              splajn = splajn.set(splajn, rez, 3);
              rez    = splajn.interpolateFromAtoB(splajn,0,1/brKorakaInterpolacije,1);
          end
        end
        %% Calculates the Hausdorff Distance between A and B 
        %  H (A, B) = max { h (A, B), h (B, A) }
        % REF http://cgm.cs.mcgill.ca/~godfried/teaching/cg-projects/98/normand/main.html
        %INPUTS 
            % A, B - niz tacaka
        %OUTPUTS
            % hd - HausdorffDist celog point seta
            % D  - HD rastojanja po cvorivma niz a
        function [hd D idCvoraB] = HausdorffDist(A, B)
            for i = 1 : numel(A(:,1))
                [hd(i) idCvoraB(i)] = min(AngioIvusMath.arsNorm(AngioIvusMath.arsMinus(A(i,:), B)));
            end
            D = min(hd);
        end
	end%methods
end