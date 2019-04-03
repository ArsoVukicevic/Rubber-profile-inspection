%%%
%INPUTS
    %konturaProfila     - Detektovana kontura profila
    %img                - Putanja do img
    %scaleFaktor        - Scale faktor za prebacivanje distanci iz pixela mm
    %debugMode          - iscrtavanje medjuresultata true/false
%OUTPUTS
    %skeletonImg           - centralne linije
    %geodesicImg           - matrica rastojanja pixela od centralnih linija
    %koordinateBifurkacija - koordinate kljucnih tacki
    %koordianateCentralnihLinija
function [geodesicImg skeletonImg koordinateBifurkacija koordianateCentralnihLinija] = skripta4_OdradiEuclidianDistanceNadjiSkeletoneKontura_b(konturaProfila, img, scaleFaktor, debugMode)
if isstr(img)
    img = imread(img);
    img = img(:,:,1);
end
%#1 Napravi poly mask od konture
    [x y]=size(img);
    imgPolyMask = poly2mask(konturaProfila(:,1), konturaProfila(:,2), x, y);
    imgPolyMask = int16(~imgPolyMask);
    if debugMode
        figure;
        imshow(mat2gray(imgPolyMask));
    end
    
%#2 Nadji Centralne Linije
    skeletonImg                      = bwmorph(~imgPolyMask,'skel',inf);
    koordianateCentralnihLinija      = arsIMG.pxelIDS2Koordinate(skeletonImg, find(skeletonImg==1));
    koordianateCentralnihLinija(:,3) = 0;
    koordianateCentralnihLinija      = AngioIvusMath.arsPlus(koordianateCentralnihLinija, [1 1 0]);

%#3 Nadji Bifurkacije - one ce da sluze za alignement kontura
    C                          = conv2(double(skeletonImg), [1 1 1; 1 1 1; 1 1 1])  ;
    idBifurkacija              = find (C>3)                                 ;
    koordinateBifurkacija      = arsIMG.pxelIDS2Koordinate(C, idBifurkacija);
    koordinateBifurkacija(:,3) = 0;
   
 %#4 Nadji geodesic distance         
    skeletonImg = bwmorph(logical(~imgPolyMask),'thin',33); 
    [geodesicImg] = bwdistgeodesic(logical(~imgPolyMask), logical(skeletonImg)); 
    geodesicImg = scaleFaktor * geodesicImg;
    
    
% debug    
     if debugMode
        figure;
        imshow(mat2gray(C)); hold on;
        scatter(koordianateCentralnihLinija(:,2), koordianateCentralnihLinija(:,1));
        scatter(koordinateBifurkacija(:,2), koordinateBifurkacija(:,1));
        plotLine(konturaProfila);
        
        figure;
        contourf(geodesicImg*2);
        colorbar
        plotLine(konturaProfila);
        axis equal;
    end
end

%% 
function h_out = colorscale(varargin)
%COLORSCALE Color scale bar using part of figure's colormap
%
%   COLORSCALE creates a color scale for use with an integer-valued
%   image which has been encoded and quantized via a monotonic
%   transformation from a real-valued array.  The original array has
%   physically meaningful units, but the image does not.  The
%   transformed image has been displayed with a direct mapping to the
%   figure's colormap. (That is, the IMAGE object has a CDATAMAPPING of
%   'direct'.)  The image and/or color scale may not use the entire
%   colormap.  The colormap may contain special values for null data
%   areas, and these special values may be excluded from the scale bar.
%   Or, the colormap may be deliberately partitioned and shared with a
%   different image in the same figure, and that image may have its own,
%   entirely distinct, color scale. In these cases, unlike the
%   standard COLORBAR function, COLORSCALE can limit the displayed
%   colors to the appropriate interval within the figure's colormap.
%
%   H = COLORSCALE(CMAPLIM,DATALIM,DATAINC,ORIENTATION,'Position',RECT)
%   adds a color scale to the current figure.  The color scale uses the
%   part of the figure's colormap specified by limits in the 1-by-2
%   vector CMAPLIM.  The ticks on the vertical axis are in physical data
%   units which may be scaled and shifted with respect to the colormap
%   indices and are controlled by DATALIM and DATAINC.  DATALIM is a
%   1-by-2 vector such that DATALIM(1) and DATALIM(2) indicate the
%   physical values corresponding to the colors in rows CMAPLIM(1) and
%   CMAPLIM(2), respectively, of the colormap.  DATALIM and DATAINC are
%   in physical data units.  DATALIM(2) must exceed DATALIM(1), but
%   CMAPLIM(2) may be less than CMAPLIM(1). ORIENTATION should be either
%   'horiz' for horizontal or 'vert' for vertical. RECT = [left, bottom,
%   width, height] specifies the location and size of the axis box
%   containing the color scale (i.e., the 'Position' property of the
%   axes containing the color scale).  Unlike COLORBAR, COLORSCALE does
%   not automatically decide where to place the axes.  COLORSCALE
%   optionally returns the axes handle.
%
%   H = COLORSCALE(...,PROP1,VAL1,PROP2,VAL2,...) sets additional
%   properties of the axes containing the color scale.
%
%   Example
%   -------
%   % Create a figure with a jet colormap
%   figure
%   colormap jet
% 
%   % Create a standard colorbar
%   set(axes,'Visible','off')
%   bar = colorbar;
%   set(get(bar,'Title'),'String','Standard Colorbar')
%   pos = get(bar,'Position');
% 
%   % Create a terrain height scale using part of the colormap
%   colorscale([31 54], [0 9000], 1000, 'vert',...
%              'Position',pos - [0.7 0 0 0])
%   title('Terrain Height')
%   ylabel('elevation above sea level, meters')
% 
%   % Create a bathymetry scale that uses a different part of the
%   % colormap and runs in reverse (i.e., CMAPLIM(1) > CMAPLIM(2)).
%   colorscale([24 1], [0 12000], 1000, 'vert',...
%              'Position',pos - [0.35 0 0 0],'YDir','reverse')
%   title('Bathymetry')
%   ylabel('depth below sea level, meters')
%
%   See also COLORBAR, COLORMAP, COLORMAPEDITOR, AXES.
% Written by: Rob Comer
% Copyright 2002-2011 The MathWorks, Inc. 
error(nargchk(6,inf,nargin,'struct'))
cmaplim = varargin{1};
datalim = varargin{2};
datainc = varargin{3};
orientation = varargin{4};
% Determine tick labels
tickvalues = datainc * (ceil(datalim(1)/datainc) : floor(datalim(2)/datainc));
ticklabel = cell(size(tickvalues));
for i = 1 : length(tickvalues)
  ticklabel{i} = num2str(tickvalues(i));
end
% Determine tick positions and image array
if cmaplim(1) < cmaplim(2)
    tick = cmaplim(1) + (tickvalues - datalim(1)) * diff(cmaplim) / diff(datalim);
    I = cmaplim(1):cmaplim(2);
    ilim = cmaplim;
else
    tick = cmaplim(2) - (tickvalues - datalim(1)) * diff(cmaplim) / diff(datalim);
    I = cmaplim(1):-1:cmaplim(2);
    ilim = [cmaplim(2) cmaplim(1)];
end
% Create the color bar axes and image object
h = axes('Visible','off');
if strncmpi(orientation,'vertical',length(orientation))
    % Vertical scale
    image(I','CDataMapping','direct','YData',ilim);
    set(h,'YTick',tick,...
          'XTick',[],...
          'YTickLabel',ticklabel);
else
    % Horizontal scale
    image(I,'CDataMapping','direct','XData',ilim);
    set(h,'XTick',tick,...
          'YTick',[],...
          'XTickLabel',ticklabel);
end
% Set axes additional properties
set(h,'YDir','normal'); % Override the default YDir set by IMAGE
set(h,varargin{5:end}); % Separate in case we have {...,'YDir','reverse',...}
set(h,'Visible','On');  % Always make it visible
% Set return value, if any
if nargout == 1
    h_out = h;
end
end