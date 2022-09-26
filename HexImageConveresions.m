% hex-coordinate functions based on redblobgames
% https://www.redblobgames.com/grids/hexagons
% Artistic style inspired by two artists residing in Rotterdam,
% Netherlands: Michael Bom and Boris Pas. 
%
% Written by Aaron B. Wong (2022)
%
% Written in Matlab R2019a.

%% 
clear; close all
ImgMat = imread('Kek-Look-Tong-square.JPG'); hexWidthPx = 25;

%% Coordinate of image pixels
xDim = size(ImgMat,2); yDim = size(ImgMat,1);
nPixel = xDim * yDim;
xx = 1:xDim; yy = 1:yDim;
% ImgPx = struct([]);
ImgPx.pxy = combvec(xx,yy); % xy coordinates of all pixel in image. Upper left pixel: (1,1)
ImgPx.qrs = p2h(ImgPx.pxy); % pixel positions in cubic coordinates
% ImgPx = struct( 'x',num2cell(ImgPx.pxy(1,:)),...
%                 'y',num2cell(ImgPx.pxy(2,:)),...
%                 'q',num2cell(qrs(1,:)),...
%                 'r',num2cell(qrs(2,:)),...
%                 's',num2cell(qrs(3,:)) ...
%                 );
clear('xx','yy');
%% calculate downsampled grid that spans the image

% grid.qrs = unique(cube_round(ImgPx.qrs)','row')';
% nHexel = size(grid.qrs,2);
% grid.pxy = h2p(grid.qrs);
disp(['Number of pixels: ',num2str(nPixel), ' (',num2str(xDim),' x ',num2str(yDim),')'])
% disp(['Number of hexels: ',num2str(nHexel)])
% expand grid
samp.qrs = unique(cube_round(ImgPx.qrs ./ hexWidthPx)','row')';
samp.pxy = h2p(samp.qrs)*hexWidthPx;
% exclude out of bound hexels
sel = samp.pxy(1,:) <= xDim & samp.pxy(2,:) <= yDim & samp.pxy(1,:) > 0 & samp.pxy(2,:) > 0;
samp.qrs = samp.qrs(:,sel);samp.pxy = samp.pxy(:,sel);
% get column-row (doublewidth) representation
samp.cr = ax2dw(samp.qrs);
nSamp = size(samp.qrs,2);
disp(['Number of downsampled hexels: ',num2str(nSamp)])
clear('grid','sel');
%% calculate output size
% specify canvas size
canvasWidth = 50; % cm
canvasHeight = 50; % cm
minCol = min(samp.cr(1,:)); maxCol = max(samp.cr(1,:));
minRow = min(samp.cr(2,:)); maxRow = max(samp.cr(2,:));
nCols = 0.5 * (maxCol - minCol) + 1;
nRows = (maxRow - minRow) + 1;

% distance between centers of adjacent hexels
margin = 0.25; %
a = sqrt(3) / 2;
HexDist1 = canvasWidth / (nCols+2*margin); % 2 * margin on both sides
HexDist2 =  canvasHeight / (nRows+2*margin) /  a; %distance between centers of adjacent hexels
HexDist = min(HexDist1,HexDist2);

canvasXOffset = 0.5 * canvasWidth - 0.5*(nCols-1)*HexDist;
canvasYOffset = 0.5 * canvasHeight - 0.5*(nRows-1)*a*HexDist;
minPntXY = HexDist*h2p(dw2ax([minCol;minRow]));
canvasOffset = [canvasXOffset;canvasYOffset] - minPntXY;

samp.xy = HexDist*h2p(samp.qrs)+canvasOffset;

%% sampling pixel intensity
% & generate pixelation image
nCh = size(ImgMat,3);
CData_samp = nan(size(samp.pxy,2),nCh);
CData_mean = nan(size(samp.pxy,2),nCh);
HexelateImg_samp = nan(size(ImgMat));
HexelateImg_mean = nan(size(ImgMat));

% nearest central pixel
tic
for ii = 1:size(samp.pxy,2)
    CData_samp(ii,:) = squeeze(ImgMat(round(samp.pxy(2,ii)),round(samp.pxy(1,ii)),:));
end
disp('Sampling done.')
toc

% mean 
tic

mapHex = cube_round(ImgPx.qrs ./ hexWidthPx);
for ii = 1:size(samp.pxy,2)
    sel = all(samp.qrs(:,ii) == mapHex,1);
%     subplot(1,3,2)
%     scatter(ImgPx.pxy(1,sel),ImgPx.pxy(2,sel))
    nn = sum(sel);
    for ch = 1:nCh
        linearInd = sub2ind([yDim,xDim,nCh], ImgPx.pxy(2,sel), ImgPx.pxy(1,sel),ch*ones(1,nn));
        CData_mean(ii,ch) = mean(ImgMat(linearInd),'all');
        HexelateImg_samp(linearInd) =CData_samp(ii,ch);
        HexelateImg_mean(linearInd) =CData_mean(ii,ch);
    end
end
disp('Averaging done.')
toc



%% output to JSON

ch = 2;
fileID = fopen('Hexelated.json','w');
JSON = jsonencode(table(samp.pxy(1,:)',samp.pxy(2,:)',...
                samp.qrs(1,:)',samp.qrs(2,:)',samp.qrs(3,:)',...
                samp.cr(1,:)',samp.cr(2,:)',...
                CData_mean(:,ch),...
                'VariableNames',{'px','py','q','r','s','col','row','c'}));
fwrite(fileID,JSON);
fclose(fileID);
%% show image on Canvas dimension
fig = figure;
screenScale = 7; % "points" / cm
figSize = screenScale*[canvasWidth,canvasHeight]; figMargin = [50,50];
fig.Units = 'points';
fig.Position = [300,50,(figSize+2*figMargin)];
nSteps = 8;
ax = gca;
ax.Units = 'points';
ax.Position = [figMargin,figSize];
ch = 2;
maxMrkrSize = pi * (0.5*(HexDist*screenScale))^2; % somehow pi*r^2 is not correct
maxMrkrSize = (HexDist*screenScale)^2; % using the area of a square works 
sz = maxMrkrSize*szLU(-CData_mean(:,ch),nSteps);
scatter(samp.xy(1,:),samp.xy(2,:),sz,'filled','Marker','o');
set(gca,'YDir','reverse');
ylim([0,canvasHeight]);xlim([0,canvasWidth]);
ylabel('cm')
xlabel('cm')
%% Show images (R,G,B, RGB)
fig = figure;
fig.Position = [300,50,750,750];
Colors = {'r','g','b'};%{'c','m','y'};
nSteps = 8;
for ch = 1:nCh
    subplot(2,2,ch)
        sz = 10*szLU(-CData_mean(:,ch),nSteps);%round((280-CData_mean(:,ch))./64)+eps;

%     scatter(samp.pxy(1,:),samp.pxy(2,:),sz,Colors{ch},'filled');
    scatter(samp.pxy(1,:),samp.pxy(2,:),sz,'filled');
    set(gca,'YDir','reverse');
%     axis square
        ylim([0,yDim+1]);xlim([0,xDim+1]);
    title(Colors{ch});
end
subplot(2,2,4)
scatter(samp.pxy(1,:),samp.pxy(2,:),9,CData_mean./255,'filled');
set(gca,'YDir','reverse');
    ylim([0,yDim+1]);xlim([0,xDim+1]);

% axis square

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
% saveas(fig,'sampOutput.pdf')
%% demo

% calculate axis labels
sel = samp.cr(2,:) == min(samp.cr(2,:));
ColLbl_x = samp.pxy(1,sel);
ColLbls = samp.cr(1,sel);
sel = samp.cr(1,:) == min(samp.cr(1,:));
RowLbl_x = samp.pxy(2,sel);
RowLbls = samp.cr(2,sel);
[RowLbls,idx] = sort(RowLbls);
RowLbl_x = RowLbl_x(idx);

clear('sel','idx');


fig = figure;
fig.Position = [200,300,900,300];

nSteps = 8;

% original image
subplot(1,3,1)
    imagesc(ImgMat);
    axis square
    
% scaled
    ch = 3; % channel to use
    dotScale = 2;
    useMean = 0; % 1 = mean, 0 = decimated samples
    sz_mean = (nSteps-1)*szLU(-CData_mean(:,ch),nSteps);%round((280-CData_mean(:,ch))./32);
    sz_samp = (nSteps-1)*szLU(-CData_samp(:,ch),nSteps);%round((280-CData_samp(:,ch))./32);
    if(useMean); sz = sz_mean; else; sz = sz_samp; end
    [uSz,ia,ic] = unique(sz);
    nColors = length(uSz);
%     CScale = colorcube(nColors);
    CScale = 0.8.*jet(nColors);
%     CScale = lines(nColors);
    c = CScale(ic,:);
    subplot(1,3,2)
    scatter(samp.pxy(1,:),samp.pxy(2,:),dotScale*sz,c,'filled');
    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square

    subplot(1,3,3)
    scatter(samp.pxy(1,:),samp.pxy(2,:),15,c,'filled');
%         scatter(samp.pxy(1,:),samp.pxy(2,:),sz_mean,'filled');

    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square
colormap(CScale);
caxis([min(uSz)-.5,max(uSz)+.5])
colorbar('Ticks',uSz,'Position',[.92,.3,.02,.5]);

% ticks
xticks(ColLbl_x);xticklabels(ColLbls);
yticks(RowLbl_x);yticklabels(RowLbls);

% label extremes
locExtremes = nan(4,1);
[~,locExtremes(1)] = min([1,1]*samp.pxy,[],2);
[~,locExtremes(2)] = min([1,-1]*samp.pxy,[],2);
[~,locExtremes(3)] = max([1,1]*samp.pxy,[],2);
[~,locExtremes(4)] = max([1,-1]*samp.pxy,[],2);

% locExtremes = union(locMinCR,locMaxCR);
for ii = 1:length(locExtremes)
    idx = locExtremes(ii);
    text(samp.pxy(1,idx),samp.pxy(2,idx),['(',num2str(samp.cr(1,idx)),',',num2str(samp.cr(2,idx)),')'],...
        'FontSize',8);
end

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
% saveas(fig,'sampOutput.pdf')
saveas(fig,'sampOutput_mean.png')

%% Instruction
fig = figure;
fig.Position = [100,100,1200,500];
    dotScale = 12;
    subplot(1,2,1)
    scatter(samp.pxy(1,:),samp.pxy(2,:),dotScale*sz,'filled');
    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square
    subplot(1,2,2)
    scatter(samp.pxy(1,:),samp.pxy(2,:),dotScale*nSteps/2,c,'filled');

    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    
    
    axis square
colormap(CScale);
caxis([min(uSz)-.5,max(uSz)+.5])
colorbar('Ticks',uSz,'Position',[.92,.3,.02,.5]);

% label extremes
locExtremes = nan(4,1);
[~,locExtremes(1)] = min([1,1]*samp.pxy,[],2);
[~,locExtremes(2)] = min([1,-1]*samp.pxy,[],2);
[~,locExtremes(3)] = max([1,1]*samp.pxy,[],2);
[~,locExtremes(4)] = max([1,-1]*samp.pxy,[],2);

% locExtremes = union(locMinCR,locMaxCR);
for ii = 1:length(locExtremes)
    idx = locExtremes(ii);
    text(samp.pxy(1,idx),samp.pxy(2,idx),['(',num2str(samp.cr(1,idx)),',',num2str(samp.cr(2,idx)),')'],...
        'FontSize',8);
end

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
% saveas(fig,'instruction.pdf')
%% Hexelate Images
fig = figure;
fig.Position = [200,300,900,250];
    subplot(1,3,1)
    imagesc(ImgMat); axis square;
    subplot(1,3,2)
    image(HexelateImg_samp./255); axis square;
    subplot(1,3,3)
    image(HexelateImg_mean./255); axis square;

%%
% here, size of 1 is defined as distance between centers of adjacent
% hexagons
function [xy] = h2p(qrs)
    % pointy top hex
    % qrs: 3x1 vector containing [q;r;s]
%     A = [   sqrt(3),    0.5*sqrt(3);...
%             0,          1.5         ];
    A = [   1,    0.5;...
            0,    0.5*sqrt(3)];
    xy = A * qrs(1:2,:);
end

function qrs = p2h(xy)
     % pointy top hex
     % xy: 2x1 vector containing [px;py]
%     A = [   sqrt(3)/3,   -1/3;...
%             0,          2/3        ];
    A = [   1,   -1/sqrt(3);...
            0,   2/sqrt(3) ];
    qr = A * xy;
    qrs = [qr;-sum(qr,1)];
end

function dist = cube_distance(a, b)
    dist = sum(abs(a-b),1) / 2;
end

function cr = ax2oddr(hex)
    col = hex(1,:) + (hex(2,:) - bitand(hex(2,:),1)) / 2;
    row = hex(2,:);
    cr = [col;row];
end
function cr = ax2dw(hex) % axial_to_doublewidth
    col = 2 * hex(1,:) + hex(2,:);
    row = hex(2,:);
    cr = [col;row];
end
function hex = oddr2ax(cr)
    q = cr(1,:) - (cr(2,:) - bitand(cr(2,:),1)) / 2;
    r = cr(2,:);
    hex = [q;r;-q-r;];
end
function hex = dw2ax(cr) % axial_to_doublewidth
    q = (cr(1,:) - cr(2,:)) / 2;
    r = cr(2,:);
    hex = [q;r;-q-r;];
end

function qrs = cube_round(frac)
%     var q = round(frac.q)
%     var r = round(frac.r)
%     var s = round(frac.s)
    qrs = round(frac);

%     var q_diff = abs(q - frac.q)
%     var r_diff = abs(r - frac.r)
%     var s_diff = abs(s - frac.s)
    qrs_diff = abs(qrs - frac);
    
%     if q_diff > r_diff and q_diff > s_diff:
%         q = -r-s
%     else if r_diff > s_diff:
%         r = -q-s
%     else:
%         s = -q-r
    [~,idx] = max(qrs_diff,[],1,'linear');
    qrs(idx) = 0;
    qrs(idx) = -sum(qrs,1);

%     return Cube(q, r, s)
end

function sz = szLU(vals,nSteps)
    minV = min(vals,[],'all');
    maxV = max(vals,[],'all');
    
    sz = (round ( (nSteps - 1) .* (vals - minV) ./ (maxV-minV ) ) ./ (nSteps-1)) + eps;
end