% hex-coordinate functions based on redblobgames
% https://www.redblobgames.com/grids/hexagons
% Artistic style inspired by two artists residing in Rotterdam,
% Netherlands: Michael Bom and Boris Pas. 
%
% Written by Aaron B. Wong (2022)
%
% Written in Matlab R2019a.

%% Clear workspace and load the image
clear; close all
ImgMat = imread('sample\Leonardo-Mona-Lisa.jpg');hexWidthPx = 30;


%%
nRow_target = 40;
hexImgAvg_double = hexSamp(ImgMat,nRow_target);

% Create a linear array of pixel values and their x and y coordinates
pixelArray = hexDouble2Array(hexImgAvg_double);
%% display the image

subplot(2,3,1)
image(ImgMat)
subplot(2,3,2)
image(uint8(hexImgAvg_double))
nRows = size(hexImgAvg_double,1);
nCols = size(hexImgAvg_double,2);
for ch = 1:3
    ax = subplot(2,3,3+ch);
    nSteps = 7;
    sz = szLU(-double(pixelArray(:,ch+2)),nSteps);
    maxMrkrSize = 20;
    Color = 'k';
    scatter(pixelArray(:,1),pixelArray(:,2),maxMrkrSize.*sz,Color,'filled','Marker','o');
    ax.YDir = 'reverse';
    ax.YLim = [0.5,nRows+0.5];
    ax.XLim = [0.5,nCols+0.5];
end
%%
% %% calculate downsampled grid that spans the image
% samp.cr = unique(dw_round(p2dw(ImgPx.pxy)./ hexWidthPx)','row','stable')';
% samp.pxy = dw2p(samp.cr)*hexWidthPx;
% % exclude out of bound hexels
% sel = samp.pxy(1,:) <= xDim & samp.pxy(2,:) <= yDim & samp.pxy(1,:) > 0 & samp.pxy(2,:) > 0;
% samp.pxy = samp.pxy(:,sel);
% samp.cr = samp.cr(:,sel);samp.qrs = dw2cube(samp.cr);
% nSamp = size(samp.pxy,2);

% disp(['Number of pixels: ',num2str(nPixel), ' (',num2str(xDim),' x ',num2str(yDim),')'])
% disp(['Number of downsampled hexels: ',num2str(nSamp)])
% clear('sel');

% %% calculate output size
% % specify canvas size
% canvasWidth = 50; % cm
% canvasHeight = 50; % cm
% minCol = min(samp.cr(1,:)); maxCol = max(samp.cr(1,:));
% minRow = min(samp.cr(2,:)); maxRow = max(samp.cr(2,:));
% nCols = 0.5 * (maxCol - minCol) + 1;
% nRows = (maxRow - minRow) + 1;

% disp(['Number of columns: ',num2str(nCols),...
%       '   Number of rows: ',num2str(nRows)]);

% % distance between centers of adjacent hexels
% margin = 0.25; % fraction of 1 Hex
% a = sqrt(3) / 2; % vertical distance between adjacent rows
% HexDist1 = canvasWidth / (nCols+2*margin); % 2 * margin on both sides
% HexDist2 =  canvasHeight / (nRows+2*margin) /  a; % 2 * margin on both sides
% HexDist = min(HexDist1,HexDist2);

% canvasXOffset = 0.5 * canvasWidth - 0.5*(nCols-1)*HexDist;
% canvasYOffset = 0.5 * canvasHeight - 0.5*(nRows-1)*a*HexDist;
% minPntXY = HexDist*h2p(dw2cube([minCol;minRow]));
% canvasOffset = [canvasXOffset;canvasYOffset] - minPntXY;

% samp.xy = HexDist*dw2p(samp.cr)+canvasOffset;

% %% sampling pixel intensity
% % & generate pixelation image
% nCh = size(ImgMat,3);
% samp.c_samp = nan(size(samp.pxy,2),nCh);
% samp.c_mean = nan(size(samp.pxy,2),nCh);
% HexelateImg_samp = nan(size(ImgMat));
% HexelateImg_mean = nan(size(ImgMat));

% % nearest central pixel
% tic
% for ii = 1:size(samp.pxy,2)
%     samp.c_samp(ii,:) = squeeze(ImgMat(round(samp.pxy(2,ii)),round(samp.pxy(1,ii)),:));
% end
% disp('Sampling done.')
% toc

% % mean 
% tic

% % mapHex = cube_round(ImgPx.qrs ./ hexWidthPx);
% mapHex = dw_round(p2dw(ImgPx.pxy) ./ hexWidthPx);

% for ii = 1:size(samp.pxy,2)
%     sel = all(samp.cr(:,ii) == mapHex,1);
%     nn = sum(sel);
%     for ch = 1:nCh
%         linearInd = sub2ind([yDim,xDim,nCh], ImgPx.pxy(2,sel), ImgPx.pxy(1,sel),ch*ones(1,nn));
%         samp.c_mean(ii,ch) = mean(ImgMat(linearInd),'all');
% %         HexelateImg_samp(linearInd) =samp.c_samp(ii,ch);
% %         HexelateImg_mean(linearInd) =samp.c_mean(ii,ch);
%     end
% end
% disp('Averaging done.')
% toc

clear('ii','nn','sel')

%% outputs

ch2Use = 2;
ch = max(ch2Use,nCh);
outTable = table(samp.pxy(1,:)',samp.pxy(2,:)',...
                samp.qrs(1,:)',samp.qrs(2,:)',samp.qrs(3,:)',...
                samp.cr(1,:)',samp.cr(2,:)',...
                samp.c_mean(:,ch),samp.c_samp(:,ch),...
                'VariableNames',{'px','py','q','r','s','col','row','c_mean','c_samp'});

% output to JSON

fileID = fopen([ImgName,'_hex.json'],'w');
JSON = jsonencode(outTable);
fwrite(fileID,JSON);
fclose(fileID);

% output as CSV
writetable(outTable,[ImgName,'_hex.csv']);

%% Show Hexelated Images
fig = figure;
fig.Position = [200,300,900,250];
    subplot(1,3,1)
    imagesc(ImgMat); axis square; title('original')
    subplot(1,3,2)
    image(HexelateImg_samp./255); axis square; title('sampled (decimate)')
    subplot(1,3,3)
    image(HexelateImg_mean./255); axis square; title('mean')
%% Show RGB Images
if nCh > 1
fig = figure;
fig.Position = [200,300,900,250];
    for ii = 1:nCh
        subplot(1,nCh,ii)
        imagesc(ImgMat(:,:,ii)); axis square; title(num2str(ii,'Channel %d'))
        colormap("gray")
    end
end
%% show image on Canvas dimension
% parameters
screenScale = 7; % points / cm; for display purpose
nSteps = 256;
ch2Use = 2;
ch = max(ch2Use,nCh); % which channel to use (1=R,2=G,3=B)
Color = [105,81,146]./255;
Color = [96,155,213]./255; % ~koningsblauw blue
% Color = [102,47,31]./255; % ~Transparent oxide brown


% --- setup figure ---
    fig = figure;
    figSize = screenScale*[canvasWidth,canvasHeight]; figMargin = [50,50];
    fig.Units = 'points';
    fig.Position = [300,50,(figSize+2*figMargin)];
    
    ax = gca;
    ax.Units = 'points';
    ax.Position = [figMargin,figSize];

% --- calculate marker size ---
    % maxMrkrSize = pi * (0.5*(HexDist*screenScale))^2; % somehow pi*r^2 is not correct
    maxMrkrSize = (HexDist*screenScale)^2; % using the area of a square works 
    sz = szLU(max(-samp.c_mean(:,ch),-200),nSteps); % reversed scale so darker hexels have larger dots

% --- plot dots ---
    scatter(samp.xy(1,:),samp.xy(2,:),maxMrkrSize.*sz,Color,'filled','Marker','o');

    set(gca,'YDir','reverse');
    ylim([0,canvasHeight]);xlim([0,canvasWidth]);
    ylabel('cm');xlabel('cm')
    

%     scatter(samp.pxy(1,:),samp.pxy(2,:),dotScale*nSteps/2,c,'filled');
%% Instruction for Dot artwork
screenScale = 4; % points / cm; for display purpose
nSteps = 8;
ch = max(ch2Use,nCh); % which channel to use (1=R,2=G,3=B)
Color = [0,0,0]; % color of dots
useMean = 1;





fig = figure;
    axSize = screenScale*[canvasWidth,canvasHeight];
    figSize = [2,1].*axSize; figMargin = [50,50];
    fig.Units = 'points';
    fig.Position = [10,10,(figSize+[2+1,2].*figMargin)];
    
% output picture
    ax = subplot(1,2,1);
    ax.Units = 'points';
    ax.Position = [figMargin,axSize];
    maxMrkrSize = (HexDist*screenScale)^2; % using the area of a square works 
    if(useMean)
        sz = szLU(-samp.c_mean(:,ch),nSteps);
    else
        sz = szLU(-samp.c_samp(:,ch),nSteps);
    end
    % --- create color lookup
        [uSz,ia,ic] = unique((nSteps-1).*sz);
        nColors = length(uSz);
        CScale = 0.8.*jet(nColors);
        c = CScale(ic,:);
    
    scatter(samp.xy(1,:),samp.xy(2,:),maxMrkrSize*sz,'filled','Marker','o','MarkerFaceColor',Color);
    set(gca,'YDir','reverse');
    ylim([0,canvasHeight]);xlim([0,canvasWidth]);
    ylabel('cm');xlabel('cm')

    ax = subplot(1,2,2);
    ax.Units = 'points';
    ax.Position = [[2,1].*figMargin+[1,0].*axSize,axSize];
    scatter(samp.xy(1,:),samp.xy(2,:),maxMrkrSize,c,'filled');

    set(gca,'YDir','reverse');
    ylim([0,canvasHeight]);xlim([0,canvasWidth]);
    ylabel('cm');xlabel('cm')
    
    colormap(CScale);
    caxis([min(uSz)-.5,max(uSz)+.5])
    colorbar('Ticks',uSz,'Position',[.92,.3,.02,.5]);

% label extremes
    locExtremes = nan(4,1);
    [~,locExtremes(1)] = min([1,1]*samp.pxy,[],2);
    [~,locExtremes(2)] = min([1,-1]*samp.pxy,[],2);
    [~,locExtremes(3)] = max([1,1]*samp.pxy,[],2);
    [~,locExtremes(4)] = max([1,-1]*samp.pxy,[],2);
    HorAlign = {'left','left','left','left'};
    VerAlign = {'bottom','top','top','bottom'};

    for ii = 1:length(locExtremes)
        idx = locExtremes(ii);
        text(samp.xy(1,idx),samp.xy(2,idx),['(',num2str(samp.cr(1,idx)),',',num2str(samp.cr(2,idx)),')'],...
            'FontSize',8,'HorizontalAlignment',HorAlign{ii},'VerticalAlignment',VerAlign{ii});
    end
% setup and save output file
fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./72; %72 points per inch
saveas(fig,[ImgName,'_instruction.pdf'])

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

function cr = p2dw(xy)
	% doublewidth
	% xy: 2xn matrix containing [px;py]
    A = [   2,   0;...
            0,   2/sqrt(3) ];
    cr = A * xy;
end

function xy = dw2p(cr)
     % doublewidth
     % cr: 2xn matrix containing [cols;rows]
    A = [   0.5,   0;...
            0,   0.5*sqrt(3) ];
    xy = A * cr;
end

function dist = cube_distance(a, b)
    dist = sum(abs(a-b),1) / 2;
end

function cr = ax2oddr(hex)
    col = hex(1,:) + (hex(2,:) - bitand(hex(2,:),1)) / 2;
    row = hex(2,:);
    cr = [col;row];
end
function cr = cube2dw(hex) % axial_to_doublewidth
    col = 2 * hex(1,:) + hex(2,:);
    row = hex(2,:);
    cr = [col;row];
end
function hex = oddr2cube(cr)
    q = cr(1,:) - (cr(2,:) - bitand(cr(2,:),1)) / 2;
    r = cr(2,:);
    hex = [q;r;-q-r;];
end
function hex = dw2cube(cr) % axial_to_doublewidth
    q = (cr(1,:) - cr(2,:)) / 2;
    r = cr(2,:);
    hex = [q;r;-q-r;];
end

function cr = dw_round(frac)
    cr = cube2dw(cube_round(dw2cube(frac)));
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

function hexImgAvg_double = hexSamp(ImgMat,nRow_target)
    ImgMat = double(ImgMat);
    %% Dimension of Image
    xDim = size(ImgMat,2); yDim = size(ImgMat,1);
    
    % Pointy-top orientation
    hexVertPitch_px = yDim/nRow_target; % an estimate
    hexSize_pix = 2/3 * hexVertPitch_px;
    hexHorizPitch_px = (sqrt(3) * hexSize_pix); % even number
    
    nRows = nRow_target;%floor(yDim / hexVertPitch_px);
    nCols = round(2 * xDim / hexHorizPitch_px);
    
    hexImgAvg_double = imresize(ImgMat,[nRows,nCols],"box");
%     hexImgSamp_double = imresize(ImgMat,[nRows,nCols],"nearest");
    
    %% averaging every two pixels (vectorized)

    odd_nCols = floor(nCols/2);
    odd_2nCols = 2*odd_nCols;
    even_nCols = floor((nCols-1)/2);
    even_2nCols = 2*even_nCols;
    odd_end = 2*odd_nCols;
    even_end = 2*even_nCols+1;
    for row = 1:nRows
        if mod(row, 2) == 1
            % Odd row: average samples 1 and 2, 3 and 4, etc.
            avgValue = repmat(mean(reshape(hexImgAvg_double(row,1:odd_end, :),[2,odd_nCols,3]), 1),2,1);
            hexImgAvg_double(row,1:odd_end, :) = reshape(avgValue,[1,odd_2nCols,3]);
        else
            % Even row: average samples 2 and 3, 4 and 5, etc.
            avgValue = repmat(mean(reshape(hexImgAvg_double(row,2:even_end, :),[2,even_nCols,3]), 1),2,1);
            hexImgAvg_double(row,2:even_end, :) = reshape(avgValue,[1,(even_2nCols),3]);
        end
    end
end

function pixelArray = hexDouble2Array(hexImgAvg_double)
    nRows = size(hexImgAvg_double,1);
    nCols = size(hexImgAvg_double,2);
    
    % Create a linear array of pixel values and their x and y coordinates
    pixelArray = [];
    for row = 1:nRows
        if mod(row, 2) == 1
            % Odd row: take odd columns
            cols = 1:2:nCols;
        else
            % Even row: take even columns
            cols = 2:2:nCols;
        end
        for col = cols
            pixelValue = squeeze(hexImgAvg_double(row, col, :));
            pixelArray = [pixelArray; col, row, pixelValue'];
        end
    end
end