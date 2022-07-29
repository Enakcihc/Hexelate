% hex-coordinate functions based on redblobgames
% https://www.redblobgames.com/grids/hexagons
% Artistic style inspired by two artists residing in Rotterdam, Netherlands
% whose names I need to figure out to give proper credit.
%
% Written by Aaron B. Wong (2022)
%
% Written in Matlab R2019a.

%% 
% Img = imread('private\sampImage2.JPG');
% Img = imread('private\Family2.JPG');

Img = imread('private\373_NZ68521-bewerkt_cropped.jpg'); downsamp = 20;

% Img = imread('Colourful_Flower_01.jpg'); downsamp = 100;


%% Coordinate of image pixels
xDim = size(Img,2); yDim = size(Img,1);
nPixel = xDim * yDim;
xx = 1:xDim; yy = 1:yDim;
xxyy = combvec(xx,yy); % xy coordinates of all pixel in image. Upper left pixel: (1,1)
qrs = p2h(xxyy); % pixel positions in cubic coordinates

%% calculate grid that spans the image

gridHex = unique(cube_round(qrs)','row')';
nHexel = size(gridHex,2);
gridXY = h2p(gridHex);
% figure; 
% imagesc(Img);hold on; 
% scatter(gridXY(1,:),gridXY(2,:),'.')
disp(['Number of pixels: ',num2str(nPixel), ' (',num2str(xDim),' x ',num2str(yDim),')'])
disp(['Number of hexels: ',num2str(nHexel)])
%% calculate grid for down sampling

% expand grid
sampHex = (gridHex) * downsamp; 
sampXY = h2p(sampHex);
% exclude out of bound hexels
sel = sampXY(1,:) <= xDim & sampXY(2,:) <= yDim & sampXY(1,:) > 0 & sampXY(2,:) > 0;
sampHex = sampHex(:,sel);sampXY = sampXY(:,sel);
% get column-row (odd-r) representation
sampCR = ax2oddr(sampHex./downsamp);
nSamp = size(sampHex,2);
disp(['Number of downsampled hexels: ',num2str(nSamp)])

% calculate axis labels
sel = sampCR(2,:) == min(sampCR(2,:));
ColLbl_x = sampXY(1,sel);
ColLbls = sampCR(1,sel);
sel = sampCR(1,:) == min(sampCR(1,:));
RowLbl_x = sampXY(2,sel);
RowLbls = sampCR(2,sel);
[RowLbls,idx] = sort(RowLbls);
RowLbl_x = RowLbl_x(idx);

clear('sel','idx');


%% sampling pixel intensity
% & generate pixelation image
nCh = size(Img,3);
CData_samp = nan(size(sampXY,2),nCh);
CData_mean = nan(size(sampXY,2),nCh);
HexelateImg_samp = nan(size(Img));
HexelateImg_mean = nan(size(Img));

% nearest central pixel
tic
for ii = 1:size(sampXY,2)
    CData_samp(ii,:) = squeeze(Img(round(sampXY(2,ii)),round(sampXY(1,ii)),:));
end
disp('Sampling done.')
toc

% mean 
tic

mapHex = downsamp*cube_round(qrs./downsamp);
for ii = 1:size(sampXY,2)
    sel = all(sampHex(:,ii) == mapHex,1);
%     subplot(1,3,2)
%     scatter(xxyy(1,sel),xxyy(2,sel))
    nn = sum(sel);
    for ch = 1:nCh
        linearInd = sub2ind([yDim,xDim,nCh], xxyy(2,sel), xxyy(1,sel),ch*ones(1,nn));
        CData_mean(ii,ch) = mean(Img(linearInd),'all');
        HexelateImg_samp(linearInd) =CData_samp(ii,ch);
        HexelateImg_mean(linearInd) =CData_mean(ii,ch);
    end
end
disp('Averaging done.')
toc

%% output to JSON

ch = 2;
fileID = fopen('Hexelated.json','w');
JSON = jsonencode(table(sampXY(1,:)',sampXY(2,:)',CData_mean(:,ch),...
                'VariableNames',{'x','y','c'}));
fwrite(fileID,JSON);
fclose(fileID);
%% Show images (R,G,B, RGB)
fig = figure;
fig.Position = [300,50,750,750];
Colors = {'r','g','b'};%{'c','m','y'};
nSteps = 8;
for ch = 1:nCh
    subplot(2,2,ch)
        sz = 10*szLU(-CData_mean(:,ch),nSteps);%round((280-CData_mean(:,ch))./64)+eps;

%     scatter(sampXY(1,:),sampXY(2,:),sz,Colors{ch},'filled');
    scatter(sampXY(1,:),sampXY(2,:),sz,'filled');
    set(gca,'YDir','reverse');
%     axis square
        ylim([0,yDim+1]);xlim([0,xDim+1]);
    title(Colors{ch});
end
subplot(2,2,4)
scatter(sampXY(1,:),sampXY(2,:),9,CData_mean./255,'filled');
set(gca,'YDir','reverse');
    ylim([0,yDim+1]);xlim([0,xDim+1]);

% axis square

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
% saveas(fig,'sampOutput.pdf')
%% demo
fig = figure;
fig.Position = [200,300,900,300];

nSteps = 8;

% original image
subplot(1,3,1)
    imagesc(Img);
    axis square
    
% scaled
    ch = 3; % channel to use
    dotScale = 2;
    useMean = 1;
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
    scatter(sampXY(1,:),sampXY(2,:),dotScale*sz,c,'filled');
    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square

    subplot(1,3,3)
    scatter(sampXY(1,:),sampXY(2,:),15,c,'filled');
%         scatter(sampXY(1,:),sampXY(2,:),sz_mean,'filled');

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
[~,locExtremes(1)] = min([1,1]*sampXY,[],2);
[~,locExtremes(2)] = min([1,-1]*sampXY,[],2);
[~,locExtremes(3)] = max([1,1]*sampXY,[],2);
[~,locExtremes(4)] = max([1,-1]*sampXY,[],2);

% locExtremes = union(locMinCR,locMaxCR);
for ii = 1:length(locExtremes)
    idx = locExtremes(ii);
    text(sampXY(1,idx),sampXY(2,idx),['(',num2str(sampCR(1,idx)),',',num2str(sampCR(2,idx)),')'],...
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
    scatter(sampXY(1,:),sampXY(2,:),dotScale*sz,'filled');
    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square
    subplot(1,2,2)
    scatter(sampXY(1,:),sampXY(2,:),dotScale*nSteps/2,c,'filled');

    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    
    
    axis square
colormap(CScale);
caxis([min(uSz)-.5,max(uSz)+.5])
colorbar('Ticks',uSz,'Position',[.92,.3,.02,.5]);

% label extremes
locExtremes = nan(4,1);
[~,locExtremes(1)] = min([1,1]*sampXY,[],2);
[~,locExtremes(2)] = min([1,-1]*sampXY,[],2);
[~,locExtremes(3)] = max([1,1]*sampXY,[],2);
[~,locExtremes(4)] = max([1,-1]*sampXY,[],2);

% locExtremes = union(locMinCR,locMaxCR);
for ii = 1:length(locExtremes)
    idx = locExtremes(ii);
    text(sampXY(1,idx),sampXY(2,idx),['(',num2str(sampCR(1,idx)),',',num2str(sampCR(2,idx)),')'],...
        'FontSize',8);
end

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
% saveas(fig,'instruction.pdf')
%% Hexelate Images
fig = figure;
fig.Position = [200,300,900,250];
    subplot(1,3,1)
    imagesc(Img); axis square;
    subplot(1,3,2)
    image(HexelateImg_samp./255); axis square;
    subplot(1,3,3)
    image(HexelateImg_mean./255); axis square;

%%
function [xy] = h2p(qrs)
    % pointy top hex
    % hex: 2x1 vector containing [q,r]
    A = [   sqrt(3),    0.5*sqrt(3);...
            0,          1.5         ];
    xy = A * qrs(1:2,:);
end

function qrs = p2h(xy)
     % pointy top hex
    A = [   sqrt(3)/3,   -1/3;...
            0,          2/3        ];
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
function hex = oddr2ax(cr)
    q = cr(1,:) - (cr(2,:) - bitand(cr(2,:),1)) / 2;
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