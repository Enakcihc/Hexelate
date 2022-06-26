% hex-coordinate functions based on redblobgames
% https://www.redblobgames.com/grids/hexagons
% Artistic style inspired by two artists residing in Rotterdam, Netherlands
% whose names I need to figure out to give proper credit.
%
% written by Aaron B. Wong (2022)

%% 
% Img = imread('sampImage2.JPG');
Img = imread('Family2.JPG');
%% Coordinate of image pixels
xDim = size(Img,2); yDim = size(Img,1);
xx = 1:xDim; yy = 1:yDim;
xxyy = combvec(xx,yy); % xy coordinates of all pixel in image. Upper left pixel: (1,1)
qrs = p2h(xxyy); % pixel positions in cubic coordinates

%% calculate grid
cr = ax2oddr(round(qrs));
colMin =min(cr(1,:));
colMax =max(cr(1,:));
rowMin =min(cr(2,:));
rowMax =max(cr(2,:));
%%
ccc = 0:colMax;%10;
rrr = 0:rowMax;%10;
gridCR = combvec(ccc,rrr);
gridHex = oddr2ax(gridCR);
gridXY = h2p(gridHex);

% down sampling
downsamp = 10;
sampHex = (gridHex) * downsamp;
sampXY = h2p(sampHex);
sel = sampXY(1,:) < xDim & sampXY(2,:) < yDim & sampXY(1,:) > 0 & sampXY(2,:) > 0;
sampHex = sampHex(:,sel);
sampCR = ax2oddr(sampHex);
sampXY = sampXY(:,sel);
%% sampling pixel intensity
nCh = size(Img,3);
CData_samp = nan(size(sampXY,2),nCh);
CData_mean = nan(size(sampXY,2),nCh);
% nearest central pixel
for ii = 1:size(sampXY,2)
    CData_samp(ii,:) = squeeze(Img(round(sampXY(2,ii)),round(sampXY(1,ii)),:));
end

% mean (to be done)
mapHex = downsamp*cube_round(qrs./downsamp);
for ii = 1:size(sampXY,2)
    sel = all(sampHex(:,ii) == mapHex,1);
    for ch = 1:nCh
        CData_mean(ii,ch) = mean(Img(xxyy(2,sel),xxyy(1,sel),ch),'all');
    end
end

%% Show images (R,G,B, RGB)
fig = figure;
fig.Position = [300,50,750,750];
Colors = {'r','g','b'};%{'c','m','y'};
nSteps = 8;
for ch = 1:3
    subplot(2,2,ch)
        sz = 10*szLU(-CData_mean(:,ch),nSteps);%round((280-CData_mean(:,ch))./64)+eps;

%     scatter(sampXY(1,:),sampXY(2,:),sz,Colors{ch},'filled');
    scatter(sampXY(1,:),sampXY(2,:),sz,'filled');
    set(gca,'YDir','reverse');
    axis square
        ylim([0,yDim+1]);xlim([0,xDim+1]);
    title(Colors{ch});
end
subplot(2,2,4)
scatter(sampXY(1,:),sampXY(2,:),9,CData_mean./255,'filled');
set(gca,'YDir','reverse');
    ylim([0,yDim+1]);xlim([0,xDim+1]);

axis square

fig.PaperUnits = 'inches';
fig.PaperSize = fig.Position(3:4)./96; %96 dpi
saveas(fig,'sampOutput.pdf')
%% demo
fig = figure;
fig.Position = [200,300,900,300];


    subplot(1,3,1)
    imagesc(Img);
    axis square
    
ch = 3;
    sz = round((280-CData_samp(:,ch))./32);
    nColors = max(unique(sz));
    CScale = jet(nColors);
    c = CScale(sz,:);
    subplot(1,3,2)
    scatter(sampXY(1,:),sampXY(2,:),sz,'filled');
    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square

    sz_mean = round((280-CData_mean(:,ch))./32);

    subplot(1,3,3)
%     scatter(sampXY(1,:),sampXY(2,:),sz,c,'filled');
        scatter(sampXY(1,:),sampXY(2,:),sz_mean,'filled');

    set(gca,'YDir','reverse');
    ylim([0,yDim]);xlim([0,xDim]);
    axis square
%% test functions
% for ii = 1:size(gridXY,2)
ii = 1;
    text(sampXY(1,ii),sampXY(2,ii),...
        sprintf('(%.1f,%.1f,%.1f)',sampHex(1,ii),sampHex(2,ii),sampHex(3,ii))...
        )
% end
%%

for jj = (5*xDim)+(1:xDim*20)
    
    text(xxyy(1,jj),xxyy(2,jj),...
        sprintf('(%.0f,%.0f,%.0f)',mapHex(1,jj),mapHex(2,jj),mapHex(3,jj)),...
        ...sprintf('(%.1f,%.1f,%.1f)',qrs(1,jj),qrs(2,jj),qrs(3,jj)),...
        'Color','w')
end

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
    
    sz = round ( (nSteps - 1) .* (vals - minV) ./ (maxV-minV ) ) ./ (nSteps-1) + eps;
end