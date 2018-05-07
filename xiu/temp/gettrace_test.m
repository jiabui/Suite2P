
%% manually examine 1 trace
L = load('C:\Users\xiuye\Documents\2P_processed\m011\120717\1\F_m011_120717_plane1.mat');
T = struct2table(L.stat); 
I = T(:,23); % field: 'std'
a = table2array(I);
find(a>10.21 & a<10.23)
%%
ix = 383;
ts = L.Fcell{1}(ix,:);
figure;plot(ts)



%% get params for experiment
% C:\Users\xiuye\Documents\2P_rawdata\m011\120717\exp1_1.3x\ThorImage_005
imInfo = xml2struct('Experiment.xml');

I = imInfo.ThorImageExperiment;
S = [];
% laser
S.gain = str2double(I.PMT.Attributes.gainA);
S.power = str2double(I.Pockels{1}.Attributes.start);
S.power2 = str2double(I.Pockels{1}.Attributes.stop);
S.pockelsMaxV = str2double(I.Pockels{1}.Attributes.pockelsMaxV);
S.pockelsMinV = str2double(I.Pockels{1}.Attributes.pockelsMinV);

% frame rate
S.frameRate = str2double(I.LSM.Attributes.frameRate);
% zoom 
S.mag = str2double(I.Magnification.Attributes.indexOfRefraction);
S.pixelSizeUM = str2double(I.LSM.Attributes.pixelSizeUM);
S.widthUM = str2double(I.LSM.Attributes.widthUM);
S.heightUM = str2double(I.LSM.Attributes.heightUM);
% image size
S.pixelX = str2double(I.LSM.Attributes.pixelX);
S.pixelY = str2double(I.LSM.Attributes.pixelY);
% ? extra
S.extClockRate = str2double(I.LSM.Attributes.extClockRate);
S.horizontalFlip = str2double(I.LSM.Attributes.horizontalFlip);
S.notes = imInfo.ThorImageExperiment.ExperimentNotes;

%% load ThorSync data
LoadSyncEpisode;
%%
% User selected C:\Users\xiuye\Documents\2P_rawdata\m011\120717\exp1_1.3x\SyncData005\Episode001.h5
x = int32(FrameTrigger);
x_start = find(diff(FrameTrigger)>0); % start / step up

% to compute find(diff(FrameTrigger)==-1), because uint32 is positive only,
% manually 'flip' and shift up as a hack
% x = FrameTrigger;
% x(x==0) = 2;
% x_stops = find(diff(x)==1);
% x_stop = x_stops(end);

x_stops = find(diff(x)==1);

%% Strobe
strb_start = find(diff(Strobe)>0);
temp = diff(strb_start);
interval_frames = mean(temp(7:end));
interval = time(strb_start(end))-time(strb_start(end-1));

%% k-means of all ROI traces
numK = 10;
M = h.dat.Fcell{1};
% M = h.dat.FcellNeu{1};
groupIX = kmeans(M,numK,'Replicates',5);
% [groupIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);

% sort data matrix based on clustering
[~,I] = sort(groupIX);
im = M;
im = im(I,:);

figure;imagesc(im)

%% Load .mat
dat = load('C:\Users\xiuye\Documents\2P_processed\F_m011_120717_plane1.mat');

% load reg file -> ops1

%%
badframes = find(ops1{1}.badframes);
im_bad = im;
im_bad(:,badframes) = 0;
figure;imagesc(im_bad)

%% find indices of selected ROI's

I = T(:,28); % field: 'std'
a = table2array(I);
IX_iscell = find(a);

%% play video - selected ROI/time-range, from (unregistered) raw images

range_f = 401:600;
ichosen = h.dat.F.ichosen;

[iclust1, iclust2, V1, V2] = ...
    getviclust(h.dat.stat, h.dat.cl.Ly,  h.dat.cl.Lx, h.dat.cl.vmap, h.dat.F.ichosen);

% make mask of chosen ROI
ipix    = h.dat.stat(ichosen).ipix;
[I,J] = ind2sub([1000,1005],ipix);

x_ = max(I)-min(I);
y_ = max(J)-min(J);
x1 = min(I)-x_;
x2 = max(I)+x_;
y1 = min(J)-y_;
y2 = max(J)+y_;
s1 = x2-x1+1;
s2 = y2-y1+1;

rectMask = zeros(1000,1005);
rectMask(x1:x2,y1:y2) = 1; 
rectMask = logical(rectMask);

mask_pad = zeros(size(ops.mimg));
mask_pad(ops.yrange,ops.xrange) = rectMask;
mask_pad = logical(mask_pad);

%%
im = ops.mimg;%iclust1;
im2 = im(mask_pad);
im3 = reshape(im2,[s1,s2]);
figure;imagesc(im3);axis equal; axis off

%% get raw images
figure;
imdir = 'C:\Users\xiuye\Documents\2P_inputdata\m011\120717\1';

cmap = colormap(gray);
nFrames = length(range_f);
IM = zeros(s1,s2,nFrames);
for i_count = 1:nFrames
    ii = range_f(i_count);
    filename = [num2str(ii),'.tiff'];
    C = imread(fullfile(imdir,filename));

    im = reshape(C(mask_pad),[s1,s2]);

    imagesc(im,[0,10843]);colormap gray
    drawnow;
    
    IM(:,:,i_count) = im;
%     f = getframe;
%     F(i_count) = f.cdata;
end

figure;imagesc(mean(IM,3));axis equal; axis off

%% write tiff stack
% display each plane and save as tif
tiffdir = 'test1_reg.tif';

I = mat2gray(IM);
% h = figure;colormap(gray)
for i_frame = 1:nFrames
    im = I(:,:,i_frame);
%     image(im);axis equal; axis off
%     drawnow;
    % save tiff
    if (i_frame == 1)
        imwrite(im, tiffdir, 'compression','none','writemode','overwrite')
    else
        imwrite(im, tiffdir, 'compression','none','writemode','append')
    end
    %     pause(0.2)
end
% close(h)

%% get registered images
imfile = 'C:\Users\xiuye\Documents\2P_reg\m011\120717\1\Plane1\120717_1_m011_2P_plane1_1.tif';
% C2 = imread(imfile,2);
cmap = colormap(gray);

nFrames = length(range_f);
IM = zeros(s1,s2,nFrames);
for i_count = 1:nFrames
    ii = range_f(i_count);
    filename = [num2str(ii),'.tiff'];
    C = imread(imfile,ii);

    im = reshape(C(mask_pad),[s1,s2]);

    imagesc(im,[0,10843]);colormap gray
    drawnow;
    
    IM(:,:,i_count) = im;
%     f = getframe;
%     F(i_count) = f.cdata;
end

figure;imagesc(mean(IM,3));axis equal; axis off