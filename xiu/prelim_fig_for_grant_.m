clear all; close all; clc

%% set path
addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\ThorLabs MatlabScripts'));
addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\imgio-matlab-master'));
addpath(genpath('C:\Users\xiuye\Dropbox\Suite2P'));
% rmpath(genpath('C:\Users\xiuye\Dropbox\Suite2P_'));
addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\Suite2P_xiu'));

%% open GUI
% new_main; % manually load the file below, press button 'export' => h in workspace
h = load('C:\Users\xiuye\Documents\2P_processed\890C\122117\3\F_890C_122117_plane1_proc_020118.mat');

%% plot mean image with ROI's drawn on top.
% mean image
im1_0 = squeeze(h.dat.mimg(:,:,2));

% ROI image
Sat1     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
Sat2     =  ones(h.dat.cl.Ly, h.dat.cl.Lx);
H1              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);
H2              = zeros(h.dat.cl.Ly, h.dat.cl.Lx);


% chose cells
for i = 1:length(h.dat.stat)
    h.dat.stat(i).iscell = (ismember(i,cIX));
end
%
[iclust1, iclust2, V1, V2] = ...
    getviclust(h.dat.stat, h.dat.cl.Ly,  h.dat.cl.Lx, h.dat.cl.vmap, h.dat.F.ichosen);

% iselect     = iclust1==h.dat.F.ichosen;
% Sat1(iselect)= 0;
% 
% iselect     = iclust2==h.dat.F.ichosen;
% Sat2(iselect)= 0;

% given cIX;
nColors = length(cIX);
c = linspace(0.05,0.9,nColors);
cl = zeros(1,length(h.dat.stat));
cl(cIX) = c;
H1(iclust1>0)   = cl(iclust1(iclust1>0));
% H1(iclust1>0)   = h.dat.cl.rands(iclust1(iclust1>0));
% H2(iclust2>0)   = h.dat.cl.rands(iclust2(iclust2>0));

I = hsv2rgb(cat(3, H1, Sat1, V1));
im2 = min(I, 1);

inew = find(H1>0);
% combine and draw
im1_1 = mat2gray(im1_0);
im1 = repmat(im1_1,1,1,3);
Low_High  = [0,0.5];
im1 = imadjust(im1,Low_High);

im3 = (im1+im2)/2;
im3(inew) = im2(inew);
numpix = numel(H1);
im3(inew+numpix) = im2(inew+numpix);
im3(inew+2*numpix) = im2(inew+2*numpix);

figure;
imagesc(im3);
% Low_High  = [0,0.5];
% im4 = imadjust(im3,Low_High);
% imshow(im4) 
axis equal; axis off

for i = 1:length(cIX)
    ichosen = cIX(i);
    x0 = mean(h.dat.stat(ichosen).xpix);
    y0 = mean(h.dat.stat(ichosen).ypix);
    clr = squeeze(hsv2rgb(cl(ichosen),1,1));
    text(x0-10,y0,num2str(i),'color','w','HorizontalAlignment','right');%clr);
end
i = 5;
    ichosen = cIX(i);
    x0 = mean(h.dat.stat(ichosen).xpix);
    y0 = mean(h.dat.stat(ichosen).ypix);
    clr = squeeze(hsv2rgb(cl(ichosen),1,1));
    text(x0+10,y0,num2str(i),'color','w','HorizontalAlignment','left');%clr);

%% rank traces

T = struct2table(h.dat.stat);
M = zscore(h.dat.Fcell{1},0,2);
IsCell = table2array(T(:,28));
IX_ROI = find(IsCell);

i_start = 2000;
i_stop = 3000;

% calculate skewness
sk = zeros(length(IX_ROI),1);
for i = 1:length(IX_ROI)
    ichosen = IX_ROI(i);
    F =  h.dat.Fcell{1}(ichosen, i_start:i_stop);
    Fneu = h.dat.FcellNeu{1}(ichosen, i_start:i_stop);


% F(:, ops.badframes)  = F(:,    indNoNaN(ix));
% Fneu(:, ops.badframes)  = Fneu(:, indNoNaN(ix));


coefNeu = 0.7 * ones(1, size(F,1));

dF                  = F - bsxfun(@times, Fneu, coefNeu(:));

% dF          = F - Fneu;

sd           = std(dF, [], 2);
sdN          = std(Fneu, [], 2);

sk(i, 1) = skewness(dF, [], 2);
end

[~,IX_sort] = sort(sk,'ascend');
cIX = IX_ROI(IX_sort(11:30));

%% draw traces


figure; hold on;
pad = 1;


Fs = 15.3;            % Sampling frequency
xv = (1:i_stop-i_start+1)/Fs;
for i = 1:length(cIX)
    ichosen = cIX(i);
    F = [];
    Fneu = [];
    for j = 1:numel(h.dat.Fcell)
        F    = cat(2, F, h.dat.Fcell{j}(ichosen, :));
        Fneu = cat(2, Fneu, h.dat.FcellNeu{j}(ichosen, :));
    end
    
    y = double(F(i_start:i_stop));%my_conv_local(medfilt1(double(F), 3), 3);
    y_m = y-mean(y);
    y_n = y_m/max(y_m);%zscore(y);
    
    clr = squeeze(hsv2rgb(cl(ichosen),1,0.8));
    plot(xv,y_n - pad*(i-1),'color',clr);%[0.5,0.5,0.5])
    axis tight
    text(-1, 0.2- pad*(i-1),num2str(i),'HorizontalAlignment','right','color',[0.5,0.5,0.5]);
end

% draw scale bar
base_y = -pad*length(cIX);
plot([0,5],[base_y,base_y],'linewidth',2,'color','k');
text(2.5,base_y-0.6,'5 sec','HorizontalAlignment','center')
xlabel('sec')
ylim([base_y-1,2])
% set(gca,'xcolor','w');

axis off

%% draw stimulus
% lines
for i = 1:length(midPoints)

    x = time(midPoints(i));
    plot([x,x],[-30,0],'k:');
end
% bar
plot(frameTime,stimulus+1,'color','r');

%% 

load('C:\Users\xiuye\Documents\2P_processed\890C\122117\3\signals.mat','frameTime','stimulus');

%%
imInfo = xml2struct('C:\Users\xiuye\Documents\2P_rawdata\890C\122117\3\Experiment.xml');
% S.frameRate
frameRate = 15.2530;

%%
numTimepoints = 3000;
clockRate = 20000000;%17900000;
filename = fullfile('C:\Users\xiuye\Documents\2P_rawdata\890C\122117\3', 'Episode001.h5');
info = h5info(filename);

time = double(h5read(filename, '/Global/GCtr')')./ clockRate;
frameOut = h5read(filename, '/DI/Frame Out')';
strobes = h5read(filename, '/DI/Strobe')';

% Zero out strobes during calibration sequence.
strobeUp_0 = find(diff(strobes));
d_sUp_0 = diff(strobeUp_0);
IX_cut = find(d_sUp_0<mean(d_sUp_0));
ix_start = 1+length(IX_cut);
strobeUp = strobeUp_0(ix_start:end);
% strobes(1:1.32e5) = 0;
% strobeUp = find(diff(strobes));

frameOutUp = find(diff(frameOut));
d_FOup = diff(frameOutUp);
sample_int = 1970; % mode
samplespersec = frame_int*frameRate;

% time = double(h5read(filename, '/Global/GCtr')')./ samplespersec;

%%
midPoints = zeros(1, length(strobeUp)-1);
for ii = 1:length(strobeUp)-1
    a = strobeUp(ii); b = strobeUp(ii+1);
    midPoints(ii) = uint32(a + (b-a)/2);
end
ae = zeros(1, length(frameOut));
for ii = 1:length(midPoints)-1
    a = midPoints(ii); b = midPoints(ii+1);
    m = 1/(b-a);
    X = 1:(b-a);
    Y = m*X;
    ae(a:b-1) = Y;
end

stimulus = zeros(1, numTimepoints);
frameTime = zeros(1, numTimepoints);
frameOutUp = find(diff(frameOut));
% frameOutUp = frameOut(frameOutUp_IX);
for ii = 1:numTimepoints
    stimulus(ii) = ae(frameOutUp(ii));
    frameTime(ii) = time(frameOutUp(ii));
end

% Shift 't' so t0 is at first frame capture.
% frameTime = frameTime - frameTime(1);

%%

% period_sec = 3.867; % sec
% % period = round(Fs*period_sec);
% period = Fs*period_sec;
% n_rep = floor(length(y)/period);

% x_stim = (1:n_rep)*period;

