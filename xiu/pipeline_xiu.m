%%
clear all; close all; clc

%% set path
addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\ThorLabs MatlabScripts'));
addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\imgio-matlab-master'));
addpath(genpath('C:\Users\xiuye\Dropbox\Suite2P'));
% addpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\fromScott'));
% rmpath(genpath('C:\Users\xiuye\Dropbox\!Research\2Pcode\fromScott'));

rawfile_path = 'C:\Users\xiuye\Documents\2P_rawdata';
tiff_path = 'C:\Users\xiuye\Documents\2P_inputdata';

%% convert file format
t1 = tic;
mouse_str = '890C';%'m011';
date_str = '2018-02-09'%'122117';%'010918';

% expPath = 'm011\010918';
range_session = 1:2;%1:5;%1:4;
% SessionPaths = {'890C\122117\1','890C\122117\2','890C\122117\3','890C\122117\4','890C\122117\5'};

%% process rawdata into tiff format i.e. inputdata
for i_session = 1:length(range_session)
    disp(['i_session = ',num2str(i_session)]);
    
    %%
    sessionPath = fullfile(mouse_str,date_str,num2str(range_session(i_session)));
    inputDir = fullfile(rawfile_path,sessionPath);
    
    %% parse folder
    cd(inputDir);
    folder1 = ls('Un*');
    folder2 = ls('Sy*');
    
    %% load raw data
    filename = fullfile(inputDir,folder1,'Experiment.xml');
    % filename = 'C:\Users\xiuye\Documents\2P_rawdata\890C\121917\1\test0_007\Experiment.xml';
    [imData,imInfo,imFlyback,numStoredFrames] = loadThorlabsExperimentRaw(filename);
    
    
    %% save input data as tiff series
    % t2 = tic;
    % parentDir = 'C:\Users\xiuye\Documents\2P_inputdata\890C\121917\1';
    outputDir = fullfile(tiff_path,sessionPath);
    mov = squeeze(imData);
    mov_crop  = mov(:,:,1:numStoredFrames);
    
    %%
    if false
        tiffseries_save(outputDir, mov_crop);
        clear mov; clear imData; clear mov_crop;
    else % alternatively, downsample, crop in time and maybe try to save as tiff stack
        mov_resized = downsampleMovie(mov_crop, 0.5);
        tiffseries_save(outputDir, mov_resized);
        clear mov; clear imData; clear mov_resized; clear mov_crop;
    end
    % et2 = toc(t2);
    
end

et1 = toc(t1);

%% run Suite2p defaults
t3 = tic;
db0 = make_db_xiu(mouse_str,date_str,range_session);
master_file;
et3 = toc(t3);

%% run concat version

%% run Suite2p defaults
t3 = tic;
isconcat = 1;
db0 = make_db_xiu(mouse_str,date_str,range_session,isconcat);
master_file;
et3 = toc(t3);