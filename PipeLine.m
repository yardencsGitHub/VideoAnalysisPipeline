%% 1 Choose FS movies
[FileName,PathName,FilterIndex] = uigetfile([pwd '/*.mov'],'MultiSelect','on');
%% 2 Convert those movies to Matlab structure 
FS_AV_Parse_batch(pwd,FileName);
%%  Convert files to intensity video matrix format
SourceDir = pwd;
TargetDir = [pwd '/mat'];
[SUCCESS,MESSAGE,MESSAGEID] = movefile(fullfile(TargetDir,'*.mat'),SourceDir);
%%
matlist = {};
for fnum = 1:numel(FileName)
    fname = FileName{fnum}; fname(end-2:end) = 'mat';
    matlist{fnum} = fname;
    [vidMat, vidTimes, Aud] = FS2MAT(fname);
    save(fullfile(TargetDir,fname),'vidMat','vidTimes','Aud','-v7.3');
end
%% 3 remove first dark frames
RemoveInactiveFrames(matlist,'reduced_',25);
cd reduced;
d = dir('*.mat');
for fnum = 1:numel(d)
    load(d(fnum).name);
    write_mat_2_moviefile(vidMat,[d(fnum).name(1:end-4) '.mp4'],'MPEG-4',0);
end
% 4 align

%% 5 concatenate
OrigFolder = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced';
pref = 'reduced_';
TargetDir = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced';
cd(OrigFolder);
movlist = dir([pref '*.mat']);
[vidMat, vidTimes, Aud] = ConcatenateMovies({movlist.name},0.5);

%%
addpath(genpath('/Users/yardenc/Documents/GitHub/CNMF-E'));

CNMFE_Folder = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced/cnmfe';
if exist(CNMFE_Folder) == 0
    mkdir(CNMFE_Folder);
end
%factor = 0.5;
Y = uint8(vidMat);
%Y = uint8(imresize(double(Y),factor));
Ysiz = size(Y)';
save(fullfile(CNMFE_Folder,'cnmfe_source.mat'),'Y','Ysiz','Aud','vidTimes','-v7.3');

%% after running the cnmfe .. manual roi choosing (for now)
CaRaster = [neuron_full.C_raw(3,:)+neuron_full.C_raw(1,:); ...
neuron_full.C_raw(6,:); ...     
neuron_full.C_raw(7,:); ... 
neuron_full.C_raw(8,:)+neuron_full.C_raw(10,:); ...
neuron_full.C_raw(11,:); ... 
neuron_full.C_raw(12,:); ... 
neuron_full.C_raw(14,:); ... 
neuron_full.C_raw(15,:); ... 
neuron_full.C_raw(5,:)+neuron_full.C_raw(17,:); ...
neuron_full.C_raw(25,:); ... 
neuron_full.C_raw(29,:); ... 
neuron_full.C_raw(16,:)+neuron_full.C_raw(31,:); ...
neuron_full.C_raw(34,:); ... 
% sparser
neuron_full.C_raw(4,:); ... 
neuron_full.C_raw(2,:); ... 
neuron_full.C_raw(23,:); ... 
];
OrigFolder = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced';
save(fullfile(OrigFolder,'CaRaster.mat'),'CaRaster');
%% Now choose the template
[template,fs]=zftftb_select_template('/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced',@AudioLoad_func);
%% find all template segments
addpath(enpath('/Users/yardenc/Documents/GitHub/find-audio'));
cd('/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced/cnmfe');
audiowrite('template.wav',template,fs);
audiowrite('Aud.wav',Aud.data,fs);
threshold = threshold_for_find_audio;
[starts, ends, scores] = find_audio(Aud.data, template, fs, 'threshold_score', threshold);

%%
%% plot the mean and SEM of some ROIs signal aligne to template
figure;
n_neurons = 14;% size(neuron_full.C,1);
locs = 1:14; %[1 4 5 6 7 8 11 12 29 14]; %setdiff(1:n_neurons,[17 18 34 36 37 38]);
n_neurons = numel(locs);
for plotnum = n_neurons:-1:1
    h = axes;
    %%%
    ca_series = CaRaster(locs(plotnum),:);%+neuron_full.C_raw(3,:);
    raster = [];
    for i=1:85
        idx = find((vidTimes>=starts(i)-0.5) & (vidTimes <= ends(i)+1));
        raster = [raster; (zscore(ca_series(idx))) zeros(1,65-numel(idx))];
    end
    mn = mean(raster);
    se = std(raster)/sqrt(85);
    fill([1:65 65:-1:1],[mn-se fliplr(mn+se)],[0.5 0.5 0.5],'EdgeColor','none');
    hold on;
    plot(mean(raster),'b','LineWidth',2);
    line([16 16],[-1 1],'Color','r','LineWidth',1,'LineStyle','--')
    line([35 35],[-1 1],'Color','r','LineWidth',1,'LineStyle','--')
    axis tight
%%%
%     plot(neuron_full.C(locs(plotnum),:));
%     hold on;
%     plot(neuron_full.C_raw(locs(plotnum),:));
    h.XTick = [];
    h.YTick = [];
    pos = [0.1300    0.9213    0.7750    0.0037]; pos(2) = (plotnum-1)/n_neurons; pos(4) = 1/n_neurons;
    h.Position = pos;
    h.Box = 'off';
    h.Color = 'none';
    
end

%% create movie of denoised cells

temp_vid = [];
temp_aud = [];

for i=1:85
      idx = find((vidTimes>=starts(i)-0.5) & (vidTimes <= ends(i)+1));
      temp_vid = cat(3,temp_vid, ...
          reshape(neuron_full.A*neuron_full.C_raw(:,idx),240,320,numel(idx)));
      idx = max(ceil((starts(i)-0.5)*48000),1):min(floor((ends(i)+1)*48000),Aud.nrFrames);
      temp_aud = [temp_aud; Aud.data(idx)];
end
temp_aud = filtfilt(Hd.Numerator,1,double(temp_aud));
%% write the denoised movie to an AVI file 
cd('/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced/cnmfe');
writerObj = vision.VideoFileWriter('denoised.avi','AudioInputPort',true);
% total number of frames
nFrames   = size(temp_vid,3);
% assign FrameRate (by default it is 30)

writerObj.FrameRate =  30;

% length of the audio to be put per frame

val = size(temp_aud,1)/nFrames;

% Read one frame at a time

for k = 1 : nFrames
    % reading frames from a directory
    Frame = mat2gray(temp_vid(:,:,k)+35,[0 120]);
    % adding the audio variable in the step function
    step(writerObj,Frame,temp_aud(val*(k-1)+1:val*k,:)); % it is 2 channel that is why I have put (:)

end

% release the video

release(writerObj)

%% deconvolve ca signals
SpikeRaster = [];
for celln = 1: size(CaRaster,1)
    disp(celln);
    F = CaRaster(celln,:);
    F = resample(F,10,3);
    [c, s, options] = deconvolveCa(F,'method','thresholded','type','ar2');
    SpikeRaster = [SpikeRaster; s' > 1];
end
%% Raster Plots
figure;
ResVidTimes = [0:(size(SpikeRaster,2)-1)]/100; %resample(vidTimes,10,3);
n_neurons = 14;% size(neuron_full.C,1);
locs = 1:14; %[1 4 5 6 7 8 11 12 29 14]; %setdiff(1:n_neurons,[17 18 34 36 37 38]);
n_neurons = numel(locs);
for plotnum = n_neurons:-1:1
    h = axes;
    %%%
    
    raster = [];
    for i=1:85
        idx = find((ResVidTimes>=starts(i)-0.5) & (ResVidTimes <= ends(i)+1));
        raster = [raster; (SpikeRaster(plotnum,idx)) zeros(1,215-numel(idx))];
    end
    mn = mean(raster*100);
    se = std(raster*100)/sqrt(85);
    fill([1:215 215:-1:1],[mn-se fliplr(mn+se)],[0.5 0.5 0.5],'EdgeColor','none');
    hold on;
    plot(mn,'b','LineWidth',2);
    line([51 51],[0 max(mn)],'Color','r','LineWidth',1,'LineStyle','--')
    line([117 117],[0 max(mn)],'Color','r','LineWidth',1,'LineStyle','--')
    axis tight
    
%%%
%     plot(neuron_full.C(locs(plotnum),:));
%     hold on;
%     plot(neuron_full.C_raw(locs(plotnum),:));
    h.XTick = [];
    h.YTick = [];
    yticks(ceil(max(mn)/2))
    pos = [0.1300    0.9213    0.7750    0.0037]; pos(2) = (plotnum-1)/n_neurons; pos(4) = 1/n_neurons;
    h.Position = pos;
    h.Box = 'off';
    h.Color = 'none';
    
end

%% Create 'time' labels by dividing the template +-0.1sec to equal bins
ResVidTimes = [0:(size(SpikeRaster,2)-1)]/100;
labels = cell(6,1);
for i=1:85
    idx = find((ResVidTimes>=starts(i)-0.05) & (ResVidTimes <= ends(i)+0.05));
    split = splitidx(idx,6);
    for j = 1:6
        labels{j} = [labels{j} split{j}];
    end
end
%% plot a sample of time labels
idx = find((ResVidTimes>140.9) & (ResVidTimes<143));
a_idx = 140.9:1/48000:143;
figure; plot(a_idx,filtered(140.9*48000:143*48000));

hold on
zer = zeros(1,numel(ResVidTimes)); zer(labels{1}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'r');
zer = zeros(1,numel(ResVidTimes)); zer(labels{2}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'g');
zer = zeros(1,numel(ResVidTimes)); zer(labels{3}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'b');
zer = zeros(1,numel(ResVidTimes)); zer(labels{4}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'m');
zer = zeros(1,numel(ResVidTimes)); zer(labels{5}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'c');
zer = zeros(1,numel(ResVidTimes)); zer(labels{6}) = 0.5;
plot(ResVidTimes(idx),zer(idx),'k');

%%
model = maxent.createModel(size(SpikeRaster,1),'ising');