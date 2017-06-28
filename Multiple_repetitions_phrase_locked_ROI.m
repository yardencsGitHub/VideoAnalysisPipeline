%%
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/small-utils'));
n_del_frames = 5;
cd ('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav');
load('lrb85315auto_annotation3.mat');
keys_index = [];
for i=1:numel(keys)
    tokens = regexp(keys{i},'_','split');
    keys_index = [keys_index; str2num(tokens{2})];
end

%%
RawDIR = '/Volumes/CanaryData/DATA/lrb853_15/RawData';
FSfolder = '/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/FreedomScope/Analysis Pipeline';
addpath(FSfolder);

filename = 'lrb85315_1115_2017_03_27_06_59_53.mat';
RawName = fullfile(RawDIR,['RawData_' filename]);
load(RawName);

%%
filt_rad = 10; filt_sigma = 10; baseline_per=5; disp_band = [100 9e3];




W = convn(vidMat, single(reshape([1 1 1] / 2, 1, 1, [])), 'same');
W = imresize(W,1/2);
%smooth3(vidMat,'box', [1 1 3])*1.5;
[rows,columns,frames] = size(W);
h = fspecial('gaussian',filt_rad,filt_sigma);
baseline = imfilter(W,h,'circular','replicate');
baseline = repmat(prctile(baseline,baseline_per,3),[1 1 frames]);
dff2 = (abs((W.^2-baseline.^2)))./baseline;
h=fspecial('disk',2);
dff2=imfilter(dff2,h); 
H = prctile(max(max(dff2(:,:,:))),60);
L = 5;
clim = [double(L) double(H)];
NormIm = mat2gray(dff2, clim);    
I = imresize(std(NormIm,[],3),2);

[ROI]=FS_Image_ROI(uint8(I/max(I(:))*255)*1.5);

[rows,columns,frames] = size(vidMat);
vid_temp = reshape(vidMat,rows*columns,frames);
for roi_num = 1:length(ROI.stats)
    cr = ROI.coordinates{roi_num}(:,2)+(ROI.coordinates{roi_num}(:,1)-1)*rows;
    roi_signal = sum(vid_temp(cr,:));
    subplot(1+length(ROI.stats),1,roi_num+1);
    dff = (roi_signal - quantile(roi_signal,0.05))/quantile(roi_signal,0.05);
    plot(vidTimes,dff);
    axis tight;
end

%%
figure;
%load cr;
for fnum = 1:numel(keys)
    if ismember(5,elements{fnum}.segType)
    
        filename = fullfile(RawDIR,['RawData_' keys{fnum}(1:end-3) 'mat']);
        load(filename,'vidMat','vidTimes');
        [rows,columns,frames] = size(vidMat);
        v_dt = mean(diff(vidTimes));
        vid_shift = n_del_frames*v_dt;
    
        phrases = return_phrase_times(elements{fnum});
        phrases.phraseFileStartTimes = phrases.phraseFileStartTimes - vid_shift;
        phrases.phraseFileEndTimes = phrases.phraseFileEndTimes - vid_shift;
        vid_temp = reshape(vidMat,rows*columns,frames);
        roi_signal = sum(vid_temp(cr,:));
        dff = (roi_signal - quantile(roi_signal,0.05))/quantile(roi_signal,0.05);
        locs = find(phrases.phraseType == 5);
        for loc=1:numel(locs)
            idx = find((vidTimes >= phrases.phraseFileStartTimes(locs(loc))) & (vidTimes <= phrases.phraseFileEndTimes(locs(loc))));
            plot(dff(idx));
            hold on;
        end
    end
end
%%
prev_reg=[];
prev_times = [];
raster = [];
for fnum = 1:numel(keys)
    if ismember(303,elements{fnum}.segType)
        display(fnum);
        filename = fullfile(RawDIR,['RawData_' keys{fnum}(1:end-3) 'mat']);
        try
        load(filename,'vidMat','vidTimes');
        [rows,columns,frames] = size(vidMat);
        v_dt = mean(diff(vidTimes));
        vid_shift = n_del_frames*v_dt;
    
        phrases = return_phrase_times(elements{fnum});
        phrases.phraseFileStartTimes = phrases.phraseFileStartTimes - vid_shift;
        phrases.phraseFileEndTimes = phrases.phraseFileEndTimes - vid_shift;
        vid_temp = reshape(vidMat,rows*columns,frames);
        roi_signal = sum(vid_temp(cr,:));
        dff = (roi_signal - quantile(roi_signal,0.05))/quantile(roi_signal,0.05);
        zdff = zscore(dff);
        locs = find(phrases.phraseType == 303);
        for loc=1:numel(locs)
            idx = find((vidTimes >= phrases.phraseFileStartTimes(locs(loc))) & (vidTimes <= phrases.phraseFileEndTimes(locs(loc))));
            plot(zdff(idx));
            hold on;
            if locs(loc) > 2       
                idx_start = min(find(vidTimes >= phrases.phraseFileStartTimes(locs(loc))-10));
                pad_start = zeros(1,ceil((10 - phrases.phraseFileStartTimes(locs(loc)) + vidTimes(idx_start))/v_dt));
                idx_end = max(find(vidTimes <= phrases.phraseFileStartTimes(locs(loc))+10));
                pad_end = zeros(1,ceil((10 + phrases.phraseFileStartTimes(locs(loc)) - vidTimes(idx_end))/v_dt));
                raster = [raster; pad_start zdff(idx_start:idx_end) pad_end];
                prev_reg = [prev_reg ; phrases.phraseType(locs(loc)-2:locs(loc)-1)'];
                prev_times = [prev_times; [phrases.phraseFileStartTimes(locs(loc)-2) phrases.phraseFileEndTimes(locs(loc)-2) ...
                    phrases.phraseFileStartTimes(locs(loc)-1) phrases.phraseFileEndTimes(locs(loc)-1)] - phrases.phraseFileStartTimes(locs(loc))];
            end
        end
        catch em
            display('*');
        end
        
    end
end