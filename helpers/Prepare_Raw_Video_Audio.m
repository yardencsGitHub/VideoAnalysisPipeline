%% This script takes a list of FS videos, extracts and aligns video and audio data.
% Results are saved to the Raw data folder (on the backup disk)
% file names to process are in the variable 'keys' in the file
% 'FS_movies_list.mat'
FSfolder = '/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline';
addpath(genpath(FSfolder),'-end');
UtilsFolder = '/Users/yardenc/Documents/GitHub/small-utils';
addpath(genpath(UtilsFolder),'-end');
% OutputFolder = '/Volumes/home/Data/Imaging/lrb853_15/RawData';

% previos_annotation_file = 'lrb85315auto_annotation3';
% load(previos_annotation_file);
% old_keys = keys;
keysinfile = 'FS_movies_list';
% DIR = pwd;
load(keysinfile);
% keys = {old_keys{:} keys{:}};
% parameters:
ndel_frames = 0;

for fnum = 1:numel(keys)
    display([num2str(fnum) '/' num2str(numel(keys)) ' files'])
    tokens = regexp(keys{fnum},'_','split');
    if (str2num(tokens{4})==04) & (str2num(tokens{5})==19)
        display('corrupt date 04/19');
        continue;
    end
    filename = [keys{fnum}(1:end-3) 'mov'];
    [a_ts, a, v_ts, v] = extractmedia(filename);
    [video.width, video.height, video.channels] = size(v{1});
    video.times = v_ts;
    video.nrFramesTotal = size(v,1);
    video.FrameRate = 1/mean(diff(v_ts));
    for ii = 1: size(v,1)
        video.frames(:,:,:,ii) = v{ii};
    end
    audio.nrChannels = 1;
    audio.bits = 16;
    audio.nrFrames = length(a);
    audio.data = double(a);
    audio.rate = 48000;
    audio.TotalDurration = audio.nrFrames/48000;
    fs = 48000;
    params.audio = audio;
    params.video = video;
    clear audio; clear video;
    % Convert resample and cut audio and video to produce the structure Aud and
    % the 3D matrix (movie) vidMat
    if isfield(params.audio,'TotalDuration')
        Aud = params.audio;
    else
        Aud = struct('nrChannels',1,'bits',16,'nrFrames',params.audio.nrFrames,'data',params.audio.data,'rate',params.audio.rate,'TotalDuration',params.audio.TotalDurration);
    end
    vidTimes = 0:1/30:params.video.times(end);
    Vid_end = vidTimes(end)+1/30;
    if (Aud.nrFrames/Aud.rate > Vid_end)
        idx = find([0:(Aud.nrFrames-1)]/Aud.rate >= Vid_end);
        Aud.data(idx) = [];
        Aud.nrFrames = numel(Aud.data);
        Aud.TotalDuration = Aud.nrFrames/Aud.rate;
    end
    if (Aud.nrFrames/Aud.rate < Vid_end)
        samples_to_add = round(Aud.rate / 30 * numel(vidTimes) - Aud.nrFrames);
        Aud.data = [Aud.data; zeros(samples_to_add,1)];
        Aud.nrFrames = numel(Aud.data);
        Aud.TotalDuration = Aud.nrFrames/Aud.rate;
    end
    %Aud.annotation = elements{fnum};
    %(Aud.nrFrames/Aud.rate-1/30);
    [Y,n]=FS_Format(params.video.frames,1,1);
    [xx,yy,zz] = meshgrid(1:params.video.height,1:params.video.width,params.video.times);
    [xxx,yyy,zzz] = meshgrid(1:params.video.height,1:params.video.width,vidTimes);
    vidMat = double(Y - min(Y(:)));
    vidMat = interp3(xx,yy,zz,vidMat,xxx,yyy,zzz);
    clear params;
    % Trim video and audio
    startidx = round(Aud.rate * vidTimes(ndel_frames + 1)); 
    Aud.nrFrames = Aud.nrFrames - startidx;
    Aud.data(1:startidx) = [];
    Aud.TotalDuration = Aud.TotalDuration - startidx / Aud.rate;
%     Aud.annotation.segAbsStartTimes = Aud.annotation.segAbsStartTimes - startidx/48000;
%     Aud.annotation.segFileStartTimes = Aud.annotation.segFileStartTimes - startidx/48000;
%     Aud.annotation.segFileEndTimes = Aud.annotation.segFileEndTimes - startidx/48000;
%     Aud.annotation.phrases = return_phrase_times(Aud.annotation);
    vidMat(:,:,1:ndel_frames) = [];
    vidTimes(1:ndel_frames) = [];
    vidTimes = vidTimes - vidTimes(1);

    OutputFileName = ['RawData_' keys{fnum}(1:end-3) 'mat'];
    save(fullfile(OutputFolder,OutputFileName),'vidMat','vidTimes','Aud','-v7.3');
end