%% This function takes a single FS videos, extracts and aligns video and audio data.

%%
function [vidMat, vidTimes, Aud] = Prepare_SingleFile_Raw_Video_Audio(path_to_file)
% Input:
%   path_to_file - full path to .mov file
%   Note: the code expects a 'TweetVision' style for the '_'-separated file
%   name (birdname_#file_yyyy_mm_dd_hh_mm_ss.mov
%   e.g. lrb85315_7811_2017_06_28_06_15_32.mov
% Outputs:
%   VidMat - the 3d 30Hz aligned video
%   vidTimes - vector of video time frames
%   Aud - the audio structure:
%       -- Aud = struct('nrChannels',#channels,'bits',16,'nrFrames',#frames,'data',wav data,'rate',sampling rate,'TotalDuration',audio duration (sec));

FSfolder = '/full/path/to/finch_scope/github/repo'; 
% this should allow access to ectractmedia.m from
% https://github.com/gardner-lab/FinchScope/tree/master/Analysis%20Pipeline
% e.g. FSfolder ='/Users/yardenc/Documents/GitHub/FinchScope/Analysis Pipeline';
addpath(genpath(FSfolder),'-end');

% ndel_frames allows deleting frames.
% Currently set to 0
ndel_frames = 0;
[path_str,fname,ext] = fileparts(path_to_file);

tokens = regexp(fname,'_','split');

[a_ts, a, v_ts, v] = extractmedia(path_to_file);
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

[Y,n]=FS_Format(params.video.frames,1);
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

vidMat(:,:,1:ndel_frames) = [];
vidTimes(1:ndel_frames) = [];
vidTimes = vidTimes - vidTimes(1);
end