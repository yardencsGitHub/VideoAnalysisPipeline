function [vidMat, vidTimes, Aud] = FS2MAT(movfile)
% This also interpolates all videos to 30Hz.
audio_field = 'audio';
FSfolder = '/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline';
addpath(FSfolder);
params = load(movfile);
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
 %(Aud.nrFrames/Aud.rate-1/30);
[Y,n]=FS_Format(params.video.frames,1,1);
[xx,yy,zz] = meshgrid(1:params.video.height,1:params.video.width,params.video.times);
[xxx,yyy,zzz] = meshgrid(1:params.video.height,1:params.video.width,vidTimes);
vidMat = double(Y - min(Y(:)));
vidMat = interp3(xx,yy,zz,vidMat,xxx,yyy,zzz);
