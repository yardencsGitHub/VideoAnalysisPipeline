function [vidMat, vidTimes, Aud] = FS2MAT(movfile)
audio_field = 'audio';
FSfolder = '/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline';
addpath(FSfolder);
params = load(movfile);
if isfield(params,audio_field)
    Aud = params.audio;
end
[Y,n]=FS_Format(params.video.frames,1,1);
vidMat = Y - min(Y(:));
vidTimes = params.video.times;