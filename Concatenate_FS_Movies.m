function [vidMat, audio] = Concatenate_FS_Movies(DataFolder,WorkingFolder,factor,outname)
%%
audio_field = 'audio';
FSfolder = '/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline';
% DataFolder = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/';
% WorkingFolder = '/Users/yardenc/Documents/GitHub/CNMF-E/demos';
addpath(FSfolder);
%%
% This script will concatenate FreedomScope movies in .mov and .mat forms
% Inputs:
%   DataFolder    - The full path to the movie files. Default = pwd.
%   WorkingFolder - The full path for the output file. Default = pwd. If
%                   set to 'none' then no file will be saved.
%   factor        - Scaling factor for the video frames (for memory issues)
%   outname       - Name of output file.
% Outputs:
%   vidMat - The 3d movies [width]x[height]x[frames] in uint8 format
%   audio  - The concatenated audio, if exists
% Dependencies:
% Requires the function FS_Format from ..[GitHubRoot]/FreedomScope/Analysis Pipeline
%%
if isempty(DataFolder)
    DataFolder = pwd;
end
if isempty(WorkingFolder)
    WorkingFolder = pwd;
end
if isempty(factor)
    factor = 0.5;
end
if isempty(outname)
    outname = 'Ca_temp';
end

[FileName,PathName,FilterIndex] = uigetfile(fullfile(DataFolder,'*.mat;*.mov'),'MultiSelect','on');
Ytemp = [];
Aud = struct('nrChannels',1,'bits',16,'nrFrames',0,'data',[],'rate',48000,'TotalDuration',0);
for fnum = 1:numel(FileName)
    [pahstr,nm,ext] = fileparts(FileName{fnum});
    if strcmp(ext,'.mov')
        FS_AV_Parse_batch(PathName,FileName);    
        PathName = [PathName '/mat']; 
        FileName = cellfun(@(a)[a(1:end-3) 'mat'],FileName,'UniformOutput',false);
    end
    
    params = load(fullfile(PathName,FileName{fnum}));
    if isfield(params,audio_field)
        Aud.nrFrames = Aud.nrFrames + params.audio.nrFrames;
        Aud.data = [Aud.data; params.audio.data];
        if isfield(params.audio,'TotalDurration')
            Aud.TotalDuration = Aud.TotalDuration + params.audio.TotalDurration;
        end
        if isfield(params.audio,'TotalDuration')
            Aud.TotalDuration = Aud.TotalDuration + params.audio.TotalDuration;
        end
    end
    if isfield(params,'video') 
        [Y,n]=FS_Format(params.video.frames,1);
    end
    if isfield(params,'Y')
        Y = params.Y;
    end
    Y = imresize(Y,factor);
    Y = Y - min(Y(:));
    Ytemp = cat(3,Ytemp,Y);
end
Y = Ytemp;
Ysiz = size(Y)';
if (~strcmp(WorkingFolder,'none'))
    target = fullfile(WorkingFolder,outname);
    save(target,'Y','Ysiz','-v7.3');
end
vidMat = Y;
audio = Aud;