%1 Choose FS movies
[FileName,PathName,FilterIndex] = uigetfile([pwd '/*.mov'],'MultiSelect','on');
%2 Convert those movies to Matlab structure 
FS_AV_Parse_batch(pwd,FileName);
% Convert files to intensity video matrix format
cd([pwd '/mat']);
matlist = {};
for fnum = 1:numel(FileName)
    fname = FileName{fnum}; fname(end-2:end) = 'mat';
    matlist{fnum} = fname;
    [vidMat, vidTimes, Aud] = FS2MAT(fname);
    save(fname,'vidMat','vidTimes','Aud','-v7.3');
end
% 3 remove first dark frames
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
[vidMat, vidTimes, Aud] = ConcatenateMovies({movlist.name});

%%
CNMFE_Folder = '/Users/yardenc/Documents/Experiments/Imaging/PracticeData/mat/reduced/cnmfe';
if exist(CNMFE_Folder) == 0
    mkdir(CNMFE_Folder);
end
factor = 0.5;
Y = vidMat;
Y = uint8(imresize(double(Y),factor));
Ysiz = size(Y)';
save(fullfile(CNMFE_Folder,'cnmfe_source.mat'),'Y','Ysiz','-v7.3');




