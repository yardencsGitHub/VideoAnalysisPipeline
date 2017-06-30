%%This script takes a list of extracted wav files and creates gif and mat
%%spectrograms under the folders 'gif' and 'mat'
DIR = pwd;
%% repositories
% addpath(genpath('/Users/tonatiuh/Documents/GitHub'));
addpath(genpath('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub');
%% confirm
wavs = dir('*.mov');
for i=1:numel(wavs)
    nm = wavs(i).name;
    sp = regexp(nm(1:end-4),'_', 'split');
    if numel(sp) ~= 8
        disp(nm);
    end
    if ~strcmp(sp{1},bird_name)
        disp(nm);
    end
    if (max(sp{2}*1)>57 || min(sp{3}*1)<48)
        disp(nm);
    end
    if ~strcmp(sp{3},'2017')  
      disp(nm);
    end
    if (max([sp{4:8}]*1)>57 || min([sp{4:8}]*1)<48)
        disp(nm);
    end
    if numel([sp{4:8}]) ~= 10
        disp(nm);
    end
end
%% 2. ExtractWavsFromMOVs
cd(DIR);
ExtractWavsFromMOVs;
disp('done extracting wavs');
%% 3. Create spectrogram image data
% create spectrograms
% pad with 1 sec of silence in each side
% save gif,mat, and csv
cd([DIR '/wav']);
DIR = pwd;
if ~exist(fullfile(DIR,'gif'),'dir')
	mkdir(fullfile(DIR,'gif'));
end
if ~exist(fullfile(DIR,'mat'),'dir')
	mkdir(fullfile(DIR,'mat'));
end
if ~exist(fullfile(DIR,'csv'),'dir')
	mkdir(fullfile(DIR,'csv'));
end
wavs = dir('*.wav');
clipping=[-2 2];
disp_band=[1 9e3];
colors='hot';
Nwavs = numel(wavs);
for file_cnt = 1:numel(wavs)
    
    wavname = wavs(file_cnt).name;
    display(wavname)
    outname = wavname(1:end-4);
    [signal,fs] = audioread(wavname);
    [s,f,t]=zftftb_pretty_sonogram(signal,fs,'len',16.7,'overlap',14,...
                'zeropad',0,'norm_amp',1,'clipping',clipping);
    
    startidx=max([find(f<=disp_band(1))]);
    stopidx=min([find(f>=disp_band(2))]);

    im=s(startidx:stopidx,:)*62;
    im=flipdim(im,1);

    imwrite(uint8(im),colormap([ colors '(63)']),fullfile(DIR,'gif',[ outname '.gif']),'gif');
    save(fullfile(DIR,'mat',[ outname '.mat']),'f','t','s');
    csvwrite(fullfile(DIR,'csv',[ outname '.txt']),s);
    disp([Nwavs file_cnt])
end
    
%%
%% 1. Rename all files - historical
% bird_name = 'lrb85315';
% wavs = dir('*.mov');
% for wav_num = 1:numel(wavs)
%     res_str = [ConvertFilenameFS2Tweet(wavs(wav_num).name,bird_name,wav_num) '.mov'];
%     [SUCCESS,MESSAGE,MESSAGEID] = movefile(wavs(wav_num).name,res_str);
% end
% disp('done renaming');
