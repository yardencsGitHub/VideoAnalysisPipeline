DIR = pwd;
fs=48000;
Fst1 = 500;  % The max frequency (Hz) of the low stop band  
Fp1 = 650;   % The min frequency (Hz) of the pass band
Fp2 = 20000;  % The max frequency (Hz) of the pass band
Fst2 = 21000; % The min frequency (Hz) of the high stop band
Ast1 = 60;   % ratio of suppression in low stop band (dB)
Ap = 1;      % Pass band ripple ratio (dB)
Ast2 = 60;   % ratio of suppression in high stop band (dB)
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2,fs);
Hd = design(d,'equiripple');
cd('/Users/yardenc/Documents/Experiments/Code and Hardware Dev/GitHub/FreedomScope/Analysis Pipeline');
addpath(genpath(pwd));
cd(DIR);
FILES = dir('*.mov');
if ~exist(fullfile(pwd,'wav'),'dir')
    mkdir('wav');
end
for nmov = 1:numel(FILES)
    [a_ts, a, v_ts, v] = extractmedia(FILES(nmov).name);
    audiowrite([DIR '/wav/' FILES(nmov).name(1:end-3) 'wav'],filtfilt(Hd.Numerator,1,double(a)),fs);
end