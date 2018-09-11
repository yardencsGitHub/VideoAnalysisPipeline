function create_baseline_subtracted_raw_data(birdnum,varargin)

bird1_params = {'lrb85315' 'lrb853_15' 'lrb85315template' 'lrb85315auto_annotation5_fix'};
bird2_params = {'lbr3022' 'lbr3022' 'lbr3022_template' 'lbr3022auto_annotation5_alexa'};
bird3_params = {'lbr3009' 'lbr3009' 'lbr3009_template_4TF' 'lbr3009_annotation_4TF'};
switch birdnum
    case 1
        bird_params = bird1_params;
    case 2
        bird_params = bird2_params;
    case 3
        bird_params = bird3_params;
end

bird_name = bird_params{1}; 
bird_folder_name = bird_params{2}; 
template_file = bird_params{3}; 
annotation_file = bird_params{4}; 

nparams=length(varargin);
temporal_smooth = 0;
for i=1:2:nparams
	switch lower(varargin{i})
		case 'temporal_smooth'
			temporal_smooth = varargin{i+1};
    end
end

SourceDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/RawData'];
TargetDir = ['/Volumes/home/Data/Imaging/' bird_folder_name '/BaselineSubtractedRawData'];

if ~exist(TargetDir,'dir')
    mkdir(TargetDir);
end

filt_rad = 150; filt_sigma = 145; % highpass filter 
h = fspecial('gaussian',filt_rad,filt_sigma);
d = dir(fullfile(SourceDir,['RawData_' bird_name '*.mat']));
for fnum = 1:numel(d)
    fname = d(fnum).name;
    load(fullfile(SourceDir,fname));
    Y = vidMat;
    base = imfilter(Y,h,'circular','replicate');
    if (temporal_smooth == 1)
        base = convn(base, reshape([1 1 1] / 3, 1, 1, []), 'same');
    end
    c = Y-base;
    c = c-min(c(:));
    %dffMat = bsxfun(@rdivide,bsxfun(@minus,c,quantile(c,0.05,3)),quantile(c,0.05,3));
    BGlessVidMat = c;
    new_fname = ['BGless_' fname];
    save(fullfile(TargetDir,new_fname),'BGlessVidMat','vidTimes','-v7.3');
    display(fname);
end


