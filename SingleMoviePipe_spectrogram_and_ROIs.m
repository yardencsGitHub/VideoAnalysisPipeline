function [vidMat, vidTimes, Aud, dff2] = SingleMoviePipe_spectrogram_and_ROIs(filename,startframe)
% parameters:
filt_rad = 10; filt_sigma = 10; baseline_per=5; disp_band = [100 9e3];
% Load a FreedomScope .mov file and arrange data into params.video and
% params.audio
    FSfolder = '/Users/yardenc/Documents/GitHub/FreedomScope/Analysis Pipeline';
    addpath(FSfolder);
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
     %(Aud.nrFrames/Aud.rate-1/30);
    [Y,n]=FS_Format(params.video.frames,1,1);
    [xx,yy,zz] = meshgrid(1:params.video.height,1:params.video.width,params.video.times);
    [xxx,yyy,zzz] = meshgrid(1:params.video.height,1:params.video.width,vidTimes);
    vidMat = double(Y - min(Y(:)));
    vidMat = interp3(xx,yy,zz,vidMat,xxx,yyy,zzz);
    clear params;
% Play movie to allow selecting frames to delete
    if (startframe < 1)
        movie = Mat2Mov(vidMat);
        mplay(movie);
        ndel_frames = input('how many frames to delete? :');
        clear movie;
    else
        ndel_frames = startframe;
    end
% Trim video and audio
    startidx = round(Aud.rate * vidTimes(ndel_frames + 1)); 
    Aud.nrFrames = Aud.nrFrames - startidx;
    Aud.data(1:startidx) = [];
    Aud.TotalDuration = Aud.TotalDuration - startidx / Aud.rate;
    audiowrite('temp.mp4',Aud.data,Aud.rate);
    vidMat(:,:,1:ndel_frames) = [];
    vidTimes(1:ndel_frames) = [];
    vidTimes = vidTimes - vidTimes(1);
% Smooth and downsample
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
   
    %figure; imshow(I);
    %caxis([0 0.2]);
% Manual ROI extraction
    [ROI]=FS_Image_ROI(uint8(I/max(I(:))*255)*1.5);
% Create df/f traces aligned to spectrogram
    [im,f,t] = zftftb_pretty_sonogram(Aud.data,Aud.rate,...
    'len',16.7,'overlap',14,'zeropad',0,'norm_amp',1,'clipping',[-2 2]);
    startidx = max([find(f <= disp_band(1))]);
    stopidx = min([find(f >= disp_band(2))]);
    im = im(startidx:stopidx,:)*62;
    im = flipdim(im,1);
    figure;
    subplot(1+length(ROI.stats),1,1);
    imshow(uint8(im),colormap(['hot(63)']));
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
       
    


