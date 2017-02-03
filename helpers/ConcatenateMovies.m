function [vidMat, vidTimes, Aud] = ConcatenateMovies(mov_files_list)
Aud = struct('nrChannels',1,'bits',16,'nrFrames',0,'data',[],'rate',48000,'TotalDuration',0);
vidMat = [];
vidTimes = [];
for fnum = 1:numel(mov_files_list) 
    params = load(mov_files_list{fnum});
    Aud.nrFrames = Aud.nrFrames + params.Aud.nrFrames;
    Aud.data = [Aud.data; params.Aud.data];
    if isfield(params.Aud,'TotalDurration')
        Aud.TotalDuration = Aud.TotalDuration + params.Aud.TotalDurration;
    end
    if isfield(params.Aud,'TotalDuration')
        Aud.TotalDuration = Aud.TotalDuration + params.Aud.TotalDuration;
    end
    vidMat = cat(3,vidMat,params.vidMat);
    if fnum == 1
        vidTimes = params.vidTimes;
    else
        vidTimes = [vidTimes; 2*vidTimes(end) - vidTimes(end-1) + params.vidTimes];
    end
end
