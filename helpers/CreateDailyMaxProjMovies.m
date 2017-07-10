%% CreateDailyMaxProjMovies
% Takes a folder full of maximum projections (in format 'fmt') and creates movies of
% concatenated max proj. frames. It also creates the folders /movies, /images and in
% them the max(max(projections)) images and movies
% Each day goes to one movie. 
% Only data the the (Tweet) date label exceeding 'first_day' are processed
DIR = pwd;
fs = 30;
fmt = 'png';
first_day = datenum('2017-02-20');
FILES = dir(['*.' fmt]);
if ~exist('max_movies','dir')
    mkdir('max_movies');
end
if ~exist('max_images','dir')
    mkdir('max_images');
end
dates = [];
for cnt = 1:numel(FILES)
    dates = [dates; TweetDateNum(FILES(cnt).name)];
end
[dates,indx] = sort(dates);
startloc = min(find(dates >= first_day));
FILES = FILES(indx);

if ~isempty(startloc)
    days = unique(dates);
    for daynum = 1:numel(days)
        locs = find(dates == days(daynum));
        tokens = regexp(FILES(locs(1)).name,'_','split');
        mov_fname = ['MaxProjMov_' tokens{1} '_' datestr(days(daynum),'yyyy_mm_dd') '.mp4'];
        im_fname = ['MaxMaxProj_' tokens{1} '_' datestr(days(daynum),'yyyy_mm_dd') '.png'];
        stack = [];
        for cnt = 1:numel(locs)
            I = imread(FILES(locs(cnt)).name);
            stack = cat(3,stack,I);
            M(cnt) = im2frame(I,gray(256));
        end
        v = VideoWriter(fullfile(DIR,'max_movies',mov_fname),'MPEG-4');
        v.FrameRate = fs;
        open(v)
        writeVideo(v,M);
        close(v)
        clear M;
        I = max(stack,[],3);
        imwrite(I,gray(256),fullfile(DIR,'max_images',im_fname));
    end
    
    
end

function res = TweetDateNum(txt)
    tokens = regexp(txt,'_','split');
    res = datenum([tokens{3} '-' tokens{4} '-' tokens{5}]);
end