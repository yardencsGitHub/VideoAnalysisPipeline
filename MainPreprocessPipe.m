%% Preprocess pipeline
% Folders on laptop:
laptop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs';
laptop_wav_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
laptop_gif_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/gif';
laptop_storage_folder = '/Volumes/CanaryData/DATA/lrb853_15/movs/';
laptop_annotated_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated';
laptop_annotated_images_dir = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images';
DamagedFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/too_large_or_damaged/';

% Folders on Data desktop:
desktop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs';
desktop_storage_folder = '/Volumes/home/Data/Imaging/lrb853_15/RawData';
desktop_max_projections_dir = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs/MaxProj';

% Parameters
last_idx = 7982;
last_date = '2017_06_29';
last_time = '08_24_46';
bird_name = 'lrb85315';


%% 1. Convert FS file names to canonical Tweet file names and save them (on work laptop)
% Make sure movies storage is mapped correctly
% Run script "Move_source_FS_files_to_Tweet_folder" after setting the correct file numbers and bird name inside that script.
SourceFolder = laptop_storage_folder; %'/Volumes/CanaryData/DATA/lrb853_15/movs/';
TargetFolder = laptop_mov_folder; %'/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/';
Move_source_FS_files_to_Tweet_folder;
%% 2. Extract audio and save wav files.
% On work laptop, cd to the movies folder.
% Run script "ExtractWavsFromMOVs" after setting the correct parameters inside.
cd(laptop_mov_folder);
startfrom = last_idx+1;
ExtractWavsFromMOVs;
%% 3. Create spectrograms 
% Create a list of wav files to process by running the scipt "CreateWavsList" in the folder with all wav files from (2)
% Remember to set all parameters in that script correctly
% Run the script "CreateSpectrogramsFromWavs" in the wav files folder.
% Remember to set parameterts in this script.
%
cd(laptop_wav_folder);
startfrom = last_idx+1;
CreateWavsList;
CreateSpectrogramsFromWavs;
%% 4. Prune out noise files
% use script "locate_noise_images" to go over the images, created in (3), and mark noise and
% song fils. Create the list of movies to process and copy it to the
% desktop computer
cd(laptop_gif_folder);
startfrom = last_idx+1;
locate_noise_images;
cd(laptop_mov_folder);
save FS_movies_list keys songfiles noisefiles;
%% 5. Now copy all movies and spectrograms to the Data desktop computer
% This should be done manually to allow for a controlled checkpoint
% The following parts run on the desktop computer (The data computer)


%% 5. Create Raw data variables and store them on the storage device
% Make sure the storage device is mapped
% Move list of movie files to process from the laptop (created in (4))
% Run script "Prepare_Raw_Video_Audio"
cd(desktop_mov_folder);
OutputFolder = desktop_storage_folder; % '/Volumes/home/Data/Imaging/lrb853_15/RawData';
Prepare_Raw_Video_Audio;

%% 6. Annotate spectrograms (in parallel to 5)
% Copy the spectrogram matlab files from (3) to the desktop 
% Set the folder names, parameters, and training checkpoint in the
% annotation script and run.
% Save results and copy to laptop
% On laptop
% Create annotation file or add newly annotated results to old one
% (create_annotation_from_auto_addition) 
% Set parameters accordingly
% Copy the file FS_movies_list.mat, created in(4) to the laptop_wav_folder.
% Needs fixing: create_annotation_from_auto_addition 
DIR = laptop_wav_folder; %'/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
annDIR = laptop_annotated_dir;
auto_file = 'Results_Aug_01_2017.mat';
old_annotation_file = 'lrb85315auto_annotation4';
new_annotation_file = 'lrb85315auto_annotation5';
template_file = 'lrb85315template.mat';
corrections = 0; % a flag that indicates that we're going to reopen old annotations and replace some syllables' annotation
syllables_to_reannotate = [300]; %  Theses syllables in the old files will trigger reannotation
trill_syllables = [0:2 4 5 8 9 200 203 208 209 300:306 308 309];
copyfile(fullfile(laptop_mov_folder,'FS_movies_list.mat'),fullfile(laptop_wav_folder,'FS_movies_list.mat'));
create_annotation_from_auto_addition;
keys = {params.keys{:} tempkeys{:}};
elements = [params.elements elements];
save(new_annotation_file,'elements','keys');
%% 7. Clean annotation results (on laptop)
% Move new annotation file to laptop_wav_folder
% Create phrase and spectrogram images using the script
% "show_spec_and_phrases"
% Use TweetVisionLite to fix mistakes
targetdir = laptop_annotated_images_dir; %'/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images';
template_file = 'lrb85315template.mat';
cd(laptop_wav_folder);
load(template_file);
syllables = [[templates.wavs.segType] -1 102 103];
show_spec_and_phrases;
%% 8. Create Maximum projection images and movies (Desktop)
% Create the maximum, background subtracted, projection for each song movie
% (using the scipt CreateMaxProjections)
% Then create daily max. proj movies using "CreateDailyMaxProjMovies"
% Set parameters accordingly
startfrom = last_idx+1;
OutputDIR = desktop_max_projections_dir;
DIR = desktop_storage_folder;
CreateMaxProjections;
cd(desktop_max_projections_dir);
first_day = datenum('2017-06-30');
CreateDailyMaxProjMovies;
%% 9. Create manual ROI daily alignments (Desktop)
% This is done by running the script "SingleDayManualROIsDataExtraction"



