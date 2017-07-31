%% Preprocess pipeline
% Folders on laptop:
laptop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs';
laptop_wav_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav';
laptop_gif_folder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/gif';
laptop_storage_folder = '/Volumes/CanaryData/DATA/lrb853_15/movs/';
DamagedFolder = '/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/too_large_or_damaged/';

% Folders on Data desktop:
desktop_mov_folder = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs';
desktop_storage_folder = '/Volumes/home/Data/Imaging/lrb853_15/RawData';

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
Prepare_Raw_Video_Audio;

%% 6. Annotate spectrograms 
