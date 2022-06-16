%BehaviorTracking_SP
% Steps to do before running behavior script. Note, should run script step
% by step. Not all at once.
% Steps:
%   1) Train a tracker in APT.
%   2) Track all videos with created tracker.
%   3) Run JAABA to label parts of behavior. MATLAB 2014 for now.
runAPT = 0;
runJAABAPrep = 0;
runJAABAGui = 0;
runPostProcessing = 1;
%% Set Paths
BASEPATH = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\';
CODE_GENERAL = [BASEPATH 'Code\General\'];

SUB = 'jcr70';
EXPER_SESSION = '20191016';
MOVIES = [BASEPATH 'ThalamusSP_Project\Data\jcr70\video\20191016\'];
APT = ['C:\MATLAB\APT-develop\APT-develop'];

addpath(genpath(CODE_GENERAL));
%% Run APT steps
if runAPT
    addpath(APT);
    %% Train APT on front & side of videos
    % Load and train a file
    % use rcb_tracker_front as starter to build on
    lObj = StartAPT;

    %% Track Front Video with APT
    % Run the following
    Front_Obj = StartAPT;
    % 1) Select and load correct trainer .lbl for front videos
    % 2) Run the following
    trackVideosAPT(MOVIES,'front',Front_Obj);

    %% Track Side Video with APT
    % Run the following
    Side_Obj = StartAPT;
    % 1) Select and load correct trainer .lbl for side videos
    % 2) Run the following
    trackVideosAPT(MOVIES,'side',Side_Obj);

    %% Run trk files are in directory under trk folder
    %     TRK = [MOVIES 'trk\'];
    %     FRONT_TRK = [TRK 'front\'];
    %     SIDE_TRK = [TRK 'side\'];
    %     if ~exist(FRONT_TRK)
    %         mkdir(FRONT_TRK);
    %         front_filelist = dir(fullfile(TRK, 'tracked_front_*.trk')); %get list of files
    %         for ifile = 1:length(front_filelist)
    %             oldFileDir = [TRK front_filelist(ifile).name];
    %             newFileDir = [FRONT_TRK front_filelist(ifile).name];
    %             movefile(oldFileDir,newFileDir);
    %         end
    %     end
    %     if ~exist(SIDE_TRK)
    %         mkdir(SIDE_TRK);
    %         side_filelist = dir(fullfile(TRK, 'tracked_side_*.trk'));
    %         for ifile = 1:length(side_filelist)
    %             oldFileDir = [TRK side_filelist(ifile).name];
    %             newFileDir = [SIDE_TRK side_filelist(ifile).name];
    %             movefile(oldFileDir,newFileDir);
    %         end
    %     end
end
%% Run JAABA prep step
if runJAABAPrep
    JAABAST = [BASEPATH 'Code\Packages\JAABAST_2021'];
    GENERAL = [BASEPATH 'Code\General'];
    addpath(GENERAL);
    addpath(genpath(JAABAST));
    %% Run JAABA
    VideoCombiningPrep(MOVIES);
end
%% Run JAABA GUI
if runJAABAGui
    JAABAST = [BASEPATH 'Code\Packages\JAABAST'];
    addpath(genpath(JAABAST));
    StartJAABA;
end
if runPostProcessing
    JAABAST = [BASEPATH 'Code\Packages\JAABAST'];
    addpath(genpath(JAABAST));
    SESSION = [MOVIES 'JAABA_Combo\'];
    SESSION = {SESSION};
    %Another instance run this
    [red,jabfiles,T0,T1] = ethogram_plot(SESSION,[],1500,'automarks',true,'exptags','readfromjabs');
     % Export to workspace and save mat file
    % save([BASEPATH 'Data_Analyzed\JAABA_Output\' SUB '_' EXPER_SESSION '_JAABA_Output.mat'],'data','stats');
end