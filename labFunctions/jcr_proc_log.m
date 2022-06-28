% jcr_proc_log.m
%
% Processing of Motor Skill Reaching Task.
%
% Jeremy Cohen 2017-
%
%% General Path and Directory Setup
%%%%%%%
% SETUP
%%%%%%%
% Must always begin with setting up paths and adding default folders
close all
clear all
addpath(genpath('D:\matlab'))
run jcr_MLstart.m


%% JAABA pre-processing Intro
%
% This process uses 'JAABA' software written by Kristin Branson and
% colleagues.
% https://github.com/kristinbranson/JAABA

% This specific version requires matlab version 2014a.

% Combine/Stitch video files using cluster
% 1) Place folder containing video files on rootdir, must be 'Y:\from_tier2\Jay\videos'
vid_dir = [vid_rootdir 'jcr71\20200115'];
setUpDirUsingCluster(vid_dir, 'frontside', true);
% 2) find the output .txt file and copy it to vid_dir
% for JCR, 'D:\hantmanlab\reaching\data\results'
% 3) using Putty, login to the cluster "login1", enter userID/pwd
% change the filename and location in the line below and run it from the putty console.
% /groups/branson/bransonlab/share/qsub_genAllFeatures_v2.pl /groups/hantman/hantmanlab/from_tier2/Jay/videos/Jeremy/JAABADirsForCluster20170725T180130.txt 1 0
% 4) copy the entire video folder back to your own directory where you want
% to work with JAABA
%% JAABA - USE MATLAB 2014a
StartJAABA
% select 'open in training mode'
% select a .jab jaaba template file
% locate folder of videos associated with this .jab template file
% load files and classifier
% change experiment list
% delete all classifiers
% start again
%% JAABA post-processing - USE MATLAB 2014a
% sets variables/directories and then runs the function ethogram_plot()
%%%%%%%%%%%
% EXPT INFO
%%%%%%%%%%%
% requires 3 input strings that describe the experiment data folder names.
% for ex:
% expt_str_1 = 'JCR1';
% expt_str_2 = '_6-30-17';
% expt_str_3 = '20170630';
%
expt_str_1 = 'JCR5';
expt_str_2 = '_7-19-17';
expt_str_3 = '20170719';
%
% expt_str_1 = 'JCR3';
% expt_str_2 = '_5-1-17';
% expt_str_3 = '20170501';
%
range=10;
group=4; % group of trials to loop through
copyfiles=1; % 0==do not copy files, 1==delete/copy new files
Nframes=2500; % number of frames to plot in ethogram
figpos = [1930 320 1000 800]; % ethogram figure position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sets variables/directories and runs the function ethogram_plot().
JAABA_postprocessing_jcr(expt_str_1,expt_str_2,expt_str_3,range,group,copyfiles,Nframes,figpos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% APT - Animal Part Tracker - create movie_list.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) NoMachine: connect to janelia cluster - login1.int.janelia.org
% open matlab using the following command to run a 48hr session:
%  bsub -W 48:00 -XF -n8 /misc/local/matlab-2020b/bin/matlab -desktop
% in matlab cd to APT directory
%  cd /groups/hantman/hantmanlab/Jeremy/software/APT
%  APT.setpath
%  lObj = Labeler

% 2) After labeling some videos and training the tracker file
% Create a movie_list.txt file for each camera and each folder of videos
% In the server terminal window run the following commands for each session:
% processed video list:

%jcr76
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200616 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200616/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200616 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200616/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200617 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200617/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200617 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200617/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200618 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200618/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200618 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200618/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200619 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200619/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200619 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200619/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200622 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200622/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200622 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr76/20200622/movie_list_side.txt

%jcr77
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200714 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200714/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200714 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200714/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200715 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200715/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200715 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200715/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200716 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200716/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200716 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200716/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200717 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200717/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200717 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr77/20200717/movie_list_side.txt

%jcr79
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200929 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200929/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200929 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200929/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200930 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200930/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200930 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr79/20200930/movie_list_side.txt

%jcr80
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201023 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201023/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201023 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201023/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201028 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201028/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201028 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201028/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201103 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201103/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201103 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201103/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201104 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201104/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201104 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201104/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201105 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201105/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201105 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201105/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201110 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201110/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201110 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201110/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201111 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201111/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201111 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr80/20201111/movie_list_side.txt

% jcr82
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20201231 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20201231/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20201231 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20201231/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210101 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210101/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210101 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210101/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210104 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210104/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210104 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210104/movie_list_side.txt

find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210105 -name *front*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210105/movie_list_front.txt
find /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210105 -name *side*.avi > /groups/hantman/hantmanlab/Jeremy/data/video/bias/jcr82/20210105/movie_list_side.txt


%% APT - Create Dirs2track.json
% 3) CREATE a .json file containing a list of all directories to track in APT
% create the output directory folder for the .json file
% mkdir('Y:\data\video\dirs2track'); % run once

% RUN the following section of code in the linux server matlab window
% directories to track
dirs2track = {[vid_rootdir 'jcr76/20200616'] [vid_rootdir 'jcr76/20200617'] [vid_rootdir 'jcr76/20200618'] [vid_rootdir 'jcr76/20200619'] [vid_rootdir 'jcr76/20200622'] ...
    [vid_rootdir 'jcr77/20200714'] [vid_rootdir 'jcr77/20200715'] [vid_rootdir 'jcr77/20200716'] [vid_rootdir 'jcr77/20200717'] ...
    [vid_rootdir 'jcr79/20200929'] [vid_rootdir 'jcr79/20200930'] ...
    [vid_rootdir 'jcr80/20201023'] [vid_rootdir 'jcr80/20201028'] [vid_rootdir 'jcr80/20201103'] [vid_rootdir 'jcr80/20201104'] [vid_rootdir 'jcr80/20201105'] [vid_rootdir 'jcr80/20201110'] [rootdir 'jcr80/20201111'] ...
    [vid_rootdir 'jcr82/20201231'] [vid_rootdir 'jcr82/20210104'] [vid_rootdir 'jcr82/20210105']};

dirs2track = {[vid_rootdir 'jcr82/20210101']};

% location of output .json file for tracking in the APT
temp_dir = '/groups/hantman/hantmanlab/Jeremy/data/video/dirs2track';
% create the .json file
make_jsons_front_side_ALL(dirs2track,temp_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Track Videos in APT
% 4)TRACK videos in dirs2track
% Open APT, load the project for first camera, select Track -> Track multiple videos.
% A window will pop up; click Load, and select the *.json for the first camera, and click Track.
% This will generate one job for the GPU cluster per movie. Wait until all the jobs have been scheduled (this can take a long time - be patient!);
% Can check on the status by using “bjobs” in the terminal, and can count the number of jobs running using
% wait 15min or longer, 30min maybe.
% server terminal command
bjobs | wc -l
% Once all the jobs are in the queue, you can close the GUI, restart it, load the project for the second camera, and repeat.
% The tracker output for each individual movie will be a *.trk file in the data directory for each session, in the folder “tracked.”

 
%% Spike Sorting - Kilosort

% script to run Kilosort Autosorter
% for 1x64 probe - H3 from cambridge neurotech
edit D:\matlab\Kilosort-2.5\jeremy_chan_maps\runKS_1x64.m
% edit D:\matlab\Kilosort-2.5\jeremy_chan_maps\runKS_4x16.m



%% Get Event Timestamps - Create Trial Structure - Create Trial Outcome 
% Extract experimental timestamps
% 1) extract pulse timestamps from recording file (eg .nidq or .imec)
%%% camera triggers, lights, lasers, tables, tones, etc...
% 2) create variables that describe the experiment trial structure
%%% control trials, light trials, light only trials, etc...
% 3) create variables that describe trial by trial outcomes
%%% use/exclude trial flag, trial time, success/fail/, lift, reach,
%%% reaction time, reach velocity,

% local data directory for Wasabi cloud storage
wasabidir = 'E:\jeremycohen\thalamus_opto_reach\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment directory list
% video directories
vid_data_dirs = {'D:\localData\jcr70\video\20191016' 'D:\localData\jcr76\video\20200616' 'D:\localData\jcr76\video\20200617' ...
    'D:\localData\jcr77\video\20200714' 'D:\localData\jcr77\video\20200715' ...
    'D:\localData\jcr79\video\20200929' 'D:\localData\jcr79\video\20200930' ...
    'D:\localData\jcr79\video\20200930' 'D:\localData\jcr80\video\20201103' ...
    'D:\localData\jcr80\video\20201104' 'D:\localData\jcr80\video\20201110' ...
    'D:\localData\jcr80\video\20201111' 'D:\localData\jcr82\video\20201231' ...
    'D:\localData\jcr82\video\20210101' 'D:\localData\jcr82\video\20210104' ...
    'D:\localData\jcr82\video\20210105' 'D:\localData\M338\video\20210705' ...
    'D:\localData\M338\video\20210707' 'D:\localData\M340\video\20210811' ...
    'D:\localData\M341\video\20210817' 'D:\localData\M341\video\20210818' ...
    'D:\localData\jcr85_2\video\20200118Laser' 'D:\localData\jcr85\video\20200120Laser'};

% corresponding .nidq event ind directories
neural_data_dirs = {'D:\localData\jcr70\ephys\jcr70_20191016_920um_g0' 'D:\localData\jcr76\ephys\jcr76_20200616_4800_900um_g0' 'D:\localData\jcr76\ephys\jcr76_20200617_4800_925um_g0' ...
    'D:\localData\jcr77\ephys\jcr77_20200714_4600_900um_g0' 'D:\localData\jcr77\ephys\jcr77_20200715_4600_915um_g0' ...
    'D:\localData\jcr79\ephys\jcr79_20200929_4600_915um_g0' 'D:\localData\jcr79\ephys\jcr79_20200930_4500_915um_g0' ...
    'D:\localData\jcr79\ephys\jcr79_20200930_4500_915um_2_g0' 'D:\localData\jcr80\ephys\jcr80_20201103_950um_g0' ...
    'D:\localData\jcr80\ephys\jcr80_20201104_920um_g0' 'D:\localData\jcr80\ephys\jcr80_20201110_4400um_g0' ...
    'D:\localData\jcr80\ephys\jcr80_20201111_4400um_g0' 'D:\localData\jcr82\ephys\jcr82_20201231_4400um_g0' ...
    'D:\localData\jcr82\ephys\jcr82_20210101_4400um_g0' 'D:\localData\jcr82\ephys\jcr82_20210104_4455um_g0' ...
    'D:\localData\jcr82\ephys\jcr82_20210105_4600um_g0' 'D:\localData\M338\ephys\M338_20210705_1000_g0_imec0' ...
    'D:\localData\M338\ephys\M338_20210707_1000ab_g0_t0' 'D:\localData\M340\ephys\M340_20210811_1000_g0' ...
    'D:\localData\M341\ephys\M341_20210817_1000_g0_imec0' 'D:\localData\M341\ephys\M341_20210818_1000_g0_t0' ...
    'D:\localData\jcr85_2\ephys\JCR85_20210118_Cx950_g0_t0' 'D:\localData\jcr85\ephys\JCR85_20210120_Cx950_g0_t0'};

% subfolder location of the nidq meta file, which states the number of channels
site = {'cortex' 'cortex' 'cortex' ...
    'cortex' 'cortex' ...
    'cortex' 'cortex' ...
    'cortex' 'cortex' ...
    'cortex' 'thalamus' ...
    'thalamus' 'thalamus' ...
    'thalamus' 'thalamus' ...
    'thalamus' 'cortex' ...
    'cortex' 'cortex' ...
    'cortex' 'cortex' ...
    'cortex' 'cortex'};

% experiment type (eg, jcr or jay)
% 1==jcr jeremy rig, 2==jay janelia new rig, 3==jay janelia old rig
exp_type = [ones(1,16) 2*ones(1,5) 3*ones(1,2)];

% dir_list = [1:2 5 7 8:numel(neural_data_dirs)]; % list of files to process


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SESSION TIMESTAMPS
%
% JCR Experiments, Jeremy Cohen
% timestamps on Jeremy's janelia rig in .nidq file
% Get data from .nidq file channel #
% for 'exp_type==1'
ch1=64;% sync pulse timestamps, jc_pulse 1hz
ch2=66;% camera timestamps
ch3=68;% table_start/cue pulse timestamps
ch4=67;% laser pulse timestamps
ch5=70;% tone pulse timestamps
ch6=69;% masking led light pulse timestamps
jcr_pulse_channels = [ch1 ch2 ch3 ch4 ch5 ch6];

% Jay Guo Experiments
% timestamps on Jay's janelia corner rig in .nidq file
% Get data from .nidq file channel #
% for 'exp_type==2'
ch1=0;% sync pulse timestamps, imec 1hz
ch2=1;% camera timestamps
ch3=2;% table_start/cue pulse timestamps
ch4=3;% laser pulse timestamps
ch5=4;% laser gate, for grab detector trials
jay_pulse_channels_1 = [ch1 ch2 ch3 ch4 ch5];

% for 'exp_type==3'
ch1=64;% camera timestamps
ch2=65;% table_start/cue pulse timestamps
ch3=67;% laser pulse timestamps
ch4=66;% laser gate, for grab detector trials
jay_pulse_channels_2 = [ch1 ch2 ch3 ch4];

% experiment type (eg, jcr or jay)
% 1==jcr, 2==jay
exp_type = [ones(1,15) 2*ones(1,5) 3];

smoothwin = 25; % smoothing window (for plotting purposes)
tic;
for i=[1]%dir_list(5:end)
    
    % set current data directory - wasabi cloud storage - local mount
    wasabi_neural = [wasabidir neural_data_dirs{i}(14:end)];
    wasabi_vid = [wasabidir vid_data_dirs{i}(14:end)];
    wasabi_calib = 'E:\jeremycohen\camera_calibration\jeremy\20191030\Calib_Results_stereo.mat';
    
    % set rootdir, neural data, and video data directories
    local_or_wasabi=0;
    if local_or_wasabi==0
        rootdir='D:\localData\'
        neural_dir = neural_data_dirs{i}
        vid_dir = vid_data_dirs{i}
    elseif local_or_wasabi==1
        rootdir=wasabidir
        neural_dir = wasabi_neural
        vid_dir = wasabi_vid
    end
    save_dir=neural_dir;
    
    % stereo calibration file for video analysis
    if i<=15
        % jeremy rig - janelia experiments
        calib_file = 'D:\localData\camera_calibration\jeremy\20191030\Calib_Results_stereo.mat';
        wasabi_calib = 'E:\jeremycohen\camera_calibration\jeremy\20191030\Calib_Results_stereo.mat';
        calib_dir = calib_file;
    elseif i>=16 && i<=22
        % jay guo corner rig - janelia experiments
        calib_file = 'D:\localData\camera_calibration\jay\camera_calibration_corner_5-22-19\Calib_Results_stereo.mat';
        wasabi_calib = 'E:\jeremycohen\camera_calibration\jay\camera_calibration_corner_5-22-19\Calib_Results_stereo.mat';
        calib_dir = calib_file;
    end
    
    % set channels for event processing
    if exp_type(i)==1
        channels=jcr_pulse_channels;
    elseif exp_type(i)==2
        channels=jay_pulse_channels_1;
    elseif exp_type(i)==3
        channels=jay_pulse_channels_2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) Extract the timestamps for every channel in .nidq file
    % stored in the local directory in a file '/ind.mat'
%     get_event_ind_jcr(neural_dir,site{i},channels,exp_type(i),save_dir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) Set up trial array and categorize trials into trial type
    % append new variables to the event ind.mat file
%     create_trial_structure_jcr(rootdir,neural_dir,vid_dir,save_dir,...
%         site{i},exp_type(i))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) Create trial outcome variables, lift, reach, grasp, handopen, 
    % supinate, at mouth, success or fail outcome, trajectories.
    if local_or_wasabi==1
       pause(10); % wait for wasabi to finish uploading file to cloud
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Experiment specific info - DO NOT CHANGE %
    % TO ALIGN LASER+CUE AND CONTROL TRIALS, align to pellet shift onset
    if exp_type(i)==1
        if i < 9
            % FOR JCR80_20201103 and earlier experiments
            frames = 1:2000; % frame range (2ms frame rate)
            LT_shift = 250; % frame shift LT trials
            CON_shift = 25;
        else
            % FOR JCR80_20201104 and later experiments
            frames = 1:2000; % frame range (2ms frame rate)
            LT_shift = 25; % frame shift LT trials
            CON_shift = 25;
        end
    elseif exp_type(i)==2 || exp_type(i)==3
        frames = 1:2000; % frame range (2ms frame rate)
        LT_shift = 0; % frame shift LT trials
        CON_shift = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotfig=0; % plot trial trajectories during processing
%     create_trial_outcome_jcr(rootdir,neural_dir,vid_dir,save_dir,...
%         calib_dir,exp_type(i),smoothwin,frames,LT_shift,CON_shift,plotfig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) create the final data structure with all processed data
    % neural data and kinematics combined.
    new_or_append = 0; %==1 update/append to current file, ==0 create new file
    if new_or_append==0
        event_ind_dir = neural_data_dirs{i};
        ctx=1; % 1==process data, 0==do not process
        thal=0; % 1==process data, 0==do not process
        create_processed_data_struct(event_ind_dir,ctx,thal)
    elseif new_or_append==1
        update_processed_data_struct(event_ind_dir)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
b=toc;
disp(b);


%% PLOT - Individual Trial Position Data
t_win=1:4000; % time window in ms
smoothwin = 25; % 25 frames, 50ms
event_dir = [neural_data_dirs{i} filesep 'event_ind']; % event_ind.m location
trial=1:5;
plot_trial_jcr(event_dir,trial,t_win,smoothwin,exp_type(i));


%% Pre-process neural data and align to behavior - OLD
% Now we create a data_structure for each experiment that combines
% event_ind.mat and the spikes from each recording site
% this is step 4) in the pre-processing section above. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each experiment, there may be a 'cortex' and/or 'thalamus' recording
% with specific electrodes and recording configurations, so we process each
% recording site separately.
% Set input variables, and set the flag for the variables cortex and thalamus
% to (1==process, 0==dont process)

% jcr76_20200616 - cortex and thalamus
savedir = 'D:\localResults\thalamus_opto_reach\dualthalctx_ChR2\jcr76_20200616_preprocessed';
dir_vid = 'D:\localData\jcr76\video\20200616\tracked';
dir_ephys = 'D:\localData\jcr76\ephys\jcr76_20200616_4800_900um_g0';
mouse_type = 'rorb';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr76_20200617 - cortex and thalamus
savedir = 'D:\localResults\thalamus_opto_reach\dualthalctx_ChR2\jcr76_20200617_preprocessed';
dir_vid = 'D:\localData\jcr76\video\20200617\tracked';
dir_ephys = 'D:\localData\jcr76\ephys\jcr76_20200617_4800_925um_g0';
mouse_type = 'rorb';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr79_20200929 - cortex and thalamus
savedir = 'D:\localResults\thalamus_opto_reach\dualthalctx_ChR2\jcr79_20200929_preprocessed';
dir_vid = 'D:\localData\jcr79\video\20200929\tracked';
dir_ephys = 'D:\localData\jcr79\ephys\jcr79_20200929_4600_915um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);


% jcr80_20201110 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr80_20201110_preprocessed';
dir_vid = 'D:\localData\jcr80\video\20201110\tracked';
dir_ephys = 'D:\localData\jcr80\ephys\jcr80_20201110_4400um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 0; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr80_20201111 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr80_20201111_preprocessed';
dir_vid = 'D:\localData\jcr80\video\20201111\tracked';
dir_ephys = 'D:\localData\jcr80\ephys\jcr80_20201111_4400um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 0; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);


% %%% GET THIS DATA SET % jcr82_20201222 - cortex
% savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20201222_preprocessed';
% dir_vid = 'D:\localData\jcr82\video\20201222\tracked';
% dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20201222_XXX_g0';
% mouse_type = 'necab';
% calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
% n_frames = 2200;
% fs_vid = 500;
% fs_ephys = 25000;
% cue_offset_c = 950; % cue offset for control trials
% cue_offset_l = 500; % cue offset for laser trials
% cortex = 1; % 1==cortex data
% thalamus = 1; % 1==thalamus data
% jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
%                         n_frames,mouse_type,...
%                         fs_vid,fs_ephys,...
%                         cue_offset_c,cue_offset_l,cortex,thalamus);
%
% %%% GET THIS DATA SET % jcr82_20201223 - cortex
% savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20201223_preprocessed';
% dir_vid = 'D:\localData\jcr82\video\20201223\tracked';
% dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20201223_XXX_g0';
% mouse_type = 'necab';
% calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
% n_frames = 2200;
% fs_vid = 500;
% fs_ephys = 25000;
% cue_offset_c = 950; % cue offset for control trials
% cue_offset_l = 500; % cue offset for laser trials
% cortex = 1; % 1==cortex data
% thalamus = 1; % 1==thalamus data
% jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
%                         n_frames,mouse_type,...
%                         fs_vid,fs_ephys,...
%                         cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr82_20201231 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20201231_preprocessed';
dir_vid = 'D:\localData\jcr82\video\20201231\tracked';
dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20201231_4400um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr82_20210101 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20210101_preprocessed';
dir_vid = 'D:\localData\jcr82\video\20210101\tracked';
dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20210101_4400um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr82_20210104 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20210104_preprocessed';
dir_vid = 'D:\localData\jcr82\video\20210104\tracked';
dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20210104_4455um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);

% jcr82_20210105 - thalamus
savedir = 'D:\localResults\thalamus_opto_reach\jcr82_20210105_preprocessed';
dir_vid = 'D:\localData\jcr82\video\20210105\tracked';
dir_ephys = 'D:\localData\jcr82\ephys\jcr82_20210105_4600um_g0';
mouse_type = 'necab';
calib_file = 'D:\localData\camera_calibration\20191030\Calib_Results_stereo';
n_frames = 2200;
fs_vid = 500;
fs_ephys = 25000;
cue_offset_c = 950; % cue offset for control trials
cue_offset_l = 500; % cue offset for laser trials
cortex = 1; % 1==cortex data
thalamus = 1; % 1==thalamus data
jcr_preprocess_ephys(dir_vid,dir_ephys,calib_file,savedir,...
    n_frames,mouse_type,...
    fs_vid,fs_ephys,...
    cue_offset_c,cue_offset_l,cortex,thalamus);



















