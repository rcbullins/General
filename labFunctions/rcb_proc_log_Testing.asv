% jcr_proc_log.m
%
% Processing of Motor Skill Reaching Task.
%
% Jeremy Cohen 2017-
%
%% General Path and Directory Setup
ProjectName = 'ThalamusSP_Project';
SUB = 'jcr70';
EXPER_SESSION = '20191016';

BASEPATH = 'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\';
PROJECT  = [BASEPATH ProjectName '\'];

CODE_SP = [PROJECT 'Code\'];
CODE_GENERAL = [BASEPATH 'Code\General\'];
CODE_CALIB = [PROJECT 'Data\Calib_Results_stereo.mat'];
CODE_BRITTON = [BASEPATH 'Code\Packages\britton_code\code\'];
CODE_BRITTON_PLOT = [BASEPATH 'Code\Packages\britton_code\other_code'];
CODE_PROBES = [BASEPATH 'Code\Packages\Probes\'];
JAABA_OUTPUT = [PROJECT 'Data_Analyzed\' SUB '\JAABA_Output\' SUB '_' EXPER_SESSION '_JAABA_Output.mat'];

addpath(genpath(CODE_GENERAL));
addpath(CODE_CALIB);
addpath(genpath(CODE_BRITTON));
addpath(genpath(CODE_BRITTON_PLOT));
addpath(genpath(CODE_PROBES));

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
vid_data_dirs = {'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\ThalamusSP_Project\jcr70\video\20191016'};

% corresponding .nidq event ind directories
neural_data_dirs = {'C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\ThalamusSP_Project\Data\jcr70\ephys\jcr70_20191016_920um_g0'};
% subfolder location of the nidq meta file, which states the number of channels
site = {'cortex'};

% experiment type (eg, jcr or jay)
% 1==jcr jeremy rig, 2==jay janelia new rig, 3==jay janelia old rig
exp_type = [1];

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
    
        rootdir='C:\Users\bullinsr\OneDrive - University of North Carolina at Chapel Hill\Hantman_Lab\ThalamusSP_Project\';
        neural_dir = neural_data_dirs{i};
        vid_dir = vid_data_dirs{i};
        
        save_dir = neural_dir;
    
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



















