function create_processed_data_struct(event_ind_dir,ctx,thal)

% script has two sections:
% 1)process event_indices
% 2) process the neural data

% find the meta file with the event indices.
% whisper system sampled timestamps occur at 25 khz, but should check
% timestamp meta file.
samp_freq = 25; % sampling frequency for the timestamps in the event_indices file.


%% Process event indices
load([event_ind_dir filesep 'event_ind.mat'])

neural_dir_ctx = [event_ind_dir filesep 'cortex'];
neural_dir_thal = [event_ind_dir filesep 'thalamus'];

% experiment name for saving
mname = ls([neural_dir_ctx '\*.nidq.meta']);
fname=mname(1:end-10);
save_file = [event_ind_dir filesep 'data_struct_' fname];

% create data struct
data_struct=struct;

%label names
track_names={'digit_2','digit_4','wrist','pellet'}; %name used for struct
track_names_v2={'dig2','dig4','wrist','pellet'};    %name used in event indices file

%other relevant data for struct
time=-4999:10000;
nTrials=size(trial_array,1);
nTPs=size(trajectory_struct.dig2_x,2);

%loop across tracked points and trials adding each to the data struct
for pt=1:4
    for tr=1:nTrials
        x_intrp=trajectory_struct.(strcat(track_names_v2{pt},'_x'))(tr,:);
        y_intrp=trajectory_struct.(strcat(track_names_v2{pt},'_y'))(tr,:);
        z_intrp=trajectory_struct.(strcat(track_names_v2{pt},'_z'))(tr,:);
        
        %pad the kinematic trajectories with nan
        x_intrp=[nan([1,5000]),x_intrp,nan([1,10000-nTPs])];
        y_intrp=[nan([1,5000]),y_intrp,nan([1,10000-nTPs])];
        z_intrp=[nan([1,5000]),z_intrp,nan([1,10000-nTPs])];
        data_struct(tr).(strcat(track_names{pt},'_x'))=x_intrp;
        data_struct(tr).(strcat(track_names{pt},'_y'))=y_intrp;
        data_struct(tr).(strcat(track_names{pt},'_z'))=z_intrp;
        data_struct(tr).time=time;
    end
end

%massaging some of the event indices for use below
tmp_LT_trial=zeros(nTrials,1);
tmp_LT_trial=[tmp_LT_trial,nan(nTrials,2)];
tmp_LT_trial([LT_trial;LT_only],1)=1;
trn_trials=setdiff(1:nTrials,LT_only);
tmp_table=zeros(nTrials,1);
tmp_table(trn_trials)=table_start;


% loop across trials defining the start of the trial, table and laser cue
for tr=1:nTrials
    
    data_struct(tr).trial_start=[1,trial_start(tr)/samp_freq]; %relative time in trial, global time relative to recording
    data_struct(tr).table=[tmp_table(tr)/samp_freq,tmp_table(tr)/samp_freq-trial_start(tr)/samp_freq]; %relative time, global time
    data_struct(tr).laser_trial=tmp_LT_trial(tr,:);
    
    %get laser pulses
    if tmp_LT_trial(tr,1)==1
        tmp_laser_pulses=laser./samp_freq;
        
        %collect up all laser pulses in the current trial
        tmp_laser_pulses=tmp_laser_pulses(tmp_laser_pulses>=trial_start(tr)/samp_freq & tmp_laser_pulses<=(trial_start(tr)/samp_freq+15000))-trial_start(tr)/samp_freq;
        
        
        data_struct(tr).laser_pulses=tmp_laser_pulses;
        data_struct(tr).laser_trial(2)=tmp_laser_pulses(1); %time of first pulse in trial
        
        %calculate frequency of laser pulses
        if length(tmp_laser_pulses)>1
            data_struct(tr).laser_trial(3)=1/(tmp_laser_pulses(2)-tmp_laser_pulses(1))*1000;
        else
            data_struct(tr).laser_trial(3)=0; % 0 Hz, single pulse
        end
    else
        data_struct(tr).laser_pulses=nan;
    end
end

%save several event_indices to the data struct
% and other additional useful variables
data_struct_tmp=data_struct;
data_struct=struct;
data_struct.trials=data_struct_tmp;

% trial type start timestamps
data_struct.laser=laser;
data_struct.laser_start=laser_start;
data_struct.LT_trials=LT_trial;
data_struct.LT_single=LT_single;
data_struct.LT_4hz=LT_4hz;
data_struct.LT_10hz=LT_10hz;
data_struct.LT_40hz=LT_40hz;
data_struct.LT_only=LT_only;

% some useful trial outcome variables
data_struct.outcome=trial_outcome_struct.trial_outcome;
data_struct.lift_flag=trial_outcome_struct.lift_flag;
data_struct.reach_flag=trial_outcome_struct.reach_flag;
data_struct.control_lift_trials=con_lift_trial;
data_struct.control_reach_trials=con_reach_trial;

% all trial outcome variables
data_struct.trial_outcome=trial_outcome_struct;

% trajectory stucture
% X, Y, Z and 3D trajectories, velocity, and acceleration vectors from
%  t=0 to 4000ms, 1 ms bins
% dist_from_con_() = euclidean distance to average success con trial
data_struct.trajectories=trajectory_struct;

clear trajectory_struct trial_outcome_struct

%% Process neural signals
% Read metadata (number of physical channels, sampling rate).
all_fils=dir(neural_dir_ctx);
for i=1:length(all_fils)
    tmp=all_fils(i).name;
    if ~isempty(strfind(tmp,'meta'))
        tmp_fil=strcat(neural_dir_ctx,'\',tmp);
        meta_text = fileread(tmp_fil);
        [foo bar] = regexp(meta_text,'fileTimeSecs=\d');
        fil_len = str2double(meta_text(bar+(0:7)));
        [foo bar] = regexp(meta_text,'niSampRate=\d');
        tmp=splitlines(meta_text(bar+(0:10)));
        samp_freq=str2num(tmp{1})/1000.0;
        
    end
end

if ctx
    % Cortex spikes.
    clear st_ctx
    [spk mua in_ctx ] = get_st_ks2(neural_dir_ctx);
    for j = 1:length(spk)
        st_ctx{j} = double(spk{j})/samp_freq;
    end
    
    % Peri-lift firing rates
    % Get cortical firing rates.
    g_fr = 20; % window size for gaussian smoothing functions
    causal_flag=0; %whether to use a causal kernel or not
    st=st_ctx;
    [rates,t_rates,tstamp_train]=convolve_time_stamps(st,g_fr,fil_len*1000,causal_flag);
    
    % cut up time-series into each trial
    r_lift_ctx=[];
    st_lift_ctx=[];
    trial_start_timestamps=[data_struct.trials(:).trial_start];
    trial_start_timestamps=trial_start_timestamps(2:2:end);
    window_lift = [-5000 10000]; % select data set here
    for t=1:length(trial_start_timestamps)
        wind_min=round(trial_start_timestamps(t)+window_lift(1)+1);
        wind_max=round(trial_start_timestamps(t)+window_lift(2));
        if wind_max<size(rates,2)
            r_lift_ctx(:,t,:)=rates(:,wind_min:wind_max);
            st_lift_ctx(:,t,:)=tstamp_train(:,wind_min:wind_max);
        else
            r_lift_ctx(:,t,:)=[rates(:,wind_min:end),nan(size(rates,1),wind_max-size(rates,2))];
            st_lift_ctx(:,t,:)=[tstamp_train(:,wind_min:end),nan(size(rates,1),wind_max-size(rates,2))];
        end
    end
    clear rates t_rates tstamp_train
    
    % Get z-scores, WITH SOFT NORMALIZATION?
    i4z = 1:400;
    c_softnorm = .5;
    c_soft_z = 5;
    z_norm_lift_ctx=[];
    for i = 1:size(r_lift_ctx,1)
        dd = squeeze(r_lift_ctx(i,:,i4z));
        dd = reshape(dd,1,numel(dd));
        mu_ctx(i) = nanmean(dd);
        std_ctx(i) = nanstd(dd);
        
        %lift aligned
        z_norm_lift_ctx(i,:,:) = (squeeze(r_lift_ctx(i,:,:))-mu_ctx(i))./(c_soft_z + std_ctx(i));
    end
    
    for i=1:length(data_struct.trials)
        data_struct.trials(i).M1_tstamps=squeeze(st_lift_ctx(:,i,:));
        data_struct.trials(i).M1_conv_20ms=squeeze(r_lift_ctx(:,i,:));
        data_struct.trials(i).M1_z_score=squeeze(z_norm_lift_ctx(:,i,:));
        data_struct.trials(i).M1_neural_time=(window_lift(1)+1):window_lift(2);
        data_struct.trials(i).M1_info=in_ctx;
    end
    
    % save data structure to experiment dir
    save(save_file,'data_struct','-v7.3')
    
    % clear variables to restore ram
    clear st_lift_ctx r_lift_ctx z_norm_lift_ctx
end

if thal
    % Thalamus spikes.
    clear st_thal
    [spk mua in_thal] = get_st_ks2(neural_dir_thal);
    for j = 1:length(spk)
        st_thal{j} = double(spk{j})/25;
    end
    
    % thalamus single units
    g_fr = 20; % window size for gaussian smoothing functions
    causal_flag=0;
    st=st_thal;
    [rates,t_rates,tstamp_train]=convolve_time_stamps(st,g_fr,fil_len*1000,causal_flag);
    
    % cut up time-series into each trial
    r_lift_thal=[];
    st_lift_thal=[];
    trial_start_timestamps=[data_struct.trials(:).trial_start];
    trial_start_timestamps=trial_start_timestamps(2:2:end);
    window_lift = [-5000 10000]; % select data set here
    for t=1:length(trial_start_timestamps)
        wind_min=round(trial_start_timestamps(t)+window_lift(1)+1);
        wind_max=round(trial_start_timestamps(t)+window_lift(2));
        r_lift_thal(:,t,:)=rates(:,wind_min:wind_max);
        st_lift_thal(:,t,:)=tstamp_train(:,wind_min:wind_max);
    end
    clear rates t_rates tstamp_train
    
    % Get z-scores, WITH SOFT NORMALIZATION?
    i4z = 1:400;
    c_softnorm = .5;
    c_soft_z = 5;
    z_norm_lift_thal=[];
    for i = 1:size(r_lift_thal,1)
        dd = squeeze(r_lift_thal(i,:,i4z));
        dd = reshape(dd,1,numel(dd));
        mu_thal(i) = nanmean(dd);
        std_thal(i) = nanstd(dd);
        
        %lift aligned
        z_norm_lift_thal(i,:,:) = (squeeze(r_lift_thal(i,:,:))-mu_thal(i))./(c_soft_z + std_thal(i));
    end
    
    for i=1:length(data_struct.trials)
        data_struct.trials(i).Thal_tstamps=squeeze(st_lift_thal(:,i,:));
        data_struct.trials(i).Thal_conv_20ms=squeeze(r_lift_thal(:,i,:));
        data_struct.trials(i).Thal_z_score=squeeze(z_norm_lift_thal(:,i,:));
        data_struct.trials(i).Thal_info=in_thal;
    end
    
    % save data structure to experiment dir
    save(save_file,'data_struct','-v7.3')
    
    clear rates t_rates tstamp_train
    clear st_lift_thal r_lift_thal z_norm_lift_thal
end


end
