%% create processed data struct

clear ALL;
close ALL;

%% load trk data and calibration file
calib_file='D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\camera_calibration\20191030\Calib_Results_stereo.mat';
trk_files='D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\jcr79\video\20200929\tracked';
load(['D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\jcr79\jcr79_20200929_4600_915um_g0\event_ind.mat']);

spk_dir_ctx='D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\jcr79\jcr79_20200929_4600_915um_g0\cortex';
spk_dir_thal='D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\jcr79\jcr79_20200929_4600_915um_g0\thalamus';


nframes=3000; %number of frames to analyze. Going to depend on how many frames in movie. set at 3000 as an upper bound
samp_freq=25; %should pull from file, but currently hard coded

%names of tracked points
track_names={'digit_2','digit_4','wrist','pellet'};

%data structure to hold all the data
data_struct=struct;

%loop over each tracked point and add to data_struct
for pt=1:4
    [traj conf] = get_traj_3D(trk_files,calib_file,nframes,pt);
     %% assume traj is 3d matrix with trial x coordinate x timepoints
    % break traj down into three components
    x=squeeze(traj(:,1,:));
    y=squeeze(traj(:,2,:));
    z=squeeze(traj(:,3,:));
    conf_side=squeeze(conf(:,1,:));
    conf_front=squeeze(conf(:,2,:));
    nTrials=size(x,1);
    nTPs=size(x,3);

    %interpolate kinematics up to 1kHz
    x_intrp=[];
    y_intrp=[];
    z_intrp=[];
    conf_side_intrp=[];
    conf_front_intrp=[];
    time=1:size(x,2);
    time_intrp=1:0.5:size(x,2)+0.5;
    for tr=1:nTrials
        tmp_x=x(tr,:);
        x_intrp(tr,:)=interp1(time,tmp_x,time_intrp);

        tmp_y=y(tr,:);
        y_intrp(tr,:)=interp1(time,tmp_y,time_intrp);

        tmp_z=z(tr,:);
        z_intrp(tr,:)=interp1(time,tmp_z,time_intrp);
        
        tmp=conf_side(tr,:);
        conf_side_intrp(tr,:)=interp1(time,tmp,time_intrp);
        
        tmp=conf_front(tr,:);
        conf_front_intrp(tr,:)=interp1(time,tmp,time_intrp);
    end

    % pad arrays with Nans to give them the same dimensions as the neural data
    x_intrp=[nan([nTrials,5000]),x_intrp,nan([nTrials,10000-nframes*2])];
    y_intrp=[nan([nTrials,5000]),y_intrp,nan([nTrials,10000-nframes*2])];
    z_intrp=[nan([nTrials,5000]),z_intrp,nan([nTrials,10000-nframes*2])];
    conf_side_intrp=[nan([nTrials,5000]),conf_side_intrp,nan([nTrials,10000-nframes*2])];
    conf_front_intrp=[nan([nTrials,5000]),conf_front_intrp,nan([nTrials,10000-nframes*2])];
    time=-4999:10000;
    for tr=1:nTrials
        data_struct(tr).(strcat(track_names{pt},'_x'))=x_intrp(tr,:);
        data_struct(tr).(strcat(track_names{pt},'_y'))=y_intrp(tr,:);
        data_struct(tr).(strcat(track_names{pt},'_z'))=z_intrp(tr,:);
        data_struct(tr).(strcat(track_names{pt},'_conf_side_tracker'))=conf_side_intrp(tr,:);
        data_struct(tr).(strcat(track_names{pt},'_conf_front_tracker'))=conf_front_intrp(tr,:);
        data_struct(tr).time=time;
    end
end

%loop over data from "event_indices" and assign them to each trial
for tr=1:nTrials
    
    %trial outcome
    data_struct(tr).trial_outcome=trial_outcome(tr);
    
    %trial start
    data_struct(tr).trial_start=[1,trial_start(tr)/samp_freq]; %relative time, global time
    
    % find table turn event
    tmp=find((table_start-trial_start(tr))<30000 & (table_start-trial_start(tr))>0);
    if length(tmp)>0
        data_struct(tr).table=[table_start(tmp)/samp_freq-trial_start(tr)/samp_freq,table_start(tmp)/samp_freq]; %relative time, global time
    else
        data_struct(tr).table=[];
    end
    
    % find tone event
    tmp=find((tone-trial_start(tr))<30000 & (tone-trial_start(tr))>0);
    if length(tmp)>0
        data_struct(tr).tone=[tone(tmp)/samp_freq-trial_start(tr)/samp_freq,tone(tmp)/samp_freq]; %relative time, global time
    else
        data_struct(tr).tone=[];
    end
    
    % find laser trials
    tmp1=find(LT_trial==tr);
    tmp2=find(LT_only==tr);
    if length(tmp1)>0 
        tmp4=find(LT_4hz==tr);
        tmp10=find(LT_10hz==tr);
        tmp40=find(LT_40hz==tr);
        if length(tmp4)>0
            data_struct(tr).laser_trial=[{1},{'4Hz'}];
        elseif length(tmp10)>0
            data_struct(tr).laser_trial=[{1},{'10Hz'}];
        elseif length(tmp40)>0
            data_struct(tr).laser_trial=[{1},{'40Hz'}];
        end
        
        tmp_las=find(laser/samp_freq>trial_start(tr)/samp_freq & laser/samp_freq<(trial_start(tr)/samp_freq+5000));
        tmp_las=laser(tmp_las)/samp_freq;
        tmp_las_rel=tmp_las-trial_start(tr)/samp_freq;
        data_struct(tr).laser_pulses=[tmp_las_rel.',tmp_las.'];

            
    elseif length(tmp2)>0
        data_struct(tr).laser_trial=[{1},{'Laser_Only'}];
        tmp_las=find(laser/samp_freq>trial_start(tr)/samp_freq & laser/samp_freq<(trial_start(tr)/samp_freq+5000));
        tmp_las=laser(tmp_las)/samp_freq;
        tmp_las_rel=tmp_las-trial_start(tr)/samp_freq;
        data_struct(tr).laser_pulses=[tmp_las_rel.',tmp_las.'];
    else
        data_struct(tr).laser_trial=0;
    end
    
    % add reach flag
    data_struct(tr).reach_flag=reach_flag(tr);
    
end

% store these as arrays just because they are useful 
data_struct_tmp=data_struct;
data_struct=struct;
data_struct.trials=data_struct_tmp;
data_struct.trial_outcome=trial_outcome;
data_struct.LT_trials=LT_trial;
data_struct.LT_4hz=LT_4hz;
data_struct.LT_10hz=LT_10hz;
data_struct.LT_40hz=LT_40hz;
data_struct.LT_only=LT_only;
data_struct.reach_flag=reach_flag;
data_struct.control_lift_trials=con_lift_trial;








%% neural data
% Read metadata (number of physical channels, sampling rate).
all_fils=dir(spk_dir_thal);
for i=1:length(all_fils)
    tmp=all_fils(i).name;
    if ~isempty(strfind(tmp,'meta'))
        tmp_fil=strcat(spk_dir_thal,'\',tmp);
        meta_text = fileread(tmp_fil);
        [foo bar] = regexp(meta_text,'fileTimeSecs=\d');
        fil_len = str2double(meta_text(bar+(0:7)));
    end
end


%% Get spikes.
% Cortical spikes.
clear st_ctx st_thal mua_thal
[spk mua in_ctx] = get_st_ks2(spk_dir_ctx);
for j = 1:length(spk)
    st_ctx{j} = double(spk{j})/25; %hard code the sample frequency (need to change)
end
[spk mua in_thal] = get_st_ks2(spk_dir_thal);
for j = 1:length(spk)
    st_thal{j} = double(spk{j})/30;%hard code the sample frequency (need to change)
end



%% Peri-lift firing rates 
% Get cortical firing rates.
g_fr = 20; % window size for gaussian smoothing functions
causal_flag=0;

[rates,t_rates,tstamp_train]=convolve_time_stamps(st_ctx,g_fr,fil_len*1000,causal_flag);

%% cut up time-series into each trial
r_lift_ctx=[];
st_lift_ctx=[];
trial_start_timestamps=[data_struct.trials(:).trial_start];
trial_start_timestamps=trial_start_timestamps(2:2:end);
window_lift = [-5000 10000]; % select data set here
for t=1:length(trial_start_timestamps)
    wind_min=round(trial_start_timestamps(t)+window_lift(1)+1);
    wind_max=round(trial_start_timestamps(t)+window_lift(2));
    r_lift_ctx(:,t,:)=rates(:,wind_min:wind_max);
    st_lift_ctx(:,t,:)=tstamp_train(:,wind_min:wind_max);
end

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




%% Get spikes.
% Cortical spikes.
clear st_ctx st_thal mua_thal
[spk mua in_thal] = get_st_ks2(spk_dir_thal);
for j = 1:length(spk)
    st_thal{j} = double(spk{j})/25;
end

%% thalamus single units
g_fr = 20; % window size for gaussian smoothing functions
causal_flag=0;
[rates,t_rates,tstamp_train]=convolve_time_stamps(st_thal,g_fr,fil_len*1000,causal_flag);

%% cut up time-series into each trial
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






assert(false)
 save('D:\Jay_Data_Kevin_Copy\VL_VA_Stim\thalamus_opto_reach\jcr79\data_struct_jcr79_20200929_4600_915um_g0_neural_kine','data_struct','-v7.3')

assert(false)


function [rates,t_rates,train]=convolve_time_stamps(st_ctx,window,sess_len,causal_flag)

    nNeurons=size(st_ctx,2);
    bin=1;
    clear rates
    for n=1:nNeurons
       tmp_spks=st_ctx{n};
       train(n,:) = histc(tmp_spks,0:bin:sess_len);
       kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
       if causal_flag
%            kernel=kernel(round(length(kernel))/2:end);
           kernel=kernel(1:round(length(kernel)/2));
       end
       rates(n,:)=conv(train(n,:),kernel,'same');
    end
    t_rates=0:1:sess_len; %in ms

end


