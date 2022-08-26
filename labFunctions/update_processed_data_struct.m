function update_processed_data_struct(event_ind_dir,vid_dir)

% experiment name for data struct
neural_dir_ctx = [event_ind_dir filesep 'cortex'];
mname = ls([neural_dir_ctx '\*.nidq.meta']);
fname=mname(1:end-10);
save_file = [event_ind_dir filesep 'data_struct_' fname];

% load data structure to experiment dir
load(save_file)

%% load event indices
load([event_ind_dir filesep 'event_ind.mat'])

%% Include JAABA output
jab_file = dir(fullfile([vid_dir filesep 'jaaba'],'*.mat'));
load([vid_dir filesep 'jaaba' filesep jab_file.name]);

% get JAABA lift, handopen, and grab times from first sequence 'GS00'
% get JAAABA supinate, atmouth, and chew times from last sequence success
% trials only 'GSSS'
for i=1:size(data,1)
    % convert from frame number reference time to 1 ms timestamp
    j_lift(i) = data(i).auto_GS00_Lift_0 * 2;
    j_handopen(i) = data(i).auto_GS00_Handopen_0 * 2;
    j_grasp(i) = data(i).auto_GS00_Grab_0 * 2;
    j_sup(i) = data(i).auto_GSSS_Sup_0 * 2;
    j_atmouth(i) = data(i).auto_GSSS_Atmouth_0 * 2;
    j_chew(i) = data(i).auto_GSSS_Chew_0 * 2;
    if ~isnan(j_chew(i))
        data_struct.jaaba.outcome(i) = 1;
    else
        data_struct.jaaba.outcome(i) = 0;
    end
end

data_struct.jaaba.lift = j_lift;
data_struct.jaaba.handopen = j_handopen;
data_struct.jaaba.grasp = j_grasp;
data_struct.jaaba.atmouth = j_atmouth;
data_struct.jaaba.chew = j_chew;


% save data structure to experiment dir
save(save_file,'data_struct','-v7.3')

disp('update complete')

end
