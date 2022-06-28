function update_processed_data_struct(event_ind_dir)

% experiment name for data struct
neural_dir_ctx = [event_ind_dir filesep 'cortex'];
mname = ls([neural_dir_ctx '\*.nidq.meta']);
fname=mname(1:end-10);
save_file = [event_ind_dir filesep 'data_struct_' fname];

% load data structure to experiment dir
load(save_file)

%% load event indices
load([event_ind_dir filesep 'event_ind.mat'])

% save data structure to experiment dir
save(save_file,'data_struct','-v7.3')



end
