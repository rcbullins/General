function [traj, conf] = get_traj_3D_RB(data_dir,calib_file,n_frames,i_obj)
% [traj conf] = get_traj_3D(data_dir,calib_file,n_frames)

% Load calibration file.
load(calib_file);

% List files.
fname_side = dir([data_dir '*side*.trk']);
fname_front = dir([data_dir '*front*.trk']);

% Check that the two cameras have the same number of frames.
if length(fname_side) ~= length(fname_front)
    error('Mismatch in the number of videos from the front and side cameras!');
end

% Set up trajectory and confidence arrays.
traj = NaN(length(fname_front),3,n_frames);
conf = NaN(length(fname_front),2,n_frames);

% Loop over trials.
for i = 1:length(fname_side)
    % Load trajectories in camera coordinates.
    tmp_front = load([fname_front(i).folder '/' fname_front(i).name],'-mat');
    tmp_side = load([fname_side(i).folder '/' fname_side(i).name],'-mat');

    if min([size(tmp_side.pTrk,3) size(tmp_front.pTrk,3)]) >= n_frames
        xy_s = squeeze(tmp_side.pTrk(i_obj,:,1:n_frames));
        xy_f = squeeze(tmp_front.pTrk(i_obj,:,1:n_frames));
    else
        xy_s = [squeeze(tmp_side.pTrk(i_obj,:,:)) NaN(2,n_frames-size(tmp_side.pTrk,3))];
        xy_f = [squeeze(tmp_front.pTrk(i_obj,:,:)) NaN(2,n_frames-size(tmp_front.pTrk,3))];
        display(['Missing frames, trial ' num2str(i)])
    end

    % Triangulate to get 3D trajectory.
    traj_tmp = stereo_triangulation(xy_s,xy_f,...
        om,T,fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,...
        kc_right,alpha_c_right);

    % Permute and invert dimensions.
    traj_tmp = traj_tmp([1 3 2],:).*[ones(1,size(traj_tmp,2)); ones(1,size(traj_tmp,2)); -1*ones(1,size(traj_tmp,2))];

    % Store values.
    traj(i,:,:) = traj_tmp;
    conf(i,1,1:min(n_frames,size(tmp_side.pTrk,3))) = tmp_side.pTrkconf(i_obj,1:min(n_frames,size(tmp_side.pTrk,3)));
    conf(i,2,1:min(n_frames,size(tmp_front.pTrk,3))) = tmp_front.pTrkconf(i_obj,1:min(n_frames,size(tmp_front.pTrk,3)));

end

end


