%f_in = 'C:\Users\Hantman Lab\Documents\Jay_Data\V1_recordings\M274_20180910_700um_g0_t0\M274_20180910_700um_g0_t0.nidq.bin';
%f_in = 'C:\Users\Hantman Lab\Documents\Jay_Data\V1_recordings\M261Slc17a7_Gtacr2\M261_20180716_950um\M261_20180716_950um_g0_t0.nidq.bin';
f_in = 'D:\rbullins\Data\M341necab1_Chr2\20210819ephys\M341_20210819_1100_g0\M340_20210807_1000_g0_t0.nidq.bin';
foo = strfind(f_in,'\');
dir_out = f_in(1:(foo(end)-1));
thresh = [.001 .002 .001 .002 .002];

% ( fname,nchan,chan2scan,thresh,invert )
% ind = get_event_ind(f_in,67,65:67,thresh,0);
cd('../probes')
ind = get_event_ind(f_in,5,[1:5],thresh,0);
cd('../')
% i_first_cam=ind{5};

save([dir_out '/ind.mat'],'ind');