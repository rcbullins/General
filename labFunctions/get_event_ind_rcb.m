function get_event_ind_rcb(neural_dir,recsite,channels,expt_type,save_dir)

% recsite=site{i};
% neural_dir=neural_data_dirs{i};
% expt_type=exp_type(i)

% for jcr experiments with from Hantman lab 2017-2021
% whisper system recording trial structure timestamps

% gets event indicies for creating trial structure timestamps

% neural data bin file
% fname = ls([neural_data_dir '\*.nidq.bin']);
fname = ls([neural_dir filesep recsite '\*.nidq.bin*']);
if isempty(fname)
    disp([fname '.bin not in folder'])
else
    
    % neural data meta file
    % mname = ls([neural_data_dir '\*.nidq.meta']);
    mname = ls([neural_dir filesep recsite filesep '*.nidq.meta*']);
    
    % file output, save to...
    % fout = [neural_data_dir(1:(end-length(site))) 'event_ind'];
    % save_dir = [neural_dir filesep 'event_ind'];
    
    % get number of channels, nchan
    metafile = [neural_dir filesep recsite filesep mname(end,:)];
    meta_text = fileread([metafile]);
    [foo bar] = regexp(meta_text,'nSavedChans=\d');
    mystr = meta_text(bar+(0:2));
    if ~isempty(str2num(mystr(end)))
        nchan = str2num(mystr);
    else
        nchan = str2num(mystr(1:(end-1)));
    end
    
    % get sample rate, sr
    [foo bar] = regexp(meta_text,'niSampRate=\d');
    mystr = meta_text(bar+(0:5));
    if ~isempty(str2num(mystr(end)))
        sr = str2num(mystr);
    else
        sr = str2num(mystr(1:(end-1)));
    end
    t_cal = 1e6/sr;
    
    
    trial_start=[];
    sync=[];
    table_start=[];
    laser=[];
    laser_start=[];
    laser_gate=[];
    laser_gate_start=[];
    
    % can get data_tmp for testing
    %     [ind, data_tmp] = get_event_ind(neural_data_dirs,fname,nchan,...
    %                                         [ch1+1 ch2+1 ch3+1 ch4+1 ch5+1],...
    %                                         [0.001 0.001 0.001 0.001 0.001],0);
    
    % NOTE: CHAN # FROM WHISPER (0'd) = CHAN+1 here
    % SESSION TIMESTAMPS
    
    %%%%%%%%%%
    % Get data
    if expt_type==1
        ch1 = channels(1);% 64;% sync pulse timestamps
        ch2 = channels(2);% 66;% camera timestamps
        ch3 = channels(3);% 68;% table_start/cue pulse timestamps
        ch4 = channels(4);% 67;% laser pulse timestamps
        ch5 = channels(5);% 70;% tone pulse timestamps
        ch6 = channels(6);% 69;% masking light pulse timestamp
        
        getSync=0;
        if getSync
            % [ind_neg, data_tmp_neg] = get_event_ind(neural_data_dirs,fname,nchan,[ch6+1],[-0.001],1);
            [ind] = get_event_ind_jcr2([neural_dir filesep recsite filesep],fname,nchan,[ch2+1 ch3+1 ch4+1 ch5+1 ch1+1],[0.001 0.001 0.001 0.001 0.001],0);
            % sync pulse timestamps
            ind_pos = [1 1+find(diff(ind{5})>100e3/t_cal)]; % greater than 200ms
            sync = ind{5}(ind_pos);
        else
            [ind] = get_event_ind_jcr2([neural_dir filesep recsite],fname,nchan,[ch2+1 ch3+1 ch4+1 ch5+1],[0.001 0.001 0.001 0.001],0);
            sync=[];
        end
        
        getLight=0;
        if getLight
            %     [ind_neg, data_tmp_neg] = get_event_ind(neural_data_dirs,fname,nchan,[ch6+1],[-0.001],1);
            [ind_neg] = get_event_ind_jcr2([neural_dir filesep recsite],fname,nchan,[ch6+1],[-0.001],1);
        else
            light=[];
        end
        
        % camera timestamps
        % first timestamp = start of trial_start
        ind_pos = [1 1+find(diff(ind{1})>50e3/t_cal)]; % greater than 50ms us
        trial_start = ind{1}(ind_pos);
        
        % table_start/cue pulse timestamps
        ind_pos = [1 1+find(diff(ind{2})>1e6/t_cal)]; % greater than 1s
        table_start = ind{2}(ind_pos);
        
        % laser pulse timestamps
        ind_pos = [1 1+find(diff(ind{3})>10e3/t_cal)]; % greater than 10ms
        ind_Tpos = [1 1+find(diff(ind{3})>5e6/t_cal)]; % greater than 5s
        laser = ind{3}(ind_pos);
        laser_start = ind{3}(ind_Tpos);
        
        % tone pulse timestamps
        ind_pos = [1 1+find(diff(ind{4})>1e4/t_cal)]; % greater than 1s
        tone = ind{4}(ind_pos);
        
        if getLight % Reagan add if statement here
        % masking led light pulse timestamps
        ind_neg_neg = [1 1+find(diff(ind_neg{1})>1e4/t_cal)]; % greater than 1s
        light = ind_neg{1}(ind_neg_neg);
        end

        disp(['trial_start: ' num2str(length(trial_start))]);
        disp(['table_start: ' num2str(length(table_start))]);
        disp(['laser: ' num2str(length(laser))]);
        disp(['laser_start: ' num2str(length(laser_start))]);
        disp(['tone: ' num2str(length(tone))]);
        if getSync
            disp(['sync: ' num2str(length(sync))]);
        else
            disp('sync: did not get');
        end
        if getLight
            disp(['light: ' num2str(length(light))]);
        else
            disp('light: did not get');
        end
        
    elseif expt_type==2
        ch1 = channels(1);% 64;% sync pulse timestamps
        ch2 = channels(2);% 66;% camera timestamps
        ch3 = channels(3);% 68;% table_start/cue pulse timestamps
        ch4 = channels(4);% 67;% laser pulse timestamps
        ch5 = channels(5);% 70;% laser gate pulse timestamps, for grab detector
        
        getSync=1;
        if getSync
            %     [ind_neg, data_tmp_neg] = get_event_ind(neural_data_dirs,fname,nchan,[ch6+1],[-0.001],1);
            [ind] = get_event_ind_jcr2([neural_dir filesep recsite],fname(end,:),nchan,[ch2+1 ch3+1 ch4+1 ch5+1 ch1+1],[0.001 0.001 0.001 0.001 0.001],0);
            % sync pulse timestamps
            ind_pos = [1 1+find(diff(ind{5})>2e2)]; % greater than 200ms
            sync = ind{5}(ind_pos);
        else
            [ind] = get_event_ind_jcr2([neural_dir filesep recsite],fname(end,:),nchan,[ch2+1 ch3+1 ch4+1 ch5+1],[0.001 0.001 0.001 0.001],0);
            sync=[];
        end
        
        
        % camera timestamps
        % first timestamp = start of trial_start
        ind_pos = [1 1+find(diff(ind{1})>50e3/t_cal)]; % greater than 50ms us
        trial_start = ind{1}(ind_pos);
        
        % table_start/cue pulse timestamps
        ind_pos = [1 1+find(diff(ind{2})>1e6/t_cal)]; % greater than 1s
        table_start = ind{2}(ind_pos);
        
        % laser pulse timestamps
        ind_pos = [1 1+find(diff(ind{3})>10e3/t_cal)]; % greater than 10ms
        ind_Tpos = [1 1+find(diff(ind{3})>5e6/t_cal)]; % greater than 5s
        laser = ind{3}(ind_pos);
        laser_start = ind{3}(ind_Tpos);
        
        % laser gate pulse timestamps
        ind_pos = [1 1+find(diff(ind{4})>10e3/t_cal)]; % greater than 10ms
        ind_Tpos = [1 1+find(diff(ind{4})>5e6/t_cal)]; % greater than 5s
        laser_gate = ind{4}(ind_pos);
        laser_gate_start = ind{4}(ind_Tpos);
        
        disp(['trial_start: ' num2str(length(trial_start))]);
        disp(['table_start: ' num2str(length(table_start))]);
        disp(['laser: ' num2str(length(laser))]);
        disp(['laser_start: ' num2str(length(laser_start))]);
        disp(['laser_gate: ' num2str(length(laser_gate))]);
        disp(['laser_gate_start: ' num2str(length(laser_gate_start))]);
        if getSync
            disp(['sync: ' num2str(length(sync))]);
        else
            disp('sync: did not get');
        end
        
    elseif expt_type==3
        ch1 = channels(1);% 224;% camera timestamps
        ch2 = channels(2);% 226;% table_start/cue pulse timestamps
        ch3 = channels(3);% 230;% laser pulse timestamps
        
        getSync=0; % no sync pulses to get
        
        [ind] = get_event_ind_jcr2([neural_dir filesep recsite],fname(end,:),nchan,[ch1+1 ch2+1 ch3+1],[0.001 0.001 0.001],0);
        sync=[];
        
        %     [ind] = get_event_ind([neural_dir filesep recsite],fname(end,:),nchan,[65],[0.001],0);
        
        % camera timestamps
        % first timestamp = start of trial_start
        ind_pos = [1 1+find(diff(ind{1})>50e3/t_cal)]; % greater than 50ms us
        trial_start = ind{1}(ind_pos);
        
        % table_start/cue pulse timestamps
        ind_pos = [1 1+find(diff(ind{2})>1e6/t_cal)]; % greater than 1s
        table_start = ind{2}(ind_pos);
        
        % laser pulse timestamps
        ind_pos = [1 1+find(diff(ind{3})>10e3/t_cal)]; % greater than 10ms
        ind_Tpos = [1 1+find(diff(ind{3})>5e6/t_cal)]; % greater than 5s
        laser = ind{3}(ind_pos);
        laser_start = ind{3}(ind_Tpos);
        
        % laser gate pulse timestamps
        %     ind_pos = [1 1+find(diff(ind{4})>0.01/t_cal)]; % greater than 10ms
        %     ind_Tpos = [1 1+find(diff(ind{4})>5/t_cal)]; % greater than 5s
        laser_gate = [];%ind{4}(ind_pos);
        laser_gate_start = [];%ind{4}(ind_Tpos);
        %
        disp(['trial_start: ' num2str(length(trial_start))]);
        disp(['table_start: ' num2str(length(table_start))]);
        disp(['laser: ' num2str(length(laser))]);
        disp(['laser_start: ' num2str(length(laser_start))]);
        %     disp(['laser_gate: ' num2str(length(laser_gate))]);
        %     disp(['laser_gate_start: ' num2str(length(laser_gate_start))]);
        if getSync
            disp(['sync: ' num2str(length(sync))]);
        else
            disp('sync: did not get');
        end
        
    end
    
    
    %%%%%%%%%%%
    % Plot data
    plotFig=1;
    if plotFig
        close all;
        figure();
        set(gcf,'Position',[1850 20 1500 1300]);
        
        subplot(6,1,1)
        %     plot(data_tmp(1,:),'Color',[0.7 0.7 0.7])
        hold on
        plot(sync,2e-3*ones(length(sync)),'k.')
        title('sync')
        
        subplot(6,1,2)
        %     plot(data_tmp(2,:),'Color',[0.7 0.7 0.7])
        hold on
        plot(trial_start,2e-3*ones(length(trial_start)),'r*')
        title('camera')
        
        subplot(6,1,3)
        %     plot(data_tmp(3,:),'Color',[0.7 0.7 0.7])
        hold on
        plot(table_start,2e-3*ones(length(table_start)),'r*')
        title('table_start')
        
        subplot(6,1,4)
        %     plot(data_tmp(4,:),'Color',[0.7 0.7 0.7])
        hold on
        plot(tone,2e-3*ones(length(tone)),'r*')
        title('tone')
        
        subplot(6,1,5)
        %     plot(data_tmp(5,:),'Color',[0.7 0.7 0.7])
        plot(laser,2e-3*ones(length(laser)),'c*')
        hold on
        plot(laser_start,2e-3*ones(length(laser_start)),'m*')
        title('laser')
        
        if getLight
            subplot(6,1,6)
            %     plot(data_tmp_neg(1,:),'Color',[0.7 0.7 0.7])
            hold on
            plot(light,-2e-3*ones(length(light)),'b*')
            title('masking light')
        end
        
        zoom on;
        
        %     linkaxesInFigure('x');
    end
end

saveData=1;
if saveData
    saveName = [save_dir filesep 'event_ind'];
    if expt_type==1
        save(saveName,'trial_start','sync','table_start','laser','laser_start','light','tone')
    elseif expt_type==2
        save(saveName,'trial_start','sync','table_start','laser','laser_start','laser_gate','laser_gate_start')
    elseif expt_type==3
        save(saveName,'trial_start','table_start','laser','laser_start')
    end
end
