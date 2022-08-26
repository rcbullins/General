function [rates,t_rates,train]=convolve_time_stamps(st,window,sess_len,causal_flag)

% convolve timestamps function 'create_processed_dat_struct.m'

    nNeurons=size(st,2);
    bin=1;
    clear rates
    for n=1:nNeurons
        tmp_spks=st{n};
        train(n,:) = histc(tmp_spks,0:bin:sess_len);
        kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
        if causal_flag
            kernel=kernel(1:round(length(kernel)/2));
        end
        rates(n,:)=conv(train(n,:),kernel,'same');
    end
    t_rates=0:1:sess_len; %in ms
    
end