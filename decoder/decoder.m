function [kine_pred_tensor,VAF] = decoder(kine_tensor,neural_tensor)
% decoder used to predict kinematics using neural activity
% inputs are kine_tensor and neural_tensor
% kine_array should be organized as a 3d tensor with trials,cartesian dimensions, time along the 1, 2, and 3, dimensions respectively
% neural_array should be organized the same as kine_array, but along time
% dimension (3) it should contain 14 more time points preceding the first
% kinematic time point.


 [nTrials,nDims,nTPs_kine]=size(kine_tensor);
 [~,~,nTPs_neur]          =size(neural_tensor);
 
 if (nTPs_neur+14)==nTPs_kine
     assert('Neural data does not have enough time points');
 end
 
 %% convert kine tensors to a 2-d matrix
 kine_array=[];
 for i=1:nTrials
    tmp=squeeze(kine_tensor(i,:,:));
    kine_array  =[kine_array tmp];
 end

 % subtract off mean activity
mean_kine_avg=mean(kine_array,2);
kine_array=kine_array-repmat(mean_kine_avg,[1,size(kine_array,2)]);

%% convert neural tensors to 2-d matrix and perform PCA
% get neural activity
neural_array_all=[];
for offset=0:2:14
    n_pc=3;
    % select time window of interest
    neural_array=[];
    
    for i=1:nTrials
        tmp=squeeze(neural_tensor(i,:,offset+1:nTPs_kine+offset));
        neural_array=[neural_array tmp];
    end

    % subtract off mean activity
    mean_z_avg=mean(neural_array,2);
    neural_array=neural_array-repmat(mean_z_avg,[1,size(neural_array,2)]);

    [u,s,v]=svd(neural_array); % all trials and timepoints
    PC=u(:,1:n_pc).';
    neural_proj=PC*neural_array;
    
    neural_array_all=[neural_array_all;neural_proj];
end


%% linear regression
kine_array=kine_array;
b=neural_array_all.'\kine_array.';

kine_array_pred=b.'*neural_array_all;
VAR=norm(kine_array,'fro');
resd=norm(kine_array-kine_array_pred,'fro');
VAF=1-resd.^2/VAR.^2;


%% convert kine_array_pred back to a tensor
kine_array_pred=kine_array_pred+repmat(mean_kine_avg,[1,size(kine_array,2)]);
kine_pred_tensor=[];
 for i=1:nTrials
    tmp=kine_array_pred(:,1+(i-1)*nTPs_kine:nTPs_kine*i);
    kine_pred_tensor(i,:,:)=tmp;
 end
 
 
end

