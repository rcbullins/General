close ALL
addpath('D:\rbullins\Code\decoder');

load('data_for_reagan')
[pred,VAF]=decoder(tmp_kine,tmp_neural);
% tmp_kine =   trial x coordinate x 
% tmp_neural = trial x neuron     x
% plot first 9 trials in the z direction
figure
hold on 
for i=1:9
    subplot(3,3,i)
    hold on 
    %plot the actual data
    tmp=squeeze(tmp_kine(i,1,:));
    plot(tmp,'r')
    
    %plot the predicted/decoded value
    tmp=squeeze(pred(i,1,:));
    plot(tmp,'b')
    legend('Observed','Predicted')
    xlabel('Time (ms)')
    ylabel('Position in Z Direction')
    
end