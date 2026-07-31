load(['Z:\Xiaoliang\mice-expdata\ANMZ577\pipeline-data\A577-20210803-06_Behav2PDataLFP.mat']);

start = 1;
range = size(Track.xMM, 1)-1;
% roi = 14;
% 
tiledlayout(3,1)
ax1 = nexttile;
plot(start:(start+range), Track.speed_MMsec(start:(start+range)))
ylim([0 inf])
xlim([start start+range])
hold on; 
plot(start:(start+range), Track.xMM(start:(start+range)))
scatter(Laps.startLfpInd,zeros(length(Laps.startLfpInd),1));
hold off; 
% 

ax2 = nexttile;
hold on; 
% plot(start:(start+range), Track.F(start:(start+range), roi))
% plot(start:(start+range), Track.Fneu(start:(start+range), roi))
plot(start:(start+range), Track.spks(start:(start+range), 14))
hold off; 
ylim([0 inf])
xlim([start start+range])

ax3 = nexttile;

plot(start:(start+range), Track.spks(start:(start+range), 1))
ylim([0 inf])
xlim([start start+range])

% image([start start+range],[0 1000],cropped_img);
linkaxes([ax1 ax2 ax3],'x')

% xlim([57620 1943205])
disp('done')