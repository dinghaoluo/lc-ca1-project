ext = load('Z:\Xiaoliang\mice-expdata\ANMZ577\pipeline-data\A577-20210803-06_DataStructure_mazeSection1_TrialType1_ext.mat');
roi = 5;

tr = 1;

has_quit = false;
while(~has_quit)
    
tiledlayout(3,1)
ax1 = nexttile;
hold on;
plot(ext.trialsExt{tr}.spikesSpeed{1});
yline(20, '--');
hold off;
ax2 = nexttile;
plot(ext.trialsExt{tr}.xMM);

ax3 = nexttile;
plot(ext.trialsExt{tr}.spikes{roi});

answer = questdlg(['Go to next plot? Curr trial: ' string(tr) ' Curr ROI: ' string(roi) ], 'Go to next plot?', ...
	'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        tr = tr + 1;
    case 'No'
        has_quit = true;
end
end