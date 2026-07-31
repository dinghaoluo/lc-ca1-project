x = load(['Z:\Xiaoliang\mice-expdata\ANMZ577\pipeline-data\A577-20210803-06_DataStructure_mazeSection1_TrialType1.mat']);

i = 10;
tr = x.trials{i};

roi1 = 12;
roi2 = 15;

tiledlayout(5, 1)
nexttile

plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.xMM, 'r.')
                    title(['Trial #' num2str(i)]);
                    xlabel('Time (s)')
                    ylabel('Distance (mm)')
                    
     nexttile          
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.speed)     
    xlabel('Time (s)')
    ylabel('Speed (mm/s)')
    
    nexttile          
    hold on;
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.F(:, roi1))
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.F(:, roi2))
hold off;
    xlabel('Time (s)')
    ylabel('F')

     nexttile        
     hold on;
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.Fneu(:, roi1))
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.Fneu(:, roi2))
hold off;
    xlabel('Time (s)')
    ylabel('Fneu')

     nexttile          
     hold on;
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.spks(:, roi1))
plot((tr.lfpIndStart:tr.lfpIndEnd)./100, tr.spks(:, roi2))
hold off;
    xlabel('Time (s)')
    ylabel('Spikes')
