Filenames = [
% ANMP214
               'A214-20221214-02'; 
               'A214-20221215-02';
               'A214-20221215-04';
               'A214-20221216-02';
               'A214-20221216-04';
               'A214-20221217-02';
               'A214-20221218-02';
               'A214-20221218-04';
               'A214-20221219-02';
               'A214-20221220-02';
               'A214-20221220-04';
               'A214-20221221-02';
               'A214-20221221-04';
               'A214-20221222-02';
               'A214-20221222-04';
               'A214-20221222-06';
               'A214-20221223-02';
               'A214-20221223-04';
               'A214-20221224-02';
               'A214-20221224-04';
               'A214-20221225-02';
               'A214-20221225-04';
               'A214-20221226-02';
               'A214-20221226-04';
               'A214-20221226-06'];
       
% Get the length of the string array using length function
nSessionsCount = length(Filenames);

% Extract info from each string
% Merge coresponding running data into images' folder
% 
numFieldArray = zeros(1,nSessionsCount);
numNeuArray = zeros(1,nSessionsCount);
for i = 1:1:nSessionsCount
    % Extract info
    anmlDate = Filenames(i,1:13);
    anmlNo = Filenames(i,2:4);
    session = Filenames(i,15:16);
    
    % Create a folder with images of the session
    targetFolder = ['Z:\Dongyan\mice-expdata\ANMP' anmlNo '\' anmlDate '\' session '\' anmlDate '-' session 'T\'];
    
    load([targetFolder Filenames(i,:) '_DataStructure_mazeSection1_TrialType1_FieldSpCorrAlignedAllTr_Run1_Run1.mat']);
    numFieldArray(i) = length(fieldSpCorrSessNonStimGood.indNeuron);
    
    load([targetFolder Filenames(i,:) '_DataStructure_mazeSection1_TrialType1_FRAlignedRun_msess1_Run1']);
    numNeuArray(i) = length(mFRStructNonStimGood.mFR);
end

figure
subplot(2,1,1)
plot(1:nSessionsCount,numFieldArray,'o');
ylabel('Num. fields')
subplot(2,1,2)
plot(1:nSessionsCount,numNeuArray,'o');
ylabel('Num. neurons')
xlabel('A214 session')