% loop over all recordings for control
rec_list = [
    'A029r-20220623-03';
    ...
    'A032r-20220726-02';
    'A032r-20220726-03';
    'A032r-20220727-01';
    'A032r-20220802-01';
    'A032r-20220802-02';
    ...
    'A045r-20221130-02';
    'A045r-20221201-02';
    'A045r-20221201-03';
    'A045r-20221205-03';
    'A045r-20221206-02';
    'A045r-20221207-02';
    'A045r-20221207-04';
    ...
    'A049r-20230103-02';
%     # 'Z:\Dinghao\MiceExp\ANMD049r\A049r-20230103\A049r-20230103-04'; bad clustering for clu 10; the only tagged clu in here
%     # 'Z:\Dinghao\MiceExp\ANMD049r\A049r-20230104\A049r-20230104-02'; similar to last session
%     # 'Z:\Dinghao\MiceExp\ANMD049r\A049r-20230104\A049r-20230104-03'; similar 
    'A049r-20230104-04';
    'A049r-20230120-04';
%     # 'Z:\Dinghao\MiceExp\ANMD049r\A049r-20230121\A049r-20230121-03'; extremely bad behaviour
    ...
%     # 'Z:\Dinghao\MiceExp\ANMD056r\A056r-20230417\A056r-20230417-04'; bad clustering with contamination 
    'A056r-20230418-02';
    'A056r-20230419-02';
    'A056r-20230420-02';
    'A056r-20230420-03';
    'A056r-20230421-02';
    'A056r-20230421-03';
    'A056r-20230422-02'
    ];

tot_rec = length(rec_list);

for i = 1:tot_rec
    filename = rec_list(i,:);
    fullfilename = [filename '_DataStructure_mazeSection1_TrialType1'];
    fullpath = ['Z:\Dinghao\MiceExp\ANMD' filename(2:5) '\' filename(1:14)];
    ProcessingMice_smTrCtrlOnly('./',fullfilename,0,0)
end