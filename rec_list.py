# -*- coding: utf-8 -*-
'''
Created on Wed Dec 21 17:53:34 2022

recording lists

@author: Dinghao Luo
'''


#%% paths for SCH pharmacology
pathSCH = [
    'Z:/Dinghao/Behav/DataAnalysis/ANMD006/A006-20210913',  # 006 1st SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD006/A006-20210915',  # 006 2nd SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD005/A005-20210921',  # 005 1st SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD005/A005-20210923',  # 005 2nd SCH session

    'Z:/Dinghao/Behav/DataAnalysis/ANMD009/A009-20211004',  # 009 1st SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD009/A009-20211011',  # 009 2nd SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD009/A009-20211013',  # 009 3rd SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD010/A010-20211014',  # 010 1st SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD011/A011-20211027',  # 011 1st SCH session
    'Z:/Dinghao/Behav/DataAnalysis/ANMD011/A011-20211103'  # 011 2nd SCH session
    ]

sessSCH = [
    [1, 2, 3, 4],  # 1: CONT, 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4],

    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4],  
    [1, 2, 3, 4]
    ]


#%% paths for tagged LC recordings
pathLC = [
    'Z:/Dinghao/MiceExp/ANMD029r/A029r-20220623/A029r-20220623-03',
    
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220726/A032r-20220726-02',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220726/A032r-20220726-03',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220802/A032r-20220802-01',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220802/A032r-20220802-02',
    
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221130/A045r-20221130-02',
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221201/A045r-20221201-02', 
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221201/A045r-20221201-03',
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221205/A045r-20221205-03',
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221206/A045r-20221206-02',
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221207/A045r-20221207-02',
    'Z:/Dinghao/MiceExp/ANMD045r/A045r-20221207/A045r-20221207-04',
    
    'Z:/Dinghao/MiceExp/ANMD049r/A049r-20230103/A049r-20230103-02',
    'Z:/Dinghao/MiceExp/ANMD049r/A049r-20230104/A049r-20230104-04',
    'Z:/Dinghao/MiceExp/ANMD049r/A049r-20230120/A049r-20230120-04',
    
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230418/A056r-20230418-02',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230419/A056r-20230419-02',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230420/A056r-20230420-02',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230420/A056r-20230420-03',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230421/A056r-20230421-02',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230421/A056r-20230421-03',
    'Z:/Dinghao/MiceExp/ANMD056r/A056r-20230422/A056r-20230422-02',
    
    'Z:/Dinghao/MiceExp/ANMD060r/A060r-20230530/A060r-20230530-02',
    'Z:/Dinghao/MiceExp/ANMD060r/A060r-20230602/A060r-20230602-01',
    'Z:/Dinghao/MiceExp/ANMD060r/A060r-20230605/A060r-20230605-01',
    
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230626/A062r-20230626-01',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230626/A062r-20230626-02',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230627/A062r-20230627-02',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230629/A062r-20230629-01',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230629/A062r-20230629-02',
    
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230726/A065r-20230726-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230727/A065r-20230727-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230728/A065r-20230728-02',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230729/A065r-20230729-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230801/A065r-20230801-01',
    
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230821/A067r-20230821-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230821/A067r-20230821-02', 
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230822/A067r-20230822-01', 
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230823/A067r-20230823-01', 
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230823/A067r-20230823-02',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230824/A067r-20230824-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230824/A067r-20230824-02',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230825/A067r-20230825-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230825/A067r-20230825-02'
    ]


#%% paths for optogenetic LC recordings
pathLCopt = [
    # old protocol: high frequency stim
    # 'Z:\Dinghao\MiceExp\ANMD056r\A056r-20230420\A056r-20230420-02',
    # 'Z:\Dinghao\MiceExp\ANMD056r\A056r-20230420\A056r-20230420-03',

    # new protocols from here onwards--.5 s, 12 Hz
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230626/A062r-20230626-01',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230626/A062r-20230626-02',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230629/A062r-20230629-01',
    'Z:/Dinghao/MiceExp/ANMD062r/A062r-20230629/A062r-20230629-02',
    
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230726/A065r-20230726-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230727/A065r-20230727-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230729/A065r-20230729-01',
    'Z:/Dinghao/MiceExp/ANMD065r/A065r-20230801/A065r-20230801-01',
    
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230821/A067r-20230821-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230821/A067r-20230821-02',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230823/A067r-20230823-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230823/A067r-20230823-02',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230824/A067r-20230824-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230824/A067r-20230824-02',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230825/A067r-20230825-01',
    'Z:/Dinghao/MiceExp/ANMD067r/A067r-20230825/A067r-20230825-02'
    ]


#%% paths for optogenetic LC behaviour
pathLCBehopt = [
    'Z:/Dinghao/Behav/DataAnalysis/ANMD014/A014-20211201',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD015/A015-20211208',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD014/A014-20211213',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD015/A015-20220118',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD015/A015-20220120',
    
    'Z:/Dinghao/Behav/DataAnalysis/ANMD052/A052-20230327',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD052/A052-20230328',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD052/A052-20230329',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD052/A052-20230330',

    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230503',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230504', 
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230505',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230510',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230512',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230518', 
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230519', 
    'Z:/Dinghao/Behav/DataAnalysis/ANMD059/A059-20230523',

    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230619',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230620',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230621',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230624',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230628', 
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230630',
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230703', 
    'Z:/Dinghao/Behav/DataAnalysis/ANMD061/A061-20230705'
    ]

sessLCBehopt = [
    [2, 2, 2],
    [2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 2, 3, 1, 2, 3],
    [1, 2, 3, 1, 2, 3],
    
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1],
    [1, 1, 1, 2, 2, 2, 3, 3, 3],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1], 

    [1, 1, 1],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1, 2, 2, 2],
    [1, 1, 1],
    [1, 1, 1],
    [1, 1, 1],
    [1, 1, 1, 2, 2, 2]
    ]

condLCBehopt = [
    [0, 3, 0],
    [0, 3, 0],
    [0, 4, 10, 0, 3, 0],
    [0, 3, 0, 0, 4, 10],
    [0, 3, 0, 0, 4, 10],
    
    [0, 3, 0],
    [0, 3, 0, 0, 4, 10],
    [0, 3, 0, 0, 4, 10],
    [0, 3, 0, 0, 4, 10],

    [0, 4, 0, 0, 4, 0], 
    [0, 2, 0],        
    [0, 4, 0, 0, 3, 0, 0, 3, 0],  
    [0, 4, 0, 0, 2, 0], 
    [0, 2, 0],         
    [0, 4, 0, 0, 2, 0], 
    [0, 11, 0, 0, 11, 0],  
    [0, 3, 0],      

    [0, 2, 0],       
    [0, 2, 0, 0, 2, 0], 
    [0, 2, 0, 0, 2, 0], 
    [0, 2, 0, 0, 1, 0], 
    [0, 2, 0],      
    [0, 2, 0],         
    [0, 2, 0],        
    [0, 2, 0, 0, 1, 0]  
    ]


#%% paths for hippocampus recording with LC stim 
pathHPCLCopt = [
    'Z:/Dinghao/MiceExp/ANMD063r/A063r-20230706/A063r-20230706-01',
    'Z:/Dinghao/MiceExp/ANMD063r/A063r-20230708/A063r-20230708-01',
    'Z:/Dinghao/MiceExp/ANMD063r/A063r-20230708/A063r-20230708-02',
    'Z:/Dinghao/MiceExp/ANMD063r/A063r-20230713/A063r-20230713-01',
    
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230905/A069r-20230905-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230905/A069r-20230905-02',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230908/A069r-20230908-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230908/A069r-20230908-02',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230909/A069r-20230909-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230909/A069r-20230909-02',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230913/A069r-20230913-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230914/A069r-20230914-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230915/A069r-20230915-01',
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230915/A069r-20230915-02',
    
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230920/A071r-20230920-02',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230921/A071r-20230921-01', 
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230921/A071r-20230921-02',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230922/A071r-20230922-01',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230922/A071r-20230922-02',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230923/A071r-20230923-01',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230923/A071r-20230923-02',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230926/A071r-20230926-01',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230928/A071r-20230928-01',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230928/A071r-20230928-02',
    
    'Z:/Dinghao/MiceExp/ANMD068r/A068r-20231024/A068r-20231024-01', 
    'Z:/Dinghao/MiceExp/ANMD068r/A068r-20231025/A068r-20231025-01', 
    'Z:/Dinghao/MiceExp/ANMD068r/A068r-20231026/A068r-20231026-01'
    ]


#%% paths for hippocampus recording with LC-HPC fibre stim 
pathHPCLCtermopt = [
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231108/A070r-20231108-01',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231109/A070r-20231109-01', 
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231110/A070r-20231110-01',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231111/A070r-20231111-01', 
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231115/A070r-20231115-01', 
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231116/A070r-20231116-01',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231117/A070r-20231117-01', 
    
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231212/A076r-20231212-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231213/A076r-20231213-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231214/A076r-20231214-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231215/A076r-20231215-01', 
    
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240124/A078r-20240124-01', 
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240125/A078r-20240125-01',
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240129/A078r-20240129-01',
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240130/A078r-20240130-01',
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240131/A078r-20240131-01', 
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240201/A078r-20240201-01', 
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240202/A078r-20240202-01', 
    
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240207/A077r-20240207-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240209/A077r-20240209-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240215/A077r-20240215-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240216/A077r-20240216-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240217/A077r-20240217-01',
    
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240307/A083r-20240307-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240308/A083r-20240308-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240309/A083r-20240309-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240311/A083r-20240311-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240312/A083r-20240312-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240314/A083r-20240314-01'
    ]


#%% HPC recordings with bad behaviour (exclude from behavioural analyses)
pathHPCbadbeh = [
    'Z:/Dinghao/MiceExp/ANMD069r/A069r-20230914/A069r-20230914-01',
    'Z:/Dinghao/MiceExp/ANMD071r/A071r-20230921/A071r-20230921-01',
    'Z:/Dinghao/MiceExp/ANMD068r/A068r-20231024/A068r-20231024-01',
    'Z:/Dinghao/MiceExp/ANMD068r/A068r-20231024/A068r-20231024-02',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231109/A070r-20231109-01',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231111/A070r-20231111-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231215/A076r-20231215-01',
    'Z:/Dinghao/MiceExp/ANMD078r/A078r-20240202/A078r-20240202-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240314/A083r-20240314-01'
    ]


#%% behaviour only
pathHPCLCtermopt_beh = [
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231108/A070r-20231108-01',
    'Z:/Dinghao/MiceExp/ANMD070r/A070r-20231111/A070r-20231111-01',

    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231213/A076r-20231213-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231214/A076r-20231214-01',
    'Z:/Dinghao/MiceExp/ANMD076r/A076r-20231215/A076r-20231215-01',
    
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240207/A077r-20240207-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240209/A077r-20240209-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240214/A077r-20240214-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240215/A077r-20240215-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240216/A077r-20240216-01',
    'Z:/Dinghao/MiceExp/ANMD077r/A077r-20240217/A077r-20240217-01',
    
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240308/A083r-20240308-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240309/A083r-20240309-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240311/A083r-20240311-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240312/A083r-20240312-01',
    'Z:/Dinghao/MiceExp/ANMD083r/A083r-20240314/A083r-20240314-01'
    ]


#%% grabne
pathHPCGRABNE = [
    
    'Z:/Dinghao/2p_recording/A093i/A093i-20240620/A093i-20240620-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240621/A093i-20240621-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240624/A093i-20240624-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240625/A093i-20240625-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240625/A093i-20240625-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240626/A093i-20240626-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240626/A093i-20240626-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240627/A093i-20240627-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240628/A093i-20240628-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240628/A093i-20240628-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240629/A093i-20240629-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240629/A093i-20240629-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240630/A093i-20240630-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240630/A093i-20240630-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240701/A093i-20240701-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240701/A093i-20240701-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240702/A093i-20240702-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240702/A093i-20240702-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240703/A093i-20240703-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240703/A093i-20240703-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240704/A093i-20240704-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240704/A093i-20240704-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240705/A093i-20240705-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240705/A093i-20240705-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240708/A093i-20240708-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240708/A093i-20240708-02',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240712/A093i-20240712-01',
    'Z:/Dinghao/2p_recording/A093i/A093i-20240712/A093i-20240712-02',
    
    'Z:/Dinghao/2p_recording/A094i/A094i-20240701/A094i-20240701-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240702/A094i-20240702-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240703/A094i-20240703-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240703/A094i-20240703-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240704/A094i-20240704-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240704/A094i-20240704-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240705/A094i-20240705-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240705/A094i-20240705-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240709/A094i-20240709-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240710/A094i-20240710-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240711/A094i-20240711-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240712/A094i-20240712-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240715/A094i-20240715-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240716/A094i-20240716-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240716/A094i-20240716-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240717/A094i-20240717-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240718/A094i-20240718-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240718/A094i-20240718-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240719/A094i-20240719-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240719/A094i-20240719-02',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240807/A094i-20240807-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240808/A094i-20240808-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240809/A094i-20240809-01',
    'Z:/Dinghao/2p_recording/A094i/A094i-20240809/A094i-20240809-02',
    
    'Z:/Dinghao/2p_recording/A097i/A097i-20240814/A097i-20240814-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240815/A097i-20240815-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240815/A097i-20240815-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240816/A097i-20240816-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240816/A097i-20240816-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240819/A097i-20240819-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240819/A097i-20240819-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240820/A097i-20240820-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240821/A097i-20240821-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240821/A097i-20240821-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240822/A097i-20240822-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240822/A097i-20240822-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240823/A097i-20240823-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240823/A097i-20240823-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240826/A097i-20240826-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240826/A097i-20240826-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240827/A097i-20240827-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240827/A097i-20240827-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240828/A097i-20240828-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240829/A097i-20240829-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240829/A097i-20240829-02',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240830/A097i-20240830-01',
    'Z:/Dinghao/2p_recording/A097i/A097i-20240830/A097i-20240830-02',
    
    'Z:/Dinghao/2p_recording/A098i/A098i-20240923/A098i-20240923-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240924/A098i-20240924-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240925/A098i-20240925-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240926/A098i-20240926-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240927/A098i-20240927-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240927/A098i-20240927-02',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240930/A098i-20240930-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20240930/A098i-20240930-02',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241002/A098i-20241002-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241003/A098i-20241003-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241004/A098i-20241004-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241007/A098i-20241007-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241008/A098i-20241008-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241011/A098i-20241011-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241014/A098i-20241014-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241015/A098i-20241015-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241016/A098i-20241016-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241017/A098i-20241017-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241018/A098i-20241018-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241021/A098i-20241021-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241022/A098i-20241022-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241023/A098i-20241023-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241024/A098i-20241024-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241024/A098i-20241024-02',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241025/A098i-20241025-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241025/A098i-20241025-02',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241028/A098i-20241028-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241028/A098i-20241028-02',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241030/A098i-20241030-01',
    'Z:/Dinghao/2p_recording/A098i/A098i-20241030/A098i-20241030-02'
    ]


#%% paths for GCaMP-axon recordings
pathLCHPCGCaMP = [
    'Z:/Dinghao/2p_recording/A101i/A101i-20241029/A101i-20241029-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241031/A101i-20241031-01',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241031/A101i-20241031-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241101/A101i-20241101-01',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241101/A101i-20241101-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241105/A101i-20241105-01',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241105/A101i-20241105-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241106/A101i-20241106-01',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241106/A101i-20241106-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241106/A101i-20241106-03',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241107/A101i-20241107-01',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241107/A101i-20241107-02',
    'Z:/Dinghao/2p_recording/A101i/A101i-20241107/A101i-20241107-03',
    
    'Z:/Dinghao/2p_recording/A105i/A105i-20241122/A105i-20241122-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241125/A105i-20241125-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241125/A105i-20241125-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241126/A105i-20241126-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241126/A105i-20241126-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241127/A105i-20241127-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241127/A105i-20241127-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241128/A105i-20241128-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241202/A105i-20241202-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241202/A105i-20241202-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241203/A105i-20241203-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241204/A105i-20241204-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241204/A105i-20241204-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241205/A105i-20241205-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241205/A105i-20241205-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241206/A105i-20241206-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241206/A105i-20241206-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241206/A105i-20241206-03',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241209/A105i-20241209-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241209/A105i-20241209-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241209/A105i-20241209-03',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241210/A105i-20241210-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241211/A105i-20241211-01',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241211/A105i-20241211-02',
    'Z:/Dinghao/2p_recording/A105i/A105i-20241211/A105i-20241211-03',
    
    'Z:/Dinghao/2p_recording/A106i/A106i-20250122/A106i-20250122-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250122/A106i-20250122-02',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250123/A106i-20250123-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250123/A106i-20250123-02',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250127/A106i-20250127-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250127/A106i-20250127-02',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250128/A106i-20250128-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250128/A106i-20250128-02',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250129/A106i-20250129-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250129/A106i-20250129-02',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250130/A106i-20250130-01',
    'Z:/Dinghao/2p_recording/A106i/A106i-20250130/A106i-20250130-02'
    ]

pathLCHPCGCaMPImmobile = [
    'Z:/Dinghao/2p_recording/A132i/A132i-20250711/A132i-20250711-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250711/A132i-20250711-02',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250714/A132i-20250714-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250714/A132i-20250714-02',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250714/A132i-20250714-03',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250715/A132i-20250715-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250715/A132i-20250715-02',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250715/A132i-20250715-03',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250716/A132i-20250716-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250717/A132i-20250717-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250717/A132i-20250717-02',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250717/A132i-20250717-03',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250718/A132i-20250718-01',
    'Z:/Dinghao/2p_recording/A132i/A132i-20250718/A132i-20250718-02',
    
    'Z:/Dinghao/2p_recording/A134i/A134i-20250714/A134i-20250714-01',
    'Z:/Dinghao/2p_recording/A134i/A134i-20250715/A134i-20250715-01'
    ]


#%% paths for VTA-CA1 GCaMP-axon
pathVTAHPCGCaMP = [
    'Z:/Dinghao/2p_recording/A102i/A102i-20250307/A102i-20250307-01',
    'Z:/Dinghao/2p_recording/A102i/A102i-20250307/A102i-20250307-02',
    ]


#%% paths for dLight + opto-LC stim.
pathdLightLCOpto = [
    'Z:/Dinghao/2p_recording/A114i/A114i-20250327/A114i-20250327-02',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250331/A114i-20250331-01',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-04',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-05',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-06',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250402/A114i-20250402-01',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250402/A114i-20250402-02',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250402/A114i-20250402-03',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250403/A114i-20250403-01',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250403/A114i-20250403-02',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250403/A114i-20250403-03',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250404/A114i-20250404-01',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250404/A114i-20250404-02',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250404/A114i-20250404-03',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250423/A114i-20250423-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A114i/A114i-20250423/A114i-20250423-02',  # with 1100-nm fibre ref.
    
    'Z:/Dinghao/2p_recording/A116i/A116i-20250409/A116i-20250409-01',
    'Z:/Dinghao/2p_recording/A116i/A116i-20250409/A116i-20250409-02',
    'Z:/Dinghao/2p_recording/A116i/A116i-20250414/A116i-20250414-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250414/A116i-20250414-02',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250414/A116i-20250414-03',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250415/A116i-20250415-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250415/A116i-20250415-02',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250415/A116i-20250415-03',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250416/A116i-20250416-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250416/A116i-20250416-02',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250417/A116i-20250417-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250418/A116i-20250418-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250418/A116i-20250418-02',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250418/A116i-20250418-03',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250421/A116i-20250421-01',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250421/A116i-20250421-02',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250421/A116i-20250421-03',  # with 1100-nm fibre ref.
    'Z:/Dinghao/2p_recording/A116i/A116i-20250422/A116i-20250422-01',  # with 1100-nm fibre ref.
    
    'Z:/Dinghao/2p_recording/A126i/A126i-20250605/A126i-20250605-02',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250605/A126i-20250605-03',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250606/A126i-20250606-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250606/A126i-20250606-02',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250609/A126i-20250609-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250609/A126i-20250609-02',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250610/A126i-20250610-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250610/A126i-20250610-02',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250611/A126i-20250611-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250612/A126i-20250612-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250616/A126i-20250616-01',
    'Z:/Dinghao/2p_recording/A126i/A126i-20250617/A126i-20250617-01',
    
    'Z:/Dinghao/2p_recording/A136i/A136i-20250724/A136i-20250724-01',
    'Z:/Dinghao/2p_recording/A136i/A136i-20250724/A136i-20250724-02', 
    'Z:/Dinghao/2p_recording/A136i/A136i-20250725/A136i-20250725-01',
    'Z:/Dinghao/2p_recording/A136i/A136i-20250728/A136i-20250728-01',
    'Z:/Dinghao/2p_recording/A136i/A136i-20250728/A136i-20250728-02', 
    'Z:/Dinghao/2p_recording/A136i/A136i-20250729/A136i-20250729-01',
    'Z:/Dinghao/2p_recording/A136i/A136i-20250729/A136i-20250729-02', 
    'Z:/Dinghao/2p_recording/A136i/A136i-20250730/A136i-20250730-02', 
    'Z:/Dinghao/2p_recording/A136i/A136i-20250731/A136i-20250731-01', 
    'Z:/Dinghao/2p_recording/A136i/A136i-20250801/A136i-20250801-01', 
    
    'Z:/Dinghao/2p_recording/A146i/A146i-20250908/A146i-20250908-01',
    'Z:/Dinghao/2p_recording/A146i/A146i-20250909/A146i-20250909-01',
    'Z:/Dinghao/2p_recording/A146i/A146i-20250910/A146i-20250910-01',
    'Z:/Dinghao/2p_recording/A146i/A146i-20250911/A146i-20250911-01',
    'Z:/Dinghao/2p_recording/A146i/A146i-20250911/A146i-20250911-02',
    'Z:/Dinghao/2p_recording/A146i/A146i-20250912/A146i-20250912-01',
        
    'Z:/Dinghao/2p_recording/A152i/A152i-20251020/A152i-20251020-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251020/A152i-20251020-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251021/A152i-20251021-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251021/A152i-20251021-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251022/A152i-20251022-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251022/A152i-20251022-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251023/A152i-20251023-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251023/A152i-20251023-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251027/A152i-20251027-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251027/A152i-20251027-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251028/A152i-20251028-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251113/A152i-20251113-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251113/A152i-20251113-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251113/A152i-20251113-03',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251114/A152i-20251114-01'
    ]

pathdLightLCOptoDbhBlock = [
    'Z:/Dinghao/2p_recording/A152i/A152i-20251029/A152i-20251029-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251029/A152i-20251029-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251030/A152i-20251030-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251030/A152i-20251030-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251030/A152i-20251030-03',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251031/A152i-20251031-01',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251031/A152i-20251031-02',
    'Z:/Dinghao/2p_recording/A152i/A152i-20251031/A152i-20251031-03'
    ]

pathdLightLCOptoInh = [
    'Z:/Dinghao/2p_recording/A127i/A127i-20250616/A127i-20250616-01',
    'Z:/Dinghao/2p_recording/A127i/A127i-20250701/A127i-20250701-01',
    'Z:/Dinghao/2p_recording/A127i/A127i-20250704/A127i-20250704-01',
    
    'Z:/Dinghao/2p_recording/A129i/A129i-20250616/A129i-20250616-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250624/A129i-20250624-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250630/A129i-20250630-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250708/A129i-20250708-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250715/A129i-20250715-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250716/A129i-20250716-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250718/A129i-20250718-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250721/A129i-20250721-01',
    'Z:/Dinghao/2p_recording/A129i/A129i-20250722/A129i-20250722-01',
    
    'Z:/Dinghao/2p_recording/A133i/A133i-20250625/A133i-20250625-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250625/A133i-20250625-02',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250701/A133i-20250701-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250702/A133i-20250702-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250703/A133i-20250703-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250703/A133i-20250703-02',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250704/A133i-20250704-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250704/A133i-20250704-02',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250707/A133i-20250707-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250707/A133i-20250707-02',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250707/A133i-20250707-03',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250708/A133i-20250708-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250708/A133i-20250708-02',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250710/A133i-20250710-01',
    'Z:/Dinghao/2p_recording/A133i/A133i-20250710/A133i-20250710-02',
    
    'Z:/Dinghao/2p_recording/A135i/A135i-20250723/A135i-20250723-01',
    'Z:/Dinghao/2p_recording/A135i/A135i-20250724/A135i-20250724-01',
    'Z:/Dinghao/2p_recording/A135i/A135i-20250725/A135i-20250725-01',
    
    'Z:/Dinghao/2p_recording/A138i/A138i-20250729/A138i-20250729-01',
    
    'Z:/Dinghao/2p_recording/A140i/A140i-20250725/A140i-20250725-03',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250728/A140i-20250728-01',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250728/A140i-20250728-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250729/A140i-20250729-01',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250729/A140i-20250729-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250730/A140i-20250730-01', 
    'Z:/Dinghao/2p_recording/A140i/A140i-20250730/A140i-20250730-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250731/A140i-20250731-01', 
    'Z:/Dinghao/2p_recording/A140i/A140i-20250731/A140i-20250731-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250731/A140i-20250731-03',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250801/A140i-20250801-01', 
    'Z:/Dinghao/2p_recording/A140i/A140i-20250801/A140i-20250801-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250801/A140i-20250801-03',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250804/A140i-20250804-01', 
    'Z:/Dinghao/2p_recording/A140i/A140i-20250804/A140i-20250804-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250804/A140i-20250804-03',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250805/A140i-20250805-01', 
    'Z:/Dinghao/2p_recording/A140i/A140i-20250805/A140i-20250805-02',
    'Z:/Dinghao/2p_recording/A140i/A140i-20250805/A140i-20250805-03',
    
    'Z:/Dinghao/2p_recording/A150i/A150i-20251119/A150i-20251119-01',
    'Z:/Dinghao/2p_recording/A150i/A150i-20251124/A150i-20251124-01',
    'Z:/Dinghao/2p_recording/A150i/A150i-20251124/A150i-20251124-02'
    ]

pathdLightLCOptoCtrl = [
    # control--sleeve blocked 
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-01',
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-02', 
    'Z:/Dinghao/2p_recording/A114i/A114i-20250401/A114i-20250401-03',
        
    'Z:/Dinghao/2p_recording/A126i/A126i-20250627/A126i-20250627-01',
    
    'Z:/Dinghao/2p_recording/A130i/A130i-20250627/A130i-20250627-01',
    'Z:/Dinghao/2p_recording/A130i/A130i-20250627/A130i-20250627-02', 
    'Z:/Dinghao/2p_recording/A130i/A130i-20250627/A130i-20250627-03',
    'Z:/Dinghao/2p_recording/A130i/A130i-20250627/A130i-20250627-04' 
    ]


#%% paths for GRABNE + opto-LC stim.
pathGRABNELCOpto = [
    'Z:/Dinghao/2p_recording/A147i/A147i-20251029/A147i-20251029-01',
    'Z:/Dinghao/2p_recording/A147i/A147i-20251029/A147i-20251029-02',
    'Z:/Dinghao/2p_recording/A147i/A147i-20251029/A147i-20251029-03',
    'Z:/Dinghao/2p_recording/A147i/A147i-20251030/A147i-20251030-01',
    'Z:/Dinghao/2p_recording/A147i/A147i-20251030/A147i-20251030-02',
    'Z:/Dinghao/2p_recording/A147i/A147i-20251030/A147i-20251030-03',
    
    'Z:/Dinghao/2p_recording/A148i/A148i-20251029/A148i-20251029-01',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251029/A148i-20251029-02',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251030/A148i-20251030-01',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251030/A148i-20251030-02',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251117/A148i-20251117-02',
    
    'Z:/Dinghao/2p_recording/A162i/A162i-20251216/A162i-20251216-01',
    'Z:/Dinghao/2p_recording/A162i/A162i-20251216/A162i-20251216-02',
    'Z:/Dinghao/2p_recording/A162i/A162i-20251218/A162i-20251218-01',
    'Z:/Dinghao/2p_recording/A162i/A162i-20251218/A162i-20251218-02',
    'Z:/Dinghao/2p_recording/A162i/A162i-20251218/A162i-20251218-03',
    'Z:/Dinghao/2p_recording/A162i/A162i-20251218/A162i-20251218-04',
    
    'Z:/Dinghao/2p_recording/A170i/A170i-20260116/A170i-20260116-01',
    'Z:/Dinghao/2p_recording/A170i/A170i-20260116/A170i-20260116-02'
    ]

pathGRABNELCOptoDbhBlock = [
    'Z:/Dinghao/2p_recording/A148i/A148i-20251117/A148i-20251117-01'
    ]


#%% paths for GRABNE + tone activation
pathGRABNETone = [
    'Z:/Dinghao/2p_recording/A148i/A148i-20251107/A148i-20251107-01'
    ]

pathGRABNEToneDbhBlock = [
    'Z:/Dinghao/2p_recording/A148i/A148i-20251110/A148i-20251110-01',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251110/A148i-20251110-02',
    'Z:/Dinghao/2p_recording/A148i/A148i-20251110/A148i-20251110-03'
    ]


#%% paths for nLight + activation LC
pathnLightLCOpto = [
    'Z:/Dinghao/2p_recording/A171i/A171i-20260116/A171i-20260116-01',
    
    'Z:/Dinghao/2p_recording/A171i/A171i-20260220/A171i-20260220-01',
    'Z:/Dinghao/2p_recording/A171i/A171i-20260220/A171i-20260220-02',
    
    'Z:/Dinghao/2p_recording/A171i/A171i-20260227/A171i-20260227-01',
    'Z:/Dinghao/2p_recording/A171i/A171i-20260227/A171i-20260227-02',

    'Z:/Dinghao/2p_recording/A191i/A191i-20260423/A191i-20260423-01',
    'Z:/Dinghao/2p_recording/A191i/A191i-20260428/A191i-20260428-01',
    'Z:/Dinghao/2p_recording/A191i/A191i-20260430/A191i-20260430-01',

    'Z:/Dinghao/2p_recording/A192i/A192i-20260508/A192i-20260508-01',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260508/A192i-20260508-02',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260508/A192i-20260508-03',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260511/A192i-20260511-01',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260511/A192i-20260511-02',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260513/A192i-20260513-01',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260514/A192i-20260514-01',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260514/A192i-20260514-02',
    'Z:/Dinghao/2p_recording/A192i/A192i-20260514/A192i-20260514-03',

    'Z:/Dinghao/2p_recording/A193i/A193i-20260429/A193i-20260429-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260429/A193i-20260429-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260430/A193i-20260430-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260430/A193i-20260430-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260506/A193i-20260506-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260506/A193i-20260506-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260506/A193i-20260506-03',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260508/A193i-20260508-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260508/A193i-20260508-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260511/A193i-20260511-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260511/A193i-20260511-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260511/A193i-20260511-03',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260512/A193i-20260512-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260512/A193i-20260512-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260512/A193i-20260512-03',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260513/A193i-20260513-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260513/A193i-20260513-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260514/A193i-20260514-01',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260514/A193i-20260514-02',
    'Z:/Dinghao/2p_recording/A193i/A193i-20260514/A193i-20260514-03'
    ]


#%% paths for tagged LC recordings
pathVTA = [
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220726/A032r-20220726-02',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220726/A032r-20220726-03',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220727/A032r-20220727-01',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220802/A032r-20220802-01',
    'Z:/Dinghao/MiceExp/ANMD032r/A032r-20220802/A032r-20220802-02',
    
    'Z:/Dinghao/MiceExp/ANMD035r/A035r-20220825/A035r-20220825-01',
    'Z:/Dinghao/MiceExp/ANMD035r/A035r-20220825/A035r-20220825-02',
    'Z:/Dinghao/MiceExp/ANMD035r/A035r-20220826/A035r-20220826-01',
    'Z:/Dinghao/MiceExp/ANMD035r/A035r-20220831/A035r-20220831-02'
    ]


#%% immobile sessions 
pathIm = [
    'Z:/Dinghao/MiceExp/ANMD107/A107-20250121-01T.txt',
    'Z:/Dinghao/MiceExp/ANMD108/A108-20250121-01T.txt']


#%% NE blocker sessions 
pathAlphaBlocker = [
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250610',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250616',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250619',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250626',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250702',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250703'
    ]

sessAlphaBlocker = [
    [1, 2],  # 1: CONT., 2: alpha 0.5 hr, 3: alpha 1 hr
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3]
    ]

pathBetaBlocker = [
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250624',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250627', 
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250630',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250701',
    'Z:/Dinghao/MiceExp/ANMD122/A122-20250708'
    ]

sessBetaBlocker = [
    [1, 2, 3],  # 1: CONT., 2: beta 0.5 hr, 3: beta 1 hr
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3]
    ]

pathXiaoliangBetaBlocker = [
    'Z:/Yingxue/DataAnalysisXiaoliang/ANMZ508/A508-20190923',
    'Z:/Yingxue/DataAnalysisXiaoliang/ANMZ510/A510-20190924'
    ]

sessXiaoliangBetaBlocker = [
    [1,3,4,5],
    [2,4,5]
    ]


#%% other people's rec lists 
pathHPC_Raphi = [
    'Z:/Raphael_tests/mice_expdata/ANM011/A011-20190218/A011-20190218-01',
    'Z:/Raphael_tests/mice_expdata/ANM011/A011-20190219/A011-20190219-01',
    'Z:/Raphael_tests/mice_expdata/ANM011/A011-20190220/A011-20190220-01',

    'Z:/Raphael_tests/mice_expdata/ANM012/A012-20190221/A012-20190221-01',
    'Z:/Raphael_tests/mice_expdata/ANM012/A012-20190223/A012-20190223-01',
    'Z:/Raphael_tests/mice_expdata/ANM012/A012-20190224/A012-20190224-01',

    'Z:/Raphael_tests/mice_expdata/ANM013/A013-20190504/A013-20190504-01',
    'Z:/Raphael_tests/mice_expdata/ANM013/A013-20190505/A013-20190505-01',

    'Z:/Raphael_tests/mice_expdata/ANM016/A016-20190531/A016-20190531-01',
    'Z:/Raphael_tests/mice_expdata/ANM016/A016-20190603/A016-20190603-01',

    'Z:/Raphael_tests/mice_expdata/ANM022/A022-20191105/A022-20191105-01',
    'Z:/Raphael_tests/mice_expdata/ANM022/A022-20191106/A022-20191106-01',
    'Z:/Raphael_tests/mice_expdata/ANM022/A022-20191107/A022-20191107-01',

    'Z:/Raphael_tests/mice_expdata/ANM023/A023-20191217/A023-20191217-01',
    'Z:/Raphael_tests/mice_expdata/ANM023/A023-20191220/A023-20191220-01',

    'Z:/Raphael_tests/mice_expdata/ANM024/A024-20200119/A024-20200119-01',
    'Z:/Raphael_tests/mice_expdata/ANM024/A024-20200122/A024-20200122-01',

    'Z:/Raphael_tests/mice_expdata/ANM025/A025-20200129/A025-20200129-01',
    'Z:/Raphael_tests/mice_expdata/ANM025/A025-20200207/A025-20200207-01',

    'Z:/Raphael_tests/mice_expdata/ANM028/A028-20200711/A028-20200711-01',
    'Z:/Raphael_tests/mice_expdata/ANM028/A028-20200712/A028-20200712-01',

    'Z:/Raphael_tests/mice_expdata/ANM029/A029-20200626/A029-20200626-01',
    'Z:/Raphael_tests/mice_expdata/ANM029/A029-20200627/A029-20200627-01',

    'Z:/Raphael_tests/mice_expdata/ANM030/A030-20200811/A030-20200811-01',
    'Z:/Raphael_tests/mice_expdata/ANM030/A030-20200812/A030-20200812-01',

    'Z:/Raphael_tests/mice_expdata/ANM030/A030-20200813/A030-20200813-01',
    'Z:/Raphael_tests/mice_expdata/ANM030/A030-20200815/A030-20200815-01',

    'Z:/Raphael_tests/mice_expdata/ANM037/A037-20201221/A037-20201221-01',
    'Z:/Raphael_tests/mice_expdata/ANM037/A037-20201222/A037-20201222-01',

    'Z:/Raphael_tests/mice_expdata/ANM039/A039-20210131/A039-20210131-01',
    'Z:/Raphael_tests/mice_expdata/ANM039/A039-20210201/A039-20210201-01',
    'Z:/Raphael_tests/mice_expdata/ANM039/A039-20210205/A039-20210205-01',

    'Z:/Raphael_tests/mice_expdata/ANM040/A040-20210307/A040-20210307-01',

    'Z:/Raphael_tests/mice_expdata/ANM040/A040-20210308/A040-20210308-01',
    'Z:/Raphael_tests/mice_expdata/ANM040/A040-20210309/A040-20210309-01',

    'Z:/Raphael_tests/mice_expdata/ANM041/A041-20210215/A041-20210215-01',
    'Z:/Raphael_tests/mice_expdata/ANM041/A041-20210216/A041-20210216-01',
    'Z:/Raphael_tests/mice_expdata/ANM041/A041-20210219/A041-20210219-01',

    'Z:/Raphael_tests/mice_expdata/ANM041/A041-20210221/A041-20210221-01',
    'Z:/Raphael_tests/mice_expdata/ANM042/A042-20210313/A042-20210313-01',
    'Z:/Raphael_tests/mice_expdata/ANM042/A042-20210316/A042-20210316-01',

    'Z:/Raphael_tests/mice_expdata/ANM042/A042-20210317/A042-20210317-01',
    'Z:/Raphael_tests/mice_expdata/ANM042/A042-20210319/A042-20210319-01',

    'Z:/Raphael_tests/mice_expdata/ANM044/A044-20210412/A044-20210412-01',

    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210421/A046-20210421-01',
    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210422/A046-20210422-01',
    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210423/A046-20210423-01',

    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210425/A046-20210425-01',
    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210426/A046-20210426-01',
    'Z:/Raphael_tests/mice_expdata/ANM046/A046-20210428/A046-20210428-02',

    'Z:/Raphael_tests/mice_expdata/ANM049/A049-20210722/A049-20210722-01',
    'Z:/Raphael_tests/mice_expdata/ANM049/A049-20210731/A049-20210731-01',

    'Z:/Raphael_tests/mice_expdata/ANM049/A049-20210801/A049-20210801-01',
    'Z:/Raphael_tests/mice_expdata/ANM049/A049-20210802/A049-20210802-01',
    'Z:/Raphael_tests/mice_expdata/ANM050/A050-20210828/A050-20210828-01',

    'Z:/Raphael_tests/mice_expdata/ANM050/A050-20210829/A050-20210829-02',
    'Z:/Raphael_tests/mice_expdata/ANM050/A050-20210831/A050-20210831-01',
    'Z:/Raphael_tests/mice_expdata/ANM050/A050-20210901/A050-20210901-01',

    'Z:/Raphael_tests/mice_expdata/ANM052/A052-20210915/A052-20210915-01',

    'Z:/Raphael_tests/mice_expdata/ANM052/A052-20210916/A052-20210916-01',
    'Z:/Raphael_tests/mice_expdata/ANM052/A052-20210917/A052-20210917-01',

    'Z:/Raphael_tests/mice_expdata/ANM053/A053-20211002/A053-20211002-01',
    'Z:/Raphael_tests/mice_expdata/ANM053/A053-20211004/A053-20211004-01',

    'Z:/Raphael_tests/mice_expdata/ANM053/A053-20211005/A053-20211005-01',
    'Z:/Raphael_tests/mice_expdata/ANM054/A054-20211008/A054-20211008-01',

    'Z:/Raphael_tests/mice_expdata/ANM054/A054-20211010/A054-20211010-01',
    'Z:/Raphael_tests/mice_expdata/ANM054/A054-20211013/A054-20211013-01',

    'Z:/Raphael_tests/mice_expdata/ANM054/A054-20211015/A054-20211015-01',

    'Z:/Raphael_tests/mice_expdata/ANM055/A055-20211028/A055-20211028-01',
    'Z:/Raphael_tests/mice_expdata/ANM055/A055-20211030/A055-20211030-01',
    'Z:/Raphael_tests/mice_expdata/ANM056/A056-20211110/A056-20211110-01',

    'Z:/Raphael_tests/mice_expdata/ANM056/A056-20211111/A056-20211111-01',
    'Z:/Raphael_tests/mice_expdata/ANM056/A056-20211114/A056-20211114-01',

    'Z:/Raphael_tests/mice_expdata/ANM056/A056-20211116/A056-20211116-01',
    'Z:/Raphael_tests/mice_expdata/ANM056/A056-20211117/A056-20211117-01',
    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211203/A057-20211203-01',

    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211204/A057-20211204-01',
    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211205/A057-20211205-01',
    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211208/A057-20211208-01',

    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211209/A057-20211209-01',
    'Z:/Raphael_tests/mice_expdata/ANM057/A057-20211210/A057-20211210-01',

    'Z:/Raphael_tests/mice_expdata/ANM062/A062-20220216/A062-20220216-01',
    'Z:/Raphael_tests/mice_expdata/ANM062/A062-20220218/A062-20220218-01',
    'Z:/Raphael_tests/mice_expdata/ANM063/A063-20220122/A063-20220122-01',

    'Z:/Raphael_tests/mice_expdata/ANM063/A063-20220124/A063-20220124-01',

    'Z:/Raphael_tests/mice_expdata/ANM063/A063-20220126/A063-20220126-01',
    'Z:/Raphael_tests/mice_expdata/ANM063/A063-20220127/A063-20220127-01',

    'Z:/Raphael_tests/mice_expdata/ANM063/A063-20220128/A063-20220128-01',
    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220228/A064-20220228-01',
    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220301/A064-20220301-01',

    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220302/A064-20220302-01',
    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220306/A064-20220306-01',
    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220307/A064-20220307-01',

    'Z:/Raphael_tests/mice_expdata/ANM064/A064-20220308/A064-20220308-01',
    'Z:/Raphael_tests/mice_expdata/ANM067/A067-20220501/A067-20220501-01',

    'Z:/Raphael_tests/mice_expdata/ANM067/A067-20220502/A067-20220502-01',

    'Z:/Raphael_tests/mice_expdata/ANM067/A067-20220503/A067-20220503-01',
    'Z:/Raphael_tests/mice_expdata/ANM067/A067-20220505/A067-20220505-01',
    'Z:/Raphael_tests/mice_expdata/ANM070/A070-20220623/A070-20220623-01',

    'Z:/Raphael_tests/mice_expdata/ANM070/A070-20220624/A070-20220624-01',
    'Z:/Raphael_tests/mice_expdata/ANM070/A070-20220625/A070-20220625-01',
    'Z:/Raphael_tests/mice_expdata/ANM072/A072-20220630/A072-20220630-01',

    'Z:/Raphael_tests/mice_expdata/ANM072/A072-20220701/A072-20220701-01',
    'Z:/Raphael_tests/mice_expdata/ANM074/A074-20220729/A074-20220729-01',
    'Z:/Raphael_tests/mice_expdata/ANM074/A074-20220730/A074-20220730-01',

    'Z:/Raphael_tests/mice_expdata/ANM074/A074-20220804/A074-20220804-01',
    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220820/A075-20220820-01',
    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220821/A075-20220821-01',

    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220822/A075-20220822-01',
    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220824/A075-20220824-01',
    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220825/A075-20220825-01',
    'Z:/Raphael_tests/mice_expdata/ANM075/A075-20220826/A075-20220826-01',

    'Z:/Raphael_tests/mice_expdata/ANM077/A077-20220924/A077-20220924-01',
    'Z:/Raphael_tests/mice_expdata/ANM081/A081-20221016/A081-20221016-01',
    'Z:/Raphael_tests/mice_expdata/ANM081/A081-20221017/A081-20221017-01',

    'Z:/Raphael_tests/mice_expdata/ANM081/A081-20221020/A081-20221020-01',
    'Z:/Raphael_tests/mice_expdata/ANM081/A081-20221023/A081-20221023-01',
    'Z:/Raphael_tests/mice_expdata/ANM082/A082-20221105/A082-20221105-01',

    'Z:/Raphael_tests/mice_expdata/ANM082/A082-20221106/A082-20221106-01',
    'Z:/Raphael_tests/mice_expdata/ANM082/A082-20221112/A082-20221112-01',

    'Z:/Raphael_tests/mice_expdata/ANM086/A086-20221213/A086-20221213-01',
    'Z:/Raphael_tests/mice_expdata/ANM086/A086-20221219/A086-20221219-01',

    'Z:/Raphael_tests/mice_expdata/ANM087/A087-20221219/A087-20221219-01',

    'Z:/Raphael_tests/mice_expdata/ANM087/A087-20221220/A087-20221220-01',

    'Z:/Raphael_tests/mice_expdata/ANM087/A087-20221224/A087-20221224-01',
    'Z:/Raphael_tests/mice_expdata/ANM092/A092-20230221/A092-20230221-01',
    'Z:/Raphael_tests/mice_expdata/ANM092/A092-20230222/A092-20230222-01',

    'Z:/Raphael_tests/mice_expdata/ANM095/A095-20230331/A095-20230331-01',
    'Z:/Raphael_tests/mice_expdata/ANM095/A095-20230401/A095-20230401-01',
    'Z:/Raphael_tests/mice_expdata/ANM096/A096-20230417/A096-20230417-01',

    'Z:/Raphael_tests/mice_expdata/ANM096/A096-20230421/A096-20230421-01',
    'Z:/Raphael_tests/mice_expdata/ANM096/A096-20230422/A096-20230422-01',

    'Z:/Raphael_tests/mice_expdata/ANM099/A099-20230503/A099-20230503-01',

    'Z:/Raphael_tests/mice_expdata/ANM099/A099-20230504/A099-20230504-01',
    'Z:/Raphael_tests/mice_expdata/ANM099/A099-20230507/A099-20230507-01',
    'Z:/Raphael_tests/mice_expdata/ANM099/A099-20230508/A099-20230508-01',

    'Z:/Raphael_tests/mice_expdata/ANM099/A099-20230509/A099-20230509-01',
    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230617/A102-20230617-01',
    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230618/A102-20230618-01',

    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230620/A102-20230620-01',
    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230621/A102-20230621-01',
    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230622/A102-20230622-01',

    'Z:/Raphael_tests/mice_expdata/ANM102/A102-20230623/A102-20230623-01',
    'Z:/Raphael_tests/mice_expdata/ANM104/A104-20230628/A104-20230628-01',
    'Z:/Raphael_tests/mice_expdata/ANM107/A107-20230731/A107-20230731-01',

    'Z:/Raphael_tests/mice_expdata/ANM107/A107-20230801/A107-20230801-01',
    'Z:/Raphael_tests/mice_expdata/ANM107/A107-20230802/A107-20230802-01',
    'Z:/Raphael_tests/mice_expdata/ANM107/A107-20230803/A107-20230803-01',

    'Z:/Raphael_tests/mice_expdata/ANM106/A106-20230802/A106-20230802-01',
    'Z:/Raphael_tests/mice_expdata/ANM106/A106-20230804/A106-20230804-01',
    'Z:/Raphael_tests/mice_expdata/ANM106/A106-20230806/A106-20230806-01',

    'Z:/Raphael_tests/mice_expdata/ANM106/A106-20230807/A106-20230807-01',

    'Z:/Raphael_tests/mice_expdata/ANM106/A106-20230808/A106-20230808-01',
    'Z:/Raphael_tests/mice_expdata/ANM109/A109-20230812/A109-20230812-01',

    'Z:/Raphael_tests/mice_expdata/ANM109/A109-20230813/A109-20230813-01',
    'Z:/Raphael_tests/mice_expdata/ANM109/A109-20230815/A109-20230815-01',
    'Z:/Raphael_tests/mice_expdata/ANM108/A108-20230818/A108-20230818-01',

    'Z:/Raphael_tests/mice_expdata/ANM110/A110-20230820/A110-20230820-01',
    'Z:/Raphael_tests/mice_expdata/ANM110/A110-20230822/A110-20230822-01',
    'Z:/Raphael_tests/mice_expdata/ANM110/A110-20230825/A110-20230825-01',

    'Z:/Raphael_tests/mice_expdata/ANM110/A110-20230826/A110-20230826-01',
    'Z:/Raphael_tests/mice_expdata/ANM110/A110-20230829/A110-20230829-01',

    'Z:/Raphael_tests/mice_expdata/ANM127/A127-20231219/A127-20231219-01',
    'Z:/Raphael_tests/mice_expdata/ANM127/A127-20231220/A127-20231220-01',
    'Z:/Raphael_tests/mice_expdata/ANM127/A127-20231223/A127-20231223-01',

    'Z:/Raphael_tests/mice_expdata/ANM126/A126-20240113/A126-20240113-01',
    'Z:/Raphael_tests/mice_expdata/ANM126/A126-20240115/A126-20240115-01',

    'Z:/Raphael_tests/mice_expdata/ANM126/A126-20240118/A126-20240118-01',
    'Z:/Raphael_tests/mice_expdata/ANM126/A126-20240119/A126-20240119-01',
    'Z:/Raphael_tests/mice_expdata/ANM126/A126-20240120/A126-20240120-01',
    'Z:/Raphael_tests/mice_expdata/ANM129/A129-20240120/A129-20240120-01',
    'Z:/Raphael_tests/mice_expdata/ANM129/A129-20240121/A129-20240121-01',

    'Z:/Raphael_tests/mice_expdata/ANM129/A129-20240122/A129-20240122-01',
    'Z:/Raphael_tests/mice_expdata/ANM129/A129-20240124/A129-20240124-01',
    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250513/A167-20250513-01',

    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250514/A167-20250514-01',
    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250515/A167-20250515-01',
    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250517/A167-20250517-01',

    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250518/A167-20250518-01',
    'Z:/Raphael_tests/mice_expdata/ANM167/A167-20250519/A167-20250519-01',
    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250520/A170-20250520-01',

    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250521/A170-20250521-01',
    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250522/A170-20250522-01',
    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250524/A170-20250524-01',

    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250525/A170-20250525-01',
    'Z:/Raphael_tests/mice_expdata/ANM170/A170-20250526/A170-20250526-01',
    'Z:/Raphael_tests/mice_expdata/ANM169/A169-20250524/A169-20250524-01',

    'Z:/Raphael_tests/mice_expdata/ANM169/A169-20250525/A169-20250525-01',
    'Z:/Raphael_tests/mice_expdata/ANM169/A169-20250526/A169-20250526-01',
    'Z:/Raphael_tests/mice_expdata/ANM165/A165-20250601/A165-20250601-01',

    'Z:/Raphael_tests/mice_expdata/ANM165/A165-20250602/A165-20250602-01',
    'Z:/Raphael_tests/mice_expdata/ANM165/A165-20250603/A165-20250603-01',
    'Z:/Raphael_tests/mice_expdata/ANM165/A165-20250605/A165-20250605-01'
    ]

# all maze_sess are resolved in-script with iteration now
