folderPath = 'Z:\Dinghao\Behav\DataAnalysis\';

%% active licking
% optogenetics Dbh (stArch) (exp 2)

activeLickOptoDbhPath = [...
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD007\A007-20210903\';... % 007 1st opto, stim at reward onset
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD007\A007-20210906\';... % 007 2st opto, stim at reward onset
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211123\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211124\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211129\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211129\'; % s1
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211129\'; % s2
%     ...
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211130\'; % s2 file corrupted but saw an effect, s1 no recov control
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211130\'; % s1 noisy
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211130\'; % s2, likely underpowered
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211201\'; % s1 noisy
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211201\'; %s2 only, s1 noisy, s1 had an effect though, so s2 likely valid
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211202\'; % likely underpowered
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211202\'; % noisy
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211203\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211206\'; % why decreased power? most likely underpowered
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211208\'; %s2 only, s1 noisy, power low though
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211209\'; %underpowered
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211209\'; % noisy
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211210\'; %underpowered
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211210\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211213\'; % underpowered
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211214\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211215\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20211216\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD014\A014-20211216\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220113\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220114\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220117\'; % maybe underpowered?
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220118\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220119\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD015\A015-20220120\'; % increased power for the last time
%     ...
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220216\'; % noisy
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220217\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220218\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220222\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220223\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220224\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220228\'; %incredibly noisy
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220301\'; %noisy
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220302\'; % underpowered, raising power next day
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220303\'; % likely still underpowered 
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220304\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD021\A021-20220308\'; % final note: in general the effect for this animal was not strong, histology showed very weak expression
%     ...
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220322\'; % no recov cont
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220323\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220324\'; % was going to have an s2 but animal performing really badly with grm
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220325\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220328\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220329\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220330\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220331\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220504\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220505\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220506\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD023\A023-20220509\';
%     ...
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220511\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220513\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220518\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220519\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220520\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220523\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220524\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220526\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220601\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220602\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220603\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220606\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220607\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD025\A025-20220608\';
%     ...
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220805\';
% %     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220808\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220809\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220810\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220811\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220812\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220816\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220824\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD031\A031-20220825\';
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230224\';  % first animal using halo, 4 wks
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230309\';  % 7 wks 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230310\'; 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230323\';  % 040 increased numLicksBefRew and totStopLenTRun, but some noise at the lickport
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230324\';  % no visible effects
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230325\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230328\';  % increased licks bef reward, only during stim trials
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230329\';  % noisy
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230331\';  % noisy
    'Z:\Dinghao\Behav\DataAnalysis\ANMD053\A053-20230403\'; 
    ...
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230424\';  % control too noise 
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230425\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230426\';  % full trial 6 s
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230501\';  % full trial recovery really bad
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230502\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230503\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230505\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230509\';  % first pupil tracking
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230510\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230511\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230515\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230516\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230517\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230518\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230519\';  % controls
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230522\';  % controls 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230523\';  % controls
    'Z:\Dinghao\Behav\DataAnalysis\ANMD057\A057-20230525\';  % controls 
];  

ALRecSessionsOptoDbh = [... % recording session (file no.)
% %     {[1 1 1]};  % 40 CONT., 100 STIM., 100 Post-CONT.
% %     {[1 1 1]};
% %     {[1 1 1]};
% %     {[1 1 1]};
% %     {[1 1 1]};
% %     {[1 1 1]};
% %     {[1 1 1]};
%     ...
% %     {[1 1 1]};
% %     {[1 1 2]};
% %     {[1 1 1]};
% %     {[1 1 1]};
%     {[2 2 2]};
% %     {[2 2 2]};
% %     {[1 1 1]};
%     {[1 1 2 3 3 3]};
% %     {[1 1 1]};
% %     {[2 2 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1]};
% %     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2 3]};
%     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
%     {[1 1 1]};
%     {[1 1 1]};
%     {[1 1 1 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     ...
% %     {[1 1 1 2]};
% %     {[1 1 1 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1]};
%     {[1 1 1 2 3 3 3]};
%     {[1 1 1 2 3 3 3]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2]};
% %     {[1 1 1 2 2 2 3]};
%     ...
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 3 3 3]};
%     {[1 1 1 2]};
%     {[1 1 1 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1]};
%     {[1 1 1]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     ...
%     {[1 1 1]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1]};
% %     {[1 1 1]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 3 3 3]};
%     {[1 1 1]};
%     {[1 1 1 2 3 3 3]};
%     {[1 1 1 2 2 2]};
%     ...
%     {[1 1 1]};
% %     {[1 1 1 2 2 2]};
%     {[1 1 1]};
%     {[1 1 1]};
%     {[1 1 1]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1]};
%     {[1 1 1]};
%     {[1 1 1]};
    ...
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
%     {[1 1 1]};
%     {[1 1 1]};
    {[1 1 1 2 2 2]};
    ...
%     {[1 1 1]};
%     {[1 1 2]};
    {[1 1 1 2 2 2]};
%     {[1 1 1]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
];

ALMazeTypeOptoDbh = [... % stimulation type (6-8: 2-4 but 3 seconds, 9: uncoupled 020, 10: uncoupled 030, 11: uncoupled 040)
% %     {[0 4 0]};
% %     {[0 4 0]};
% %     {[0 2 0]};
% %     {[0 4 0]};
% %     {[0 4 0]};
% %     {[0 4 0]};
% %     {[0 4 0]};
%     ...
% %     {[0 4 0]};
% %     {[0 4 0]};
% %     {[0 4 10]};
% %     {[0 4 0]};
%     {[0 3 0]};
% %     {[0 4 10]};
% %     {[0 4 0]};
%     {[0 4 0 0 4 0]};
% %     {[0 4 10]};
% %     {[0 3 0]};
% %     {[0 4 0 0 2 0]};
% %     {[0 4 0 0 2 0]};
% %     {[0 4 0 0 5 0]};
% %     {[0 4 0]};
% %     {[0 4 10 0 3 0]};
%     {[0 4 0 0 4 10 0]};
%     {[0 4 0 0 4 0]};
% %     {[0 4 0 0 4 0]};
%     {[0 4 10]};
%     {[0 4 10]};
%     {[0 4 10 0]};
%     {[0 3 0 0 4 10]};
%     {[0 3 0 0 4 10]};
% %     {[0 2 0 0 4 10]};
%     {[0 3 0 0 4 10]};
%     ...
% %     {[0 4 0 0]};
% %     {[0 4 0 0]};
% %     {[0 2 0 0 4 0]};
% %     {[0 3 0]};
%     {[0 4 10 0 0 3 0]};
%     {[0 4 10 0 0 2 11]};
% %     {[0 4 10 0 3 0]};
% %     {[0 4 10]};
% %     {[0 3 0 0 4 10]};
% %     {[0 3 0 0 4 10]};
% %     {[0 4 10 0 4 10]};
% %     {[0 4 10 0 2 0 0]};
%     ...
%     {[0 4 10 0 4 10]};
%     {[0 4 10 0 0 2 0]};
%     {[0 2 11 0]};
%     {[0 2 11 0]};
%     {[0 3 0 0 2 0]};
%     {[0 4 0 0 4 0]};
%     {[0 4 0 0 4 0]};
%     {[0 4 0]};
%     {[0 4 0]};
%     {[0 4 0 0 3 0]};
%     {[0 4 0 0 2 0]};
%     {[0 4 0 0 3 0]};
%     ...
%     {[0 4 0]};
%     {[0 4 0 0 3 0]};
%     {[0 4 0]};
% %     {[0 4 0]};
%     {[0 4 0 0 2 0]};
%     {[0 4 0 0 2 0]};
%     {[0 4 0 0 3 0]};
%     {[0 4 0 0 2 0]};
%     {[0 4 10 0]};
%     {[0 2 0 0 4 0]};
%     {[0 2 0 0 0 3 0]};
%     {[0 2 0]};
%     {[0 2 0 0 0 3 0]};
%     {[0 4 0 0 4 0]};
%     ...
%     {[0 4 0]};
% %     {[0 4 0 0 4 0]};
%     {[0 4 0]};
%     {[0 4 0]};
%     {[0 4 0]};
%     {[0 2 0 0 3 0]};
%     {[0 4 0]};
%     {[0 4 0]};
%     {[0 4 0]};
    ...
    {[0 8 0 0 6 0 0 7 0]};
    {[0 8 0 0 7 0]};
    {[0 8 0]};
    {[0 8 0 0 8 0 0 7 0]};
    {[0 8 0 0 7 0 0 9 0]};
    {[0 8 0 0 7 0]};
    {[0 8 0 0 7 0]};
%     {[0 4 0]};
%     {[0 4 0]};
    {[0 8 0 0 7 0]};
    ...
%     {[0 4 0]};
%     {[0 4 0]};
    {[0 4 0 0 4 0]};
%     {[0 4 0]};
    {[0 4 0 0 2 0]};
    {[0 4 0 0 4 0]};
    {[0 4 0 0 3 0 0 2 0]};
    {[0 4 0 0 2 0]};
    {[0 4 0 0 2 0]};
    {[0 4 0 0 3 0]};
    {[0 2 0 0 2 0 0 4 0]};
    {[0 2 0 0 4 0]};
    {[0 3 0 0 2 0 0 4 0]};
    {[0 3 0 0 2 0 0 9 0]};
    {[0 9 0 0 9 0]};  % uncoupled
    {[0 9 0 0 9 0 0 9 0]};
    {[0 11 0 0 9 0 0 9 0]};
    {[0 9 0 0 11 0]};
];

ALMazeSessionOptoDbh = [ ... % subsession number within one session (or file)
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
%     ...
% %     {[1 2 3]};
% %     {[1 2 1]};
% %     {[1 2 3]};
% %     {[1 2 3]};
%     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
%     {[1 2 1 1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3 1]};
%     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
%     {[1 2 3 1]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     ...
% %     {[1 2 3 1]};
% %     {[1 2 3 1]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3]};
%     {[1 2 3 1 1 2 3]};
%     {[1 2 3 1 1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3]};
% %     {[1 2 3 1 2 3 1]};
%     ...
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 1 2 3]};
%     {[1 2 3 1]};
%     {[1 2 3 1]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     ...
%     {[1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3]};
% %     {[1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3 1 1 2 3]};
%     {[1 2 3 1 2 3]};
%     ...
% %     {[1 2 3]};
% %     {[1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
    ...
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3]};
    {[1 2 3 1 2 3]};
    ...
%     {[1 2 3]};
%     {[1 2 3]};
    {[1 2 3 1 2 3]};
%     {[1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
];

ALOptoStimPulseWDbh = [ ... % the length of the stimulation pulse
% %     {[0 5 0]};
% %     {[0 6 0]};
% %     {[0 5 0]};
% %     {[0 5 0]};
% %     {[0 5 0]};
% %     {[0 5 0]};
% %     {[0 3 0]};
%     ...
% %     {[0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0]};
%     {[0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0]};
%     {[0 3 0 0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0 0]};
%     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
%     {[0 3 0 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     ...
% %     {[0 5 0 0]};
% %     {[0 5 0 0]};
% %     {[0 5 0 0 5 0]};
% %     {[0 3 0]};
%     {[0 3 0 0 0 3 0]};
%     {[0 3 0 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0]};
% %     {[0 3 0 0 3 0 0]};
%     ...
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 0 3 0]};
%     {[0 3 0 0]};
%     {[0 3 0 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     ...
%     {[0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0]};
% %     {[0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0 0 0 3 0]};
%     {[0 3 0 0 3 0]};
%     ...
%     {[0 3 0]};
% %     {[0 3 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
%     {[0 3 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
    ...
    {[0 3 0 0 3 0 0 3 0]};
    {[0 3 0 0 3 0]};
    {[0 3 0]};
    {[0 3 0 0 3 0 0 3 0]};
    {[0 3 0 0 3 0 0 3 0]};
    {[0 3 0 0 3 0]};
    {[0 3 0 0 3 0]};
%     {[0 3 0]};
%     {[0 3 0]};
    {[0 3 0 0 3 0]};
    ...
%     {[0 3 0]};
%     {[0 3 0]};
    {[0 6 0 0 6 0]};
%     {[0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};    
    {[0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};
    {[0 5 0 0 5 0 0 5 0]};
    {[0 5 0 0 5 0]};
];

ALOptoPowerDbh = [... % power of optogenetics
% %     3.25; % 5.8
% %     4.5;
% %     4.5;
% %     4.5;
% %     4;
% %     2.5;
% %     2.5;
%     ...
% %     2.5;
% %     2.5;
% %     2.5;
% %     2.5;
%     2.0;
% %     2.0;
% %     2.5;
%     2.8;
% %     2.5;
% %     2.5;
% %     3.5;
% %     4.0;
% %     4.0;
% %     4.5;
% %     3.5;
%     4.2;
%     5.0;
% %     5.0;
%     5.0;
%     4.2;
%     5.0;
%     5.0;
%     5.0;
% %     4.5;
%     5.5;
%     ...
% %     1.0; %starts from v low for this animal to avoid overstimulation
% %     1.0;
% %     1.2;
% %     1.5;
%     1.5;
%     1.5;
% %     1.5;
% %     1.8;
% %     1.8;
% %     2.2;
% %     3.0;
% %     4.5;
%     ...
%     1.0;
%     1.5;
%     1.5;
%     1.5;
%     2.0;
%     3.0;
%     3.2;
%     3.5;
%     4.5;
%     4.5;
%     5.5;
%     5.5;
%     ...
%     1.0;
%     1.2;
%     1.5;
% %     2.0;
%     2.0;
%     2.0;
%     2.5;
%     2.5;
%     3.0;
%     3.0;
%     4.0;
%     5.0;
%     5.0;
%     5.0;
%     ...
%     1.0;
% %     1.2;
%     2.0;
%     2.0;
%     2.0;
%     2.4;
%     3.5;
%     4.5;
%     5.0;
    ...
    1.4;
    2.1;
    2.1;
    3.2;
    3.2;
    3.2;
    3.2;
%     3.2;
%     3.5;
    3.5;
    ...
%     1.5;
%     1.5;
    3.0;
%     3.0;
    3.0;
    3.0;
    3.0;
    4.0;
    4.0;
    4.0;
    4.5;
    4.5;
    4.5;
    5.0;
    5.0;
    5.0;
    5.0;
    5.0;
];

% Sulpiride (exp 1) (cond 1)

activeLickSulPath = [...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210901\';... % 005 1st SUL session (1mM, 200nl, pH = 5.94)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD006\A006-20210907\';... % 006 1st SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210908\';... % 005 2nd SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD006\A006-20210909\';... % 006 2nd SUL session (1mM, 200nl, pH = 5.96)
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211006\';... % 009 1st SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211015\';... % 009 2nd SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211021\';... % 009 3rd SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211027\';... % 010 1st SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211029\';... % 010 2nd SUL session (1mM, 200nl, pH = 5.96)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211101\';... % 011 1st SUL session
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211108\';... % 011 2nd SUL session
];

ALRecSessionsSul = [...
    {[1 2 3 4 5]}; % 1-3 omitted, 4 as 1: CONT., INJ. omitted, 7 as 2: SUL30mins, 3: SUL1hr, 4: SUL2hr, 5: SUL3hr
    {[1 2 3 4 5]}; % 1 CONT., INJ. omitted, 5 as 2: SUL30mins, 3: SUL1hr, 4: SUL2hr, 5: SUL3hr
    {[1 2 3 4]}; % 1-4 bad CONT., 5 as 1: CONT., INJ. omitted, 7 as 2: SUL1hr, 3: SUL2hr, 4: SUL3hr
    {[1 2 3]}; % 1-3 bad belt condition, 4 as 1: CONT., INJ. omitted, 6 belt stuck, 7/8 as 2: SUL1hr, 9 as 3: SUL2hr
    ...
    {[1 2 3 4]}; 
    {[1 2 3 4]};
    {[1 2 3 4]};
    {[1 2 3 4]};
    {[1 2 3 4]};
    {[1 2 3 4]};
    {[1 2 3 4]};
];

ALMazeTypeSul = [...
    {[0 2 1 1 1]};
    {[0 2 1 1 1]};
    {[0 2 1 1]};
    {[0 2 1]};
    ...
    {[0 2 1 1]};
    {[0 2 1 1]};
    {[0 2 1 1]};
    {[0 2 1 1]};
    {[0 2 1 1]};
    {[0 2 1 1]};
    {[0 2 1 1]};
];

ALMazeSessionSul = [ ...
    {[1 1 1 1 1]};
    {[1 1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1]};
    ...
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
];

ALOptoStimPulseWSul = [...
    {[0 0 0 0 0]};
    {[0 0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0]};
    ...
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
];

ALOptoPowerSul = [...
    0.0;
    0.0;
    0.0;
    0.0;
    ...
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
    0.0;
];

% SCH (exp 3)

activeLickSCHPath = [...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD006\A006-20210913\';... % 006 1st SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD006\A006-20210915\';... % 006 2nd SCH session
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210921\';... % 005 1st SCH session bad run
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210923\';... % 005 2nd SCH session bad run
    ...
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211004\';... % 009 1st SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211011\';... % 009 2nd SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211013\';... % 009 3rd SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211014\';... % 010 1st SCH session
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211019\';... % 010 2nd SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211022\';... % 010 3rd SCH session
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211027\';... % 011 1st SCH session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211103\';... % 011 2nd SCH session
];

ALRecSessionsSCH = [...
    {[1 2 3 4]}; % from 3 as 1: CONT., INJ. omitted, 5 as 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    {[1 2 3 4]}; % from 6 as 1: CONT., INJ. omitted, 8 as 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
%     {[1 2 3 4]}; % from 3 as 1: CONT., INJ. omitted, 5 as 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
%     {[1 2 3 4]}; % from 3 as 1: CONT., INJ. omitted, 5 as 2: SCH0.5hr, 3: SHC1hr, 4: SCH2hr
    ...
%     {[1 2 3 4]}; % 1: CONT., 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    {[1 2 3 4]}; % 1: CONT., 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    {[1 2 3 4]}; % 1: CONT., 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    {[1 2 3 4]}; % 1: CONT., 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
%     {[1 2 3 4]}; % 1: CONT., 2: SCH0.5hr, 3: SCH1hr, 4: SCH2hr
    {[1 2 3 4]};
%     {[1 2 3 4]};
    {[1 2 3 4]};
];

ALMazeTypeSCH = [...
    {[0 1 1 0]};
    {[0 1 1 0]};
%     {[0 1 1 0]};
%     {[0 1 1 0]};
    ...
%     {[0 1 1 0]};
    {[0 1 2 0]};
    {[0 1 1 0]};
    {[0 1 2 0]};
%     {[0 1 1 0]};
    {[0 1 1 0]};
%     {[0 1 1 0]};
    {[0 1 1 0]};
];

ALMazeSessionSCH = [ ...
    {[1 1 1 1]};
    {[1 1 1 1]};
%     {[1 1 1 1]};
%     {[1 1 1 1]};
    ...
%     {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
%     {[1 1 1 1]};
    {[1 1 1 1]};
%     {[1 1 1 1]};
    {[1 1 1 1]};
];

ALOptoStimPulseWSCH = [...
    {[0 0 0 0]};
    {[0 0 0 0]};
%     {[0 0 0 0]};
%     {[0 0 0 0]};
    ...
%     {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
%     {[0 0 0 0]};
    {[0 0 0 0]};
%     {[0 0 0 0]};
    {[0 0 0 0]};
];

ALOptoPowerSCH = [...
    0.0;
    0.0;
%     0.0;
%     0.0;
    ...
%     0.0;
    0.0;
    0.0;
    0.0;
%     0.0;
    0.0;
%     0.0;
    0.0;
];

% Saline (exp 4)

activeLickSalPath = [...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD006\A006-20210917\';... % 006 saline session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210925\';... % 005 saline session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD005\A005-20210927\';... % 005 2nd saline session
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211008\';... % 009 saline session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD009\A009-20211019\';... % 009 2nd saline session
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211025\';... % 010 saline session
    'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211105\';... % 011 saline session
];

ALRecSessionsSal = [...
    {[1 2 3 4]}; % s2 as 1: CONT., 3&4: INJ., 5&6: nervous, belt stuck, 7 as 2: SAL30min, 3: SAL1hr, 4: SAL2hr
    {[1 2 3 4]}; % s2 as 1: CONT., INJ. using ANMD999, 2: SAL30min, 3: SAL1hr, 4: SAL2hr
    {[1 2 3 4]}; % s2 as 1: CONT., INJ. using ANMD999, 2: SAL30min, 3: SAL1hr, 4: SAL2hr
    ...
    {[1 2 3 4]}; % s1: CONT., INJ. omitted, 2: SAL30min, 3: SAL1hr, 4: SAL2hr
    {[1 2 3 4]}; % 1: CONT., INJ. omitted, 2: SAL30min, 3: SAL1hr, 4: SAL2hr
%     {[1 2 3 4]};
    {[1 2 3]};
];

ALMazeTypeSal = [...
    {[0 1 0 0]};
    {[0 1 0 0]};
    {[0 1 0 0]};
    ...
    {[0 1 0 0]};
    {[0 1 0 0]};
%     {[0 1 0 0]};
    {[0 1 0]};
];

ALMazeSessionSal = [ ...
    {[1 1 1 1]};
    {[1 1 1 1]};
    {[1 1 1 1]};
    ...
    {[1 1 1 1]};
    {[1 1 1 1]};
%     {[1 1 1 1]};
    {[1 1 1]};
];

ALOptoStimPulseWSal = [...
    {[0 0 0 0]};
    {[0 0 0 0]};
    {[0 0 0 0]};
    ...
    {[0 0 0 0]};
    {[0 0 0 0]};
%     {[0 0 0 0]};
    {[0 0 0]};
];

ALOptoPowerSal = [...
    0.0;
    0.0;
    0.0;
    ...
    0.0;
    0.0;
%     0.0;
    0.0;
];

% VTA Dat Opto (cond 5)

activeLickOptoDatPath = [...
	'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221013\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221014\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221017\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221019\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221020\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD040\A040-20221021\';
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD043\A043-20221108\';
];

ALRecSessionsOptoDat = [... % recording session (file no.)
	{[1 1 1]};
    {[1 1 1]};
    {[1 1 1]};
    {[1 1 1]};
    {[1 1 1]};
    {[1 1 1]};
    ...
    {[1 1 1 2 2 2 3 3 3]};
];

ALMazeTypeOptoDat = [... % stimulation type (0 -- control; other numbers depend on the protocol; 10 -- recov subsess that did not see recov)
    {[0 2 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    ...
    {[0 3 0 0 4 0 0 3 0]};
];

ALMazeSessionOptoDat = [ ... % subsession number within one session (or file)
	{[1 2 3]};
    {[1 2 3]};
    {[1 2 3]};
    {[1 2 3]};
    {[1 2 3]};
    {[1 2 3]};
    ...
    {[1 2 3 1 2 3 1 2 3]};
];

ALOptoStimPulseWDat = [ ... % the length of the stimulation pulse
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    {[0 3 0]};
    ...
    {[0 3 0 0 3 0 0 3 0]};
];

ALOptoPowerDat = [... % power of optogenetics AT TIPS, OUTSIDE BRAIN
    1.5;
    1.5;
    2.0;
    2.5;
    2.5;
    2.5;
    ...
    1.0;
];

% Dbh Opto activation (cond 6)

activeLickOptoDbhActPath = [...
	'Z:\Dinghao\Behav\DataAnalysis\ANMD052\A052-20230327\';  % both 030 and 040 less licks on stim
    'Z:\Dinghao\Behav\DataAnalysis\ANMD052\A052-20230328\';  % 030 a bit less licks on stim (non-sig)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD052\A052-20230329\';  % 030 more licks on stim, 040 less licks on stim
    'Z:\Dinghao\Behav\DataAnalysis\ANMD052\A052-20230330\';  % 030 more licks on stim, 040 less licks on stim
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230503\';  % so far it seems that med 1st five licks increase
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230504\';  % first 020, med lick distance increased as well 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230505\';  % second session has some speed change
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230509\';  % a bit noisy, no recovery 020
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230510\';  % noisy
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230512\'; 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230518\';   % last session motivation issue at recovery, exclude
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230519\';   % controls 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD059\A059-20230523\';
    ...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230619\';   % low power (0.2 mW)
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230620\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230621\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230622\';  % first 010, unfortunately behaviour is not good
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230623\';  % too much noise
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230624\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230628\';  % rewrapped to isolate noise 
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230629\';  % the noise returned 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230630\';  % 2nd session too noisy
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230703\';  % 2nd session speed drop + noise 
    'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230705\';
%     'Z:\Dinghao\Behav\DataAnalysis\ANMD061\A061-20230706\';  % speed is not controlled
];

ALRecSessionsOptoDbhAct = [... % recording session (file no.)
	{[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    ...
    {[1 1 1 2 2 2]};
    {[1 1 1]};
    {[1 1 1 2 2 2 3 3 3]};
%     {[1 1 1]};
%     {[1 1 1 2 2 2]};
    {[1 1 1]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1]};
    ...
    {[1 1 1]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2 3 3 3]};
    {[1 1 1 2 2 2]};
    {[1 1 1]};
%     {[1 1 1 2 2 2]};
    {[1 1 1]};
    {[1 1 1]};
    {[1 1 1 2 2 2]};
%     {[1 1 1 2 2 2]};
];

ALMazeTypeOptoDbhAct = [... % stimulation type (6-8: 2-4 but 3 seconds, 9: uncoupled 020, 10: uncoupled 030, 11: uncoupled 040)
    {[0 3 0 0 4 0]};
    {[0 3 0 0 3 0]};
    {[0 3 0 0 4 0]};
    {[0 3 0 0 4 0]};
    ...
    {[0 4 0 0 4 0]};
    {[0 2 0]};
    {[0 4 0 0 3 0 0 3 0]};
%     {[0 2 0]};
%     {[0 4 0 0 2 0]};
    {[0 2 0]};
    {[0 4 0 0 2 0]};
    {[0 11 0 0 11 0]};
    {[0 3 0]};
    ...
    {[0 2 0]};
    {[0 2 0 0 2 0]};
    {[0 2 0 0 2 0]};
%     {[0 2 0 0 1 0]};
%     {[0 2 0 0 1 0 0 2 0]};
    {[0 2 0 0 1 0]};
    {[0 2 0]};
%     {[0 2 0 0 2 0]};
    {[0 2 0]};
    {[0 2 0]};
    {[0 2 0 0 1 0]};
%     {[0 2 0 0 2 0]};
];

ALMazeSessionOptoDbhAct = [ ... % subsession number within one session (or file)
	{[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    ...
    {[1 2 3 1 2 3]};
    {[1 2 3]};
    {[1 2 3 1 2 3 1 2 3]};
%     {[1 2 3]};
%     {[1 2 3 1 2 3]};
    {[1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3]};
    ...
    {[1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3]};
%     {[1 2 3 1 2 3]};
    {[1 2 3]};
    {[1 2 3]};
    {[1 2 3 1 2 3]};
%     {[1 2 3 1 2 3]};
];

ALOptoStimPulseWDbhAct = [ ... % the length of the stimulation pulse
    {[0 1 0 0 1 0]};
    {[0 1 0 0 1 0]};
    {[0 1 0 0 1 0]};
    {[0 1 0 0 1 0]};
    ...
    {[0 1 0 0 1 0]};
    {[0 1 0]};
    {[0 1 0 0 1 0 0 1 0]};
%     {[0 1 0]};
%     {[0 1 0 0 1 0]};
    {[0 1 0]};
    {[0 1 0 0 1 0]};
    {[0 1 0 0 1 0]};
    {[0 1 0]};
    ...
    {[0 .5 0]};
    {[0 .5 0 0 .5 0]};
    {[0 .5 0 0 .5 0]};
%     {[0 .5 0 0 .5 0]};
%     {[0 .5 0 0 .5 0 0 .5 0]};
    {[0 .5 0 0 .5 0]};
    {[0 .5 0]};
%     {[0 .5 0 0 .5 0]};
    {[0 .5 0]};
    {[0 .5 0]};
    {[0 .5 0 0 .5 0]};
%     {[0 .5 0 0 .5 0]};
];

ALOptoPowerDbhAct = [... % power of optogenetics AT TIPS, OUTSIDE BRAIN
    5.5;
    5.5;
    5.5;
    5.5;
    ...
    4.0;
    4.0;
    4.0;
%     4.0;
%     4.0;
    4.4;
    4.8;
    4.8;
    4.8;
    ...
    0.2;
    1.4;
    1.4;
%     1.4;
%     1.4;
    1.4;
    1.4;
%     1.4;
    1.9;
    1.9;
    2.4;
%     2.4;
];

% Dbh Opto activation (EPHYS - NEURONEXUS) (cond 7)

activeLickOptoDbhActEphysPath = [...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD062\A062-20230624\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD062\A062-20230626\';
    'Z:\Dinghao\Behav\DataAnalysis\ANMD062\A062-20230627\';
];

ALRecSessionsOptoDbhActEphys = [... % recording session (file no.)
	{[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
    {[1 1 1 2 2 2]};
];

ALMazeTypeOptoDbhActEphys = [... % stimulation type (6-8: 2-4 but 3 seconds, 9: uncoupled 020, 10: uncoupled 030, 11: uncoupled 040)
    {[0 2 0 0 1 0]};
    {[0 2 0 0 1 0]};
    {[0 2 0 0 2 0]};
];

ALMazeSessionOptoDbhActEphys = [ ... % subsession number within one session (or file)
	{[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
    {[1 2 3 1 2 3]};
];

ALOptoStimPulseWDbhActEphys = [ ... % the length of the stimulation pulse
    {[0 .5 0 0 .5 0]};
    {[0 .5 0 0 .5 0]};
    {[0 .5 0 0 .5 0]};
];

ALOptoPowerDbhActEphys = [... % power of optogenetics AT TIPS, OUTSIDE BRAIN
    0.5;
    0.5;
    0.5
];

% training practice

PracticePath = [...
    'Z:\Dinghao\Behav\DataAnalysis\ANMD011\A011-20211207\';... % methiothepin mesylate salt inj.
    'Z:\Dinghao\Behav\DataAnalysis\ANMD010\A010-20211208\';... % methiothepin mesylate salt inj.
];

ALRecSessionsOptoPractice = [... % recording session (file no.)
    {[1 2 3 4]};
    {[1 2 3 4]};
];

ALOptoMazeTypePractice = [... % stimulation type (0 -- control; other numbers depend on the protocol)
    {[0 1 1 1]};
    {[0 1 1 1]};
];

ALOptoMazeSessionPractice = [ ... % subsession number within one session (or file)
    {[1 1 1 1]};
    {[1 1 1 1]};
];

ALOptoStimPulseWPractice = [ ... % the length of the stimulation pulse
    {[0 0 0 0]};
    {[0 0 0 0]};
];

ALOptoPowerPractice = [... % power of optogenetics
    0.0;
    0.0;
];