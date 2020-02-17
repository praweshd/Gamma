%% Prawesh Dahal
% July 3, 2019
%Prawesh Dahal (last verified Sept 3, 2019)
% Grand Gamma Detection (works for both high and low)
% - detect Gamma on all channels, 
% - store the res into a struct 
% - save the res file 

close all
clear all
clc

% A = [56	30	82	68	8	53	123	109	33
% 60	26	81	67	10	98	103	115	37
% 61	29	78	69	91	97	104	116	34
% 63	25	77	70	92	59	126	118	38
% 62	24	79	86	12	57	125	117	41
% 58	20	80	85	14	55	102	107	45
% 54	16	64	71	17	95	101	108	48
% 50	13	15	72	19	96	42	119	40
% 47	9	18	83	90	4	44	120	127
% 43	5	22	84	89	2	46	106	111
% 39	1	23	74	65	0	49	105	112
% 36	21	27	73	66	6	99	121	114
% 31	3	32	75	87	94	100	122	113
% 35	7	28	76	88	93	51	124	110];
% 
% B = reshape(A',1,126);
% B = B+1;
% 
% NG_2019_v4.map = B;
% NG_2019_v4.dim = [9, 14];


%createbadchannels([setdiff(1:128,good_CH+1)]-1)
 
%% If LFP has ground noise, filter
% 
% 
% lfp_file = dir('*.lfp');
% f_in = lfp_file.name; 
% extention = 'fil';
% fc = [275]; 
% filter_type = 'low'; 
% Rs = 1250; 
% CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
% dat2fil(f_in,extention,fc,filter_type,Rs,CH_N)

 
%%
lfp_file = dir('*.lfp');
CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
[~, fbasename, ~] = fileparts(lfp_file.name);

gamma_file = strcat(fbasename, '_highGAM_RES'); 

%%
% resfile = dir('*_autowakestates*');
% load(resfile.name); 
% wakestates = find(states==1); 
% %%
% clear states
% resfile = dir('*-states*');
% load(resfile.name); 
% 
% states(wakestates) =1; 
% 
% %%
tic

parfor i = 1 : 128
    
[gamma] = Gamma_Detect(CH_N, i); 

gammas(i).res = gamma; 

% gamma_file = strcat(fbasename, '_', num2str(i), 'gamma_wav'); 
% save(gamma_file, 'gamma');
gamma_events = strcat(fbasename,'_', num2str(i), '.gam.evt');
channelID = i-1; 
SaveRippleEvents(gamma_events, gamma, channelID);
 
mkdir EventFiles
movefile([strcat(fbasename,'_', num2str(i), '.gam.evt')], strcat(pwd,'/','EventFiles')) 


end 

toc
save(gamma_file, 'gammas');

%%
 