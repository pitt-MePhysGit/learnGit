function [ args ] = prep_ampOddClick_function( exp, electrodeLayout, SR, doSpikes )

% mod TT 20150722:
% changes sampling rate for llp to 500 Hz and
% expanded window from -150 to 750 ms.
% no baseline correction: needs to be done in R

% new version based on pre_TRF_function


% ============================

%args.baseDir         = '/data/';
%args.baseDir         = '/Volumes/DATA/EEG/';
%args.baseDir         = '/Volumes/192.168.0.1/';
%args.baseDir         = '/Volumes/rawData/EEG/';
%args.baseDir          = find_baseDir();

args.baseDir         = '';
if exist('/Volumes/DATA/EEG/')
    args.baseDir         = '/Volumes/DATA/EEG/';
end

if exist('/Volumes/rawData/EEG/')
    args.baseDir         = '/Volumes/rawData/EEG/';
end



args.electrodeLayout = electrodeLayout;

args.IntanComp       = 'PlanetExpress';
args.monkeyLogicComp = 'Leela';
args.directory       = 'ampOddClick';

args.soundDelay      = 0.001;  % sound delay in seconds, for EK2 equal to 1ms

%args.baseFile        = '';
%args.localDir        = '/data/EEGLab/';
%args.localDir        = '/Volumes/DATA/EEG/EEGLab/';
args.localDir        = '/Volumes/Drobo5D3/EEG/EEGLab/';
%args.localDir        = find_localDir();

args.stimartefact    = 0;
args.triggerartefact    = 1;

args.exp             = exp;

args.epoch      = [-.150 .750];
if SR > 9999
    args.epoch = [-.010 .250];
end

tmp          = strsplit(args.exp{1},'_');
animal       = tmp(1);
date         = tmp(2);
args.animal  = animal{1};
args.date    = date{1};


%startFromScratch = 1;
forceRedo        = 1;
a                = exist( [args.localDir '/' args.directory '/' args.exp{1} '/amp.mat'],'file');
startFromScratch = ~( a>0 & ~forceRedo );

% -----------------------------------------
[args.baseDir args.animal '/' args.animal '_' args.date  '/' args.exp{1} '/*.rhd']
start = 1;
for i = 1:length(exp)
    [args.baseDir args.animal '/' args.animal '_' args.date  '/' args.exp{i} '/*.rhd']
    d = dir(           [args.baseDir args.animal '/' args.animal '_' args.date  '/' args.exp{i} '/*.rhd']);
    [dirList{start:(start+length(d)-1),1}] = deal(d.name);
    start = length(dirList)+1;
end
args.dirList = dirList;%(1:9);
args.dirList
%%
[electrodeNumber, analogInputNumber] = getElectrodeLayout(args);
args.analogInputNumber               = analogInputNumber;
args.electrodeNumber                 = electrodeNumber;

%% specify desired properties of read-in data
args.SRamp = SR;   % desired sampling rate amplifier
args.SRdig = -1;    % desired sampling rate digital data;
args.SRana = -1;    % desired sampling rate analog data
args.SRmua = SR;   % desired smpaling rate for mua
% if minus 1, use highest possible sampling rate

args.Namp = length(args.electrodeNumber);
args.Nana = length(args.analogInputNumber);
args.Ndig = 1;

mkdir( [args.localDir '/' args.directory '/' args.exp{1}]  );
%SLsave( [args.localDir '/' args.directory '/' args.exp{1} '/args_' num2str(SR)], 'args' );
SLsave( [args.localDir '/' args.directory '/' args.exp{1} '/'], ['args_' num2str(args.SRamp)], args );

%% check if there are any potential spike-channels
args.spkChan = [];

if sum(args.electrodeLayout ~= 'v')>0 & doSpikes
    suargs = args;
    suargs.electrodeLayout = args.electrodeLayout( find(args.electrodeLayout ~= 'e') );
    [args.spkChan,~]       = getElectrodeLayout(suargs);
end

if startFromScratch
    %% ==========================================================
    [ amp ana dig mua tamp tana tdig tmua impedance ] = read_Intan_function( args.dirList, args);
end

if ~startFromScratch
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/amp.mat']);
    amp = object;
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/ana.mat']);
    ana = object;
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/dig.mat']);
    dig = object;
    
    doMua = 0;
    if exist( [args.localDir '/' args.directory '/' args.exp{1} '/mua.mat'] )
       doMua = 1; 
    end
    
    if doMua
        load( [args.localDir '/' args.directory '/' args.exp{1} '/mua.mat']);
        mua = object;
    end
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/tamp.mat']);
    tamp = object;
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/tana.mat']);
    tana = object;
    
    load( [args.localDir '/' args.directory '/' args.exp{1} '/tdig.mat']);
    tdig = object;
    
    if doMua
        load( [args.localDir '/' args.directory '/' args.exp{1} '/tmua.mat']);
        tmua = object;
    end
    impedance = 0;
end

%% ===========================================================
[tone, trial, valSoundStr] = analyzeEvents_ampOddClick(dig, ana, tdig, tdig, args, impedance);

SLsave( [args.localDir '/' args.directory '/' args.exp{1} '/'], ['args_' num2str(args.SRamp)], args );

%% ============================================================

raw_amp       = pop_importdata('setname','raw_amp','data',amp, 'dataformat','matlab','nbchan',args.Namp+args.Nana+args.Ndig,'xmin',0,'srate',args.SRamp,'ref','recording');
raw_amp.times = tamp;

raw_mua       = pop_importdata('setname','raw_mua','data',mua, 'dataformat','matlab','nbchan',args.Namp,'xmin',0,'srate',args.SRmua,'ref','recording');
raw_mua.times = tmua;


% add event information
[raw_amp, eventnumbers] = pop_importevent(raw_amp, 'event',[args.localDir '/' args.directory '/' args.exp{1} '/events.txt'],  'fields',{'type','latency' }, 'append', 'no', 'timeunit', 1 );
[raw_mua, eventnumbers] = pop_importevent(raw_mua, 'event',[args.localDir '/' args.directory '/' args.exp{1} '/events.txt'],  'fields',{'type','latency' }, 'append', 'no', 'timeunit', 1 );


% extract epochs
raw                = pop_epoch( raw_amp, valSoundStr, args.epoch);
raw.setname        = ['raw_' num2str(SR)];

tmp                = raw.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:); % don't rereference the analog and digital input channels
raw.data           = reref(raw.data,1, 'keepref', 'on');
raw.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:) = tmp;               % restore unreferenced analog input channels

% save data
disc               = pop_saveset(raw    , 'filename',raw.setname,'filepath', [args.localDir '/' args.directory '/' args.exp{1} '/'], 'savemode','onefile');

%
mua                = pop_epoch( raw_mua, valSoundStr, args.epoch);
mua.setname        = ['MUA_' num2str(SR)];
disc               = pop_saveset(mua    , 'filename',mua.setname,'filepath', [args.localDir '/' args.directory '/' args.exp{1} '/'], 'savemode','onefile');



%% filtered 'long latency' potentials
hicutoff  = 40;   % 80
locutoff  = 2;
filtorder = ceil(0.256*args.SRamp);  % 256 milliseconds long

['transition bandwidth = ' num2str(args.SRamp* 5.5/filtorder) 'Hz']
['filter length        = ' num2str(1000/args.SRamp*     filtorder) 'ms']
m = pop_firwsord('blackman', args.SRamp, 20)
%

[raw_amp_filt, com, b] = pop_firws(raw_amp, 'fcutoff',[locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');

% allow filtering of LeftAudio
raw_amp_filt.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:) = raw_amp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:);

llp = pop_epoch( raw_amp_filt, valSoundStr, args.epoch);
%llp = lineNoiseNotchEEG( llp, [60 120], 34:57,0 );

clear raw_amp_filt

llp.setname        = 'LongLatencyPotentials';

tmp                = llp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:); % don't rereference the analog input channels
llp.data           = reref(llp.data,1, 'keepref', 'on');
llp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:) = tmp;               % restore unreferenced analog input channels

thisFileName = ['LLP_' num2str(SR)];
disc     = pop_saveset(llp,  'filename', thisFileName, 'filepath', [args.localDir '/' args.directory '/' args.exp{1} '/'], 'savemode','onefile');

%% -------------------------------
%% filtered 'FFR' potentials
locutoff  = 100;
filtorder = ceil(0.032*args.SRamp);  % 16 milliseconds long

['transition bandwidth = ' num2str(args.SRamp* 5.5/filtorder) 'Hz']
['filter length        = ' num2str(1000/args.SRamp*     filtorder) 'ms']
m = pop_firwsord('blackman', args.SRamp, 20)
%

[raw_amp_filt, com, b] = pop_firws(raw_amp, 'fcutoff',[locutoff], 'forder',filtorder,'ftype','highpass');

% allow filtering of LeftAudio
raw_amp_filt.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:) = raw_amp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:);
llp = pop_epoch( raw_amp_filt, valSoundStr, args.epoch);

clear raw_amp_filt


llp.setname        = ['FFR' num2str(SR)];
tmp                = llp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:); % don't rereference the analog input channels
llp.data           = reref(llp.data,1, 'keepref', 'on');
llp.data( (args.Namp+1):(args.Namp+args.Nana+args.Ndig),:,:) = tmp;               % restore unreferenced analog input channels

thisFileName = ['FFR_' num2str(SR)];
disc     = pop_saveset(llp,  'filename', thisFileName, 'filepath', [args.localDir '/' args.directory '/' args.exp{1} '/'], 'savemode','onefile');

%%
%disc = ampOddClick2R(exp{1}, 0); % export data
'done!'

end


