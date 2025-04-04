function [ disc ] = prep_RS_flex_rockey_function_old( exp )

 % mod TT 20150722:
 % changes sampling rate for llp to 500 Hz and 
 % expanded window from -150 to 750 ms.
 % no baseline correction: needs to be done in R

overnight       = 0; % if 1, do filtering on full-bandwidth data. 
                     % if 0, downsample first, then filter

% define paths 
%baseDir         = '~/Dropbox/Nibbler/';
baseDir         = '/Volumes/rawData/EEG/';

IntanComp       = 'PlanetExpress';
monkeyLogicComp = 'Leela';
directory       = 'ampOddClick';

soundDelay      = 0.001;  % sound delay in seconds, for EK2 equal to 1ms

%% -------------------------------------
% setting ketamine to one sets short filter setting as default
ketamine         = 1;
shortFilter      = 0;  % use short Filter variants?

%% ========================

% === define paths for the one test-data set  
%baseDir         = '~/Documents';
%IntanComp       = '';
%monkeyLogicComp = 'Leela';
%directory       = 'testTiming192kHz';
% ===================================================

baseFile        = '';
localDir        = '/Volumes/rawData/EEG/EEGLab/';

tmp          = strsplit(exp{1},'_');
animal       = tmp(1);

startFromScratch = 1;

xflag            = 0;
if shortFilter==1
    xflag = 2;
end

doLLP   = 1; % long-latency potentials
doRAW   = 1; % raw down-sampled
doRaw5k = 1; % 5kS raw data
doMLP   = 1; % mid-latency potentials
doMBP   = 1; % mid-brain potentials
doBSP   = 1; % brain-stem potentials

% -----------------------------------------

start = 1;
for i = 1:length(exp)
    
    %d = dir(           [baseDir '/' IntanComp '/' directory  '/' exp{i} '/*.rhd']);
    d = dir(           [baseDir animal{1} '/' exp{i} '/*.rhd']);
    [dirList{start:(start+length(d)-1),1}] = deal(d.name);
    start = length(dirList)+1;

end

%%
%channelFile = '~/Dropbox/eCog/rockey/positionList_ISI_surgeryAdjusted.txt';
%w = textread(channelFile);

if strcmp(animal,'Rockey')
    
    %electrode.number = [8 5 32 -1 14 19 22    6 1 29 11 13 18   4 31 28 27 15 20   23 24];
    electrode.number = [8 5 32 -1 15 19 22    6 1 29 11 13 18   4 31 28 27 14 20   23 24];  % swapped input for channel 5 and 18. tt20151123
    electrode.name   = {'1x8_0x1_z','2x8_0x1_z','3x8_0x1_z','4x8_0x1_z','5x8_0x1_z','6x8_0x1_z','7x8_0x1_z',...
        '2x8_1x2_l','3x8_1x2_l','4x8_1x3_l','4x8_2x3_l','5x8_1x2_l','6x8_1x2_l',...
        '2x8_1x2_r','3x8_1x2_r','4x8_1x3_r','4x8_2x3_r','5x8_1x2_r','6x8_1x2_r', 'EOG_l','EOG_r'};
    
    electrode.ml      = [0 0 0 0 0 0 0    -1.3 -1.78 -1.32 -2.4 -1.8 -1.5    1.3 1.78 1.32 2.4 1.8 1.5    -1.65 1.65];
    electrode.ap      = [-2.8 -1.95 -0.84 0.3 1.35 2.61 3.7    -2.22 -1.17 -0.04 0.24 1.33 2.53    -2.22 -1.17 -0.04 0.24 1.33 2.53    5.0 5.0];
    electrode.dv      = [ 2.3 3.25 3.85 4.09 4.13 4 3.17   2.9 3.5 4.1 3.45 3.85 3.85    2.9 3.5 4.1 3.45 3.85 3.85    2.75 2.75];
    
    analogInput.number = [5 3 4];
end

% 28  -> 26 ->
if strcmp(animal,'Jesse')
    
    %electrode.number = [8 -1 3 27 14 19 22    6 2 30 9 13 18   5 1 28 26 15 20    7 32 10 11 12 17    4 31 29 25 16 21   23 24];
    electrode.number = [8 -1 3 28 14 19 22    6 2 30 9 13 18   5 1 27 26 15 20    7 32 10 11 12 17    4 31 29 25 16 21   23 24];  % swapped input for channel 4 and 16. tt20151123
    
    electrode.name   = {'1x8_0x1_z','2x8_0x1_z','3x8_0x1_z','4x8_0x1_z','5x8_0x1_z','6x8_0x1_z','7x8_0x1_z',...
        '2x8_1x2_l','3x8_1x2_l','4x8_1x3_l','4x8_2x3_l','5x8_1x2_l','6x8_1x2_l',...
        '2x8_1x2_r','3x8_1x2_r','4x8_1x3_r','4x8_2x3_r','5x8_1x2_r','6x8_1x2_r',...
        '2x8_2x2_l','3x8_2x2_l','4x8_3x3_l','4x8_4x3_l','5x8_2x2_l','6x8_2x2_l',...
        '2x8_2x2_r','3x8_2x2_r','4x8_3x3_r','4x8_4x3_r','5x8_2x2_r','6x8_2x2_r',...
        'EOG_l','EOG_r'};
    
    % needs to be fixed
    electrode.ml      = [0 0 0 0 0 0 0                         -1.3 -1.78 -1.32 -2.4 -1.8 -1.5      1.3 1.78 1.32 2.4 1.8 1.5         0 0 0 0 0 0    0 0 0 0 0 0   -1.65 1.65  ];
    electrode.ap      = [-2.8 -1.95 -0.84 0.3 1.35 2.61 3.7    -2.22 -1.17 -0.04 0.24 1.33 2.53    -2.22 -1.17 -0.04 0.24 1.33 2.53   0 0 0 0 0 0    0 0 0 0 0 0    5.0 5.0  ];
    electrode.dv      = [ 2.3 3.25 3.85 4.09 4.13 4 3.17        2.9 3.5 4.1 3.45 3.85 3.85          2.9 3.5 4.1 3.45 3.85 3.85        0 0 0 0 0 0    0 0 0 0 0 0    2.75 2.75     ];
    
    analogInput.number = [5 3 4];
end

if strcmp(animal,'Walter')
    
    %electrode.number = [2 -1 16 15 17 18 19     9 10 12 11 13 14   25 24 23 22 21 20    3 4 6 5 7 8   31 30 29 28 27 26    1 32];
    electrode.number = [2 -1 16 15 17 18 19     9 10 12 11 13 14   25 24 23 22 21 20    3 4 5 6 7 8   31 30 29 28 27 26    1 32];   % swapped input for channels 22 and 23. tt20151123
    
    electrode.name   = {'1x8_0x1_z','2x8_0x1_z','3x8_0x1_z','4x8_0x1_z','5x8_0x1_z','6x8_0x1_z','7x8_0x1_z',...
        '2x8_1x2_l','3x8_1x2_l','4x8_1x3_l','4x8_2x3_l','5x8_1x2_l','6x8_1x2_l',...
        '2x8_1x2_r','3x8_1x2_r','4x8_1x3_r','4x8_2x3_r','5x8_1x2_r','6x8_1x2_r',...
        '2x8_2x2_l','3x8_2x2_l','4x8_3x3_l','4x8_4x3_l','5x8_2x2_l','6x8_2x2_l',...
        '2x8_2x2_r','3x8_2x2_r','4x8_3x3_r','4x8_4x3_r','5x8_2x2_r','6x8_2x2_r',...
        'COG_l','COG_r'};
    
    % needs to be fixed
    electrode.ml      = [0 0 0 0 0 0 0                         -1.3 -1.78 -1.32 -2.4 -1.8 -1.5      1.3 1.78 1.32 2.4 1.8 1.5         0 0 0 0 0 0    0 0 0 0 0 0   -1.65 1.65  ];
    electrode.ap      = [-2.8 -1.95 -0.84 0.3 1.35 2.61 3.7    -2.22 -1.17 -0.04 0.24 1.33 2.53    -2.22 -1.17 -0.04 0.24 1.33 2.53   0 0 0 0 0 0    0 0 0 0 0 0    5.0 5.0  ];
    electrode.dv      = [ 2.3 3.25 3.85 4.09 4.13 4 3.17        2.9 3.5 4.1 3.45 3.85 3.85          2.9 3.5 4.1 3.45 3.85 3.85        0 0 0 0 0 0    0 0 0 0 0 0    2.75 2.75     ];
    
    analogInput.number = [5 3 4];
end



if strcmp(animal,'Sam')
    
    %electrode.number = [32 -1 1 15 16 17 18     10 9 12 11 13 14   24 23 22 21 20 19    2 3 4 5 6 7   31 30 28 27 29 26    8 25];
    electrode.number = [32 -1 1 15 16 17 18     10 9 12 11 13 14   24 23 22 21 20 19    2 3 5 4 6 7   31 30 28 27 29 26    8 25];   % swapped input for channels 22 and 23. tt20151123
    
    electrode.name   = {'1x8_0x1_z','2x8_0x1_z','3x8_0x1_z','4x8_0x1_z','5x8_0x1_z','6x8_0x1_z','7x8_0x1_z',...
        '2x8_1x2_l','3x8_1x2_l','4x8_1x3_l','4x8_2x3_l','5x8_1x2_l','6x8_1x2_l',...
        '2x8_1x2_r','3x8_1x2_r','4x8_1x3_r','4x8_2x3_r','5x8_1x2_r','6x8_1x2_r',...
        '2x8_2x2_l','3x8_2x2_l','4x8_3x3_l','4x8_4x3_l','5x8_2x2_l','6x8_2x2_l',...
        '2x8_2x2_r','3x8_2x2_r','4x8_3x3_r','4x8_4x3_r','5x8_2x2_r','6x8_2x2_r',...
        'COG_l','COG_r'};
    % needs to be fixed
    electrode.ml      = [0 0 0 0 0 0 0                         -1.3 -1.78 -1.32 -2.4 -1.8 -1.5      1.3 1.78 1.32 2.4 1.8 1.5         0 0 0 0 0 0    0 0 0 0 0 0   -1.65 1.65  ];
    electrode.ap      = [-2.8 -1.95 -0.84 0.3 1.35 2.61 3.7    -2.22 -1.17 -0.04 0.24 1.33 2.53    -2.22 -1.17 -0.04 0.24 1.33 2.53   0 0 0 0 0 0    0 0 0 0 0 0    5.0 5.0  ];
    electrode.dv      = [ 2.3 3.25 3.85 4.09 4.13 4 3.17        2.9 3.5 4.1 3.45 3.85 3.85          2.9 3.5 4.1 3.45 3.85 3.85        0 0 0 0 0 0    0 0 0 0 0 0    2.75 2.75     ];
    
    analogInput.number = [5 3 4];
end


Namp = length(find(electrode.number));
Nana = length(analogInput.number);

clear chanlocs
chan = struct('labels',[], 'type',[], 'theta',[], 'radius',[], 'X',[], 'Y',[], 'Z',[],'sph_theta',[],'sph_phi',[], 'sph_radius',[], 'urchan',[], 'ref',[]);
for c = 1:length(electrode.number)
    chanlocs(c) = chan;
    chanlocs(c).labels = electrode.name(c);
    chanlocs(c).X      = electrode.ml(c);
    chanlocs(c).Y      = electrode.ap(c);
    chanlocs(c).Z      = electrode.dv(c);
    chanlocs(c).urchan = c;
end


%electrode.number = [-1  1 2 3 4 5    6 7 8 9 10   11 12 13 14 15   16 17 18 19 20   21 22 23 24 25   26 27 28 29 30   31 32];
%%
if startFromScratch
    %% ====================== start reading the Intan data
    %my_read_Intan_RHD2000_fcn( [baseDir '/' IntanComp '/' directory '/' exp{1} '/'] ,dirList{1})
    my_read_Intan_RHD2000_fcn( [baseDir animal{1} '/' exp{1} '/'] ,dirList{1})

    processAmp = exist( 'amplifier_data'    ,'var'); % process amplifier data if present
    processDig = exist( 'board_dig_in_data' ,'var'); % process digital in if present
    processAna = 1; %exist( 'board_adc_data'   ,'var'); % process analog in if present
    
    all_t_amplifier    = zeros( (size(t_amplifier       )).*[1 length(dirList)+.1] );
    if processAmp
        all_amplifier_data = zeros( [length(analogInput.number)+length(electrode.number) size(amplifier_data,2)].*[1 length(dirList)+.1] );
    end
    
    if processDig
        all_digin_data     = zeros( size(board_dig_in_data ).*[1 length(dirList)+.1] );
    end
    
    if processAna
        all_anain_data     = zeros( size(board_adc_data    ).*[1 length(dirList)+.1] );
    end
    
    startInd = 1;
    for cind = 1:length(dirList)
        startInd
        
        %my_read_Intan_RHD2000_fcn( [baseDir '/' IntanComp '/' directory '/' exp{1} '/'] ,dirList{cind})
        my_read_Intan_RHD2000_fcn( [baseDir animal{1} '/' exp{1} '/'] ,dirList{cind})

        stopInd = startInd + length(t_amplifier) - 1;
        
        addData = 1;
        max(abs(diff(t_amplifier)))
        if max(abs(diff(t_amplifier)))>.001
            %keyboard
            addData = 0;
            
            % in some cases there is just a temporary blip
            % find how many violations there are:
            
            if length( find(abs(diff(t_amplifier))>.001 ))== 2  % just one skip
                
                skipInd = find(abs(diff(t_amplifier))>.001 );
                
                diffT = diff(t_amplifier);
                diffT(skipInd) = median( diff(t_amplifier) );
                newT = t_amplifier(1) + cumsum( [0 diffT] );
                
                dtt = newT - t_amplifier;
                t_amplifier( skipInd(1):(skipInd(2)+1) ) = newT( skipInd(1):(skipInd(2)+1) );
                
                % test if problem was resolved
                skipInd = find(abs(diff(t_amplifier))>.001 );
                if isempty(skipInd)
                    addData = 1;
                end
                
                if ~isempty(skipInd) % take a closer look at what is going on
                    keyboard
                end
                
            end
            
           
        end
        
        if addData
            all_t_amplifier(1, startInd:stopInd )         = t_amplifier;         %+ max(all_t_amplifier)+t_amplifier(2) ;
            if processAmp
                all_amplifier_data( electrode.number>0,  startInd:stopInd)      = amplifier_data(electrode.number(electrode.number>0),:);
                all_amplifier_data((Namp+1):(Namp+Nana) ,startInd:stopInd)      = board_adc_data(analogInput.number,:);
            end
            
            if processDig
                all_digin_data(    :,startInd:stopInd)        = board_dig_in_data;
            end
            
            if processAna
                all_anain_data(    :,startInd:stopInd)        = board_adc_data;
            end
            
            startInd = stopInd + 1;
        end
    end
    
    if processAmp
        all_amplifier_data      = all_amplifier_data(:,1:stopInd);
        %all_amplifier_data(8,:) = 0;
    end
    
    if processDig
        all_digin_data         = all_digin_data(:,1:stopInd);
    end
    
    if processAna
        all_anain_data         = all_anain_data(:,1:stopInd);
    end
    
    all_t_amplifier        = all_t_amplifier(1,1:stopInd);
    
    tres_amp = median(diff(all_t_amplifier));
    if min(diff(all_t_amplifier))+.000000001 < max(diff(all_t_amplifier))
        'warning: tres might have changed during recording'
    end
    
    % fix bug in timing when combining two data sets
    tdiff = diff(all_t_amplifier);
    all_t_amplifier_real = all_t_amplifier;
    all_t_amplifier_fake = all_t_amplifier;
    if min(tdiff)<0
        % old correction assumed no interruptions, necessary to extract epochs
        tdiff( tdiff<0 ) = tres_amp;
        all_t_amplifier_fake  = cumsum( [0 tdiff] );
        
        tdiff = diff(all_t_amplifier);
        tdiff( tdiff<0 ) = 30; % new version: assumes 30 seconds between
        all_t_amplifier_real  = cumsum( [0 tdiff] );
    end
    
    
    num_amplifier_channels = size(all_amplifier_data,1);
    
    if processAna
        num_anain_channels     = size(all_anain_data,1);
    end
    
    %% filtering on full data set is very time-consuming.
    %  run overnight
    
    %if overnight
    %    % band-pass filter the full-bandwidth raw data
    %    locutoff  = 1 ;  % 120 HZ highpass filter
    %    hicutoff  = 1250;  %
    %    filtorder = 16502;
    %    EEG = struct('data',all_amplifier_data, 'srate',1/tres_amp, 'trials',1,'setname','raw','nbchan',size(all_amplifier_data,1),'event',[],...
    %                  'pnts',size(all_amplifier_data,2),'times',all_t_amplifier,
    %                  'xmin',min(all_t_amplifier),
    %                  'xmax',max(all_t_amplifier));
    %
    %    [filt, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    %    all_amplifier_data = filt.data;
    %end
    
    
    
    if( doLLP | doRAW)
        %% downsample from 5000 to 1000 Hz, only necessary for LLP potentials
        a =  1;
        %b =  5; % 5;
        b = 10; % mod TT 20150722
        
        if b==1  % no downsampling, just rename
            down_amplifier_data   = all_amplifier_data;
            clear all_amplifier_data
            
            tres_down             = tres_amp;
            down_t_amplifier_real = all_t_amplifier_real;
            down_t_amplifier_fake = all_t_amplifier_fake;
        end
        
        
        if b>1
            'Downsampling...'
            tst                 = resample(all_amplifier_data(1,:), a,b);
            down_amplifier_data = zeros( num_amplifier_channels, size(tst,2) ); % add one row for the reference electrode
            %down_anain_data     = zeros( num_anain_channels,       size(tst,2) );
            %down_t_amplifier    = resample(all_t_amplifier,a,b);
            tres_down           = tres_amp * b/a;
            
            % old version
            %down_t_amplifier    =min(all_t_amplifier):tres_down:max(all_t_amplifier);%
            
            %new version
            down_t_amplifier_real    = all_t_amplifier_real( round(b/2): b : length(all_t_amplifier_real) );
            down_t_amplifier_fake    = all_t_amplifier_fake( round(b/2): b : length(all_t_amplifier_fake) );
            
            
            if abs( median(diff(down_t_amplifier_fake)) - tres_down ) > 0.0000001
                'check tres after resampling'
            end
            
            print_increment = 10;
            percent_done = print_increment;
            for i=1:num_amplifier_channels
                down_amplifier_data(i,:) = resample(all_amplifier_data(i,:), a, b);
                %if(i < num_anain_channels+1)
                %down_anain_data(i,:) = resample(all_anain_data(i,:), a, b);
                %end
                
                fraction_done = 100 * (i / num_amplifier_channels);
                if (fraction_done >= percent_done)
                    fprintf(1, '%d%% done...\n', percent_done);
                    percent_done = percent_done + print_increment;
                end
                
            end
        end
        
        % save the data for import in eeglab
        if ~exist( [localDir '/' directory '/' exp{1}], 'dir')
            mkdir( [localDir '/' directory '/' exp{1}] )
        end
    end
    
    %save( [localDir '/' directory '/' exp{1} '/amp.mat'], 'down_amplifier_data', '-v7.3')
    %if processAna
    %    save([localDir '/' directory '/' exp{1} '/ana.mat'], 'all_anain_data','-v7.3')
    %end
    
    
    % keep all_amplifier_data
    % clear all_amplifier_data raw_amp
    
    %% ===============================================
    %   read the .bhv data using the bhv_read function
    
    %setpref('MonkeyLogic','Directories',struct('BaseDirectory','~/Dropbox/MonkeyLogic','RunTimeDirectory','~/Dropbox/MonkeyLogic/runtime/','ExperimentDirectory',[baseDir '/' monkeyLogicComp '/' exp '/'] ));
    %trials = bhv_read( [baseDir monkeyLogicComp '/' directory '/' exp '.bhv']  );
    
    %% ==================================
   
    % Find all events and event
    eventStrobe    = all_digin_data(1:8,:); %vsyncInd+1);
    eventCodeCnt   = eventStrobe' * [1;2;4;8; 16;32;64;128];
    
    vsyncIndTrig = diff(all_digin_data(16,:)) > .1;
    vsyncInd     = abs(diff([eventCodeCnt(1) eventCodeCnt']))>0.1;
    
    eventTime_real   = all_t_amplifier_real(vsyncInd);
    eventTime_fake   = all_t_amplifier_fake(vsyncInd);
    eventCode   = eventCodeCnt(vsyncInd);
    
    Ind45               = find(eventCode == 116);
    trialOnsetInd       = Ind45( find( diff([0 eventTime_real(Ind45)]) > 0.1));
    if(trialOnsetInd(1)>5) % first trial Onset trigger was missed
        trialOnsetInd = [1; trialOnsetInd];
    end
    
    trialOnsetTime_real = all_t_amplifier_real(vsyncInd(trialOnsetInd));
    trialOnsetTime_fake = all_t_amplifier_fake(vsyncInd(trialOnsetInd));
    
    
    % calculate trial number for each event index
    trialInd                = zeros(1,length(eventCode));
    trialInd(trialOnsetInd) = 1;
    trialNrAll              = cumsum(trialInd);
    
    % find which tones are being used in this paradigm
    for i = 1:23
        tbl(i) = length(find(eventCode==i));
    end
    valSounds = find(tbl>1);
    ['Valid Sounds: ', int2str(valSounds) ]
    soundStr = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
    valSoundStr = soundStr(valSounds);
    
    % find all tone onsets
    toneOnsetInd = find(ismember(eventCode, valSounds));
    
    % remove the trial info triggers from the list of potential toneOnsets
    tmpOnsetInd      = toneOnsetInd(find(toneOnsetInd>1)); % toneOnsetInd should never be the first index, but it can happen if ML does not get turned off after intan crash
    trialInfoTrigger = find(eventCode(tmpOnsetInd-1)==99);
    toneOnsetInd     = toneOnsetInd(  setdiff(1:length(toneOnsetInd), trialInfoTrigger )  );
    
    % there are some random events that happen between trials, but are
    % not associated with an actual tone-presentation. These can be
    % recognized by a lack of a 32 marker immediately following it
    keep = eventCode(toneOnsetInd+1)==32;
    eventCode(toneOnsetInd(keep==0)) = 999;
    toneOnsetInd = toneOnsetInd(keep);
    ['Found ', int2str(length(find(keep==0))),' non-terminated sounds']
    
    
    % trial number for tone onsets only
    tone.trialNr                 = trialNrAll(toneOnsetInd);
    if(tone.trialNr(1)==0)
        % in case the trial onset trigger was not recorded on the first trial
        tone.trialNr = tone.trialNr + 1;
    end
    
    tone.time = eventTime_real(toneOnsetInd);
    tone.ISI  = [-1 diff(tone.time)];
    
    toneInd                = zeros(1,length(eventCode));
    toneInd(toneOnsetInd)  = 1;
    toneCount              = cumsum(toneInd);
    
    % find sequence number for each tone onset
    seqNrAll   = zeros(1, length( eventCode ));
    seqCntr    = 0;
    for i = 1:length(eventCode)
        seqCntr = seqCntr + toneInd(i);
        seqNrAll(i) = seqCntr;
        
        if trialInd(i) == 1
            seqCntr = 0;
        end
    end
    tone.seqNr = seqNrAll(toneOnsetInd);
    
    % calculate tone height for each tone onset
    tone.pitch  = eventCode(toneOnsetInd)';
    tone.deltaPitch = [-99 diff(tone.pitch)];
    
    tone.ampInd     = mod(eventCode(toneOnsetInd),11)';
    
    % calcluate trial type of each tone: 1: oddball; 2: controll mean; 3: controll abs 1; 4: controll abs 2;
    tone.trial   = ones(size(tone.pitch)); %trialTypeMap(tmp);
    
    % calcluate the tone type: 0: repeat; 1:novel
    tone.type   = ones(size(tone.pitch));
    
    % ---------------------- trial
    devInd = find(tone.type==1);
    trial.deviantPosition                       = ones(1, max(tone.trialNr) );
    trial.deviantPitch                          = -1 * ones(1, max(tone.trialNr) );
    
    if(max(tone.seqNr) < 14) % makes only sense in the trial paradigms
        trial.deviantPosition                       = 14 * ones(1, max(tone.trialNr) );
        trial.deviantPosition(tone.trialNr(devInd)) = tone.seqNr(devInd);
        
        trial.deviantPitch(tone.trialNr(devInd)) = tone.pitch(devInd);
    end
    
    trial.onsetTone = find(tone.seqNr==1);
    trial.onsetTime = tone.time(trial.onsetTone);
    trial.type      = -99*ones(1,max(tone.trialNr));
    
    trialOnsetIndVal = trialOnsetInd(find(trialOnsetInd+3<length(eventCode)) );
    tmp               = eventCode(trialOnsetIndVal+3);  % these values should all be 99
    tmp2              = eventCode(trialOnsetIndVal+4);
    badInd            = find(~(tmp2>60 & tmp2<70));
    if length(badInd)>0
        tmp(badInd) = 0;
    end
  
    % standard skip is 4 in order to get to the trial information
    jump             = 4*ones(length(trialOnsetIndVal),1 );
    % in some cases there are extra or missing triggers
    missingInd       = find(tmp~=99);
    
    if( length(find(tmp~=99))>0 )
        %keyboard
        for mind = 1:length(missingInd) 
           strt  = trialOnsetIndVal(missingInd(mind)); 
           jmp   = find(  eventCode( (strt+1):(strt+7))>63 & eventCode( (strt+1):(strt+7))<69 ); 
           if length(jmp)==1
                jump(missingInd(mind)) = jmp;
           else
               jump(missingInd(mind)) = NaN;
           end
           
        end
        
        jump
        %trialOnsetInd
        %'no trialType information available in the digital codes'
        
    end
    
    if length(find(~isnan(jump)))>0
        trial.type(find(~isnan(jump)))      = eventCode(trialOnsetIndVal(find(~isnan(jump)))+jump(find(~isnan(jump)))) - 64;
    end
    
    % calculate the sequence number of the deviant for each trial; use 14 as
    % default
    tone.deviantPosition = trial.deviantPosition(tone.trialNr);
    tone.relPosition     = tone.seqNr - tone.deviantPosition;
    tone.deviantPitch    = trial.deviantPitch(tone.trialNr);
    
    %% =========
    % for all stimulus Onset events, look for the microphone onset or left audio trigger
    % If found, use actual onset time rather than predicted onset time
    
    eventTimeCorrected_real = eventTime_real;
    eventTimeCorrected_fake = eventTime_fake;
    
    %chanInd = 2;  % use chanInd 4 for mic in booth; use chanInd=5 after channel 2 broke
    chanInd = 5;
    if size(all_anain_data,1)==1
        chanInd = 1;
    end
    
    %thrsh = [0.02 0.015];  % for mic inside booth
    thrshMic = [0.025 0.25];  % for mic inside booth
    thrshLA  = 2.5*[.1 .05];   % for left audio out
    show     = 0; % by default don't show the detection routine
    
    %processAna = 0;
    if processAna
        
        % identify if left audio (ch 2 or ch 5), photo-diode (ch 3) and mic (ch 4) are
        % present or not. If present determine best threshold
        % approach: more than half the time, there should be no signal on the
        % channels in question. Use the inter-quartile range to estimate
        % noise-floor of the channel.
        
        %interQuartile = quantile(all_anain_data(chanInd,:), [0.025 0.975 0.9999] );
        
        %tst = quantile(all_anain_data(chanInd,:), 0:0.01:1 );
        %plot(tst, 0:0.01:1)
        if chanInd == 2
            %thrshLA(1)  = 3 * diff( interQuartile(1:2) )
        end
        
        
        %keyboard
        stimOnInd = toneOnsetInd;
        vsi       = find(vsyncInd);
        
        for(i = 1:length(stimOnInd))
            
            thisInd   = (vsi(stimOnInd(i))-200):(vsi(stimOnInd(i))+100);
            thisTime  = all_t_amplifier_real( thisInd  ) - eventTime_real(         stimOnInd(i));   % plus minus 10 ms
            thisMic   = all_anain_data(chanInd, thisInd );   % plus minus 10 ms
            
            
            if chanInd ==4
                dM = [0 diff(thisMic)];
                ddM = [0 diff(dM)];
                
                tmp = abs(dM);
                noiseDiff = sort(tmp(1:100));
                
                tmp = tmp - noiseDiff(100);
                tmp(tmp<0) = 0;
                tst = cumsum( tmp );
                
                thisMicOn = min(find( tst > thrshMic(1) ));
            end
            
            % the left audio channel
            if chanInd == 2 || chanInd == 5
                %thisMicOn = min(find( abs(thisMic- mean(thisMic(1:100)) )>thrshLA(1) ));
                thisMicOn = min(find( abs(thisMic- min(thisMic(1:100)) )>thrshLA(1) ));
            end
            
            micCorrection(i)                   = -999;
            eventTimeCorrected( stimOnInd(i) ) = -999;
            
            show= 1;
            if(~isempty(thisMicOn))
                thisTime(thisMicOn);
                micCorrection(i) = thisTime(thisMicOn);
                show = 0;
            end
            
            show = 0;
            if(show)
                hold off
                plot(thisTime, abs(diff([thisMic(1),thisMic])))
                hold on
                plot(thisTime, thisMic - thisMic(1),'r')
                if chanInd == 4
                    plot(thisTime, tst);
                end
                
                plot( 0*[1 1], [-.3,.3] , 'k' )
                plot(  micCorrection(i) * [1 1], [-.5,.5] , 'r' )
                xlim([-0.02 0.02]);
                pause(.5)
            end
            show=0;
        end
        
        % ============ summarize the correction procedure
        eventTimeCorrected_real(stimOnInd) =    eventTime_real(stimOnInd) + micCorrection;
        eventTimeCorrected_fake(stimOnInd) =    eventTime_fake(stimOnInd) + micCorrection;
        
        plot( micCorrection )
        
        mean(micCorrection)
        std(micCorrection)
        % ===========================
    end
    
    %%
    if ~exist( [localDir '/' directory '/' exp{1}], 'dir' )
        mkdir( [localDir '/' directory '/' exp{1}] );
    end
    
    outfile = [localDir '/' directory '/' exp{1} '/events.txt'];
    fid = fopen(outfile, 'w');
    bl = ' ';
    %for eind = 1:length(eventTime)
    %    cur_line = [ num2str(eventCode(eind)) bl num2str(eventTimeCorrected(eind)) '\n'];
    %    fprintf(fid, cur_line);
    %end
    
    
    for eind = 1:(length(toneOnsetInd))
        cur_line = [ num2str(eventCode(toneOnsetInd(eind))) bl num2str(eventTimeCorrected_fake(toneOnsetInd(eind)) + soundDelay ) '\n'];
        fprintf(fid, cur_line);
    end
    
    fclose(fid);
    
    
end % startFromScratch


%% ====================================
%   use command line eeglab methods
%   to load and preproceess data set
%
clear raw_amp raw_amp5k raw_filt

% start with the raw_amp data set
% it is based on the 1kHz downsampled data, but does not include any
% filtering other than the one automatically done during downsampling
% the raw_amp data set then gives rise to 'raw' which is epoched and
% rerefferenced

% raw_amp is also low-pass filtered to give rise to raw_amp_filt.
% raw_amp_filt is then epoched and rereferenced to give rise to llp

if doRAW | doLLP
    % import data into EEG-lab format
    %raw_amp       = pop_importdata('setname','raw_amp','data',[localDir '/' directory '/' exp{1} '/amp.mat'], 'dataformat','matlab','nbchan',Namp+Nana,'xmin',0,'srate',real(1/tres_down),'ref','4x8_0x3');
    raw_amp       = pop_importdata('setname','raw_amp','data',down_amplifier_data, 'dataformat','matlab','nbchan',Namp+Nana,'xmin',0,'srate',real(1/tres_down),'ref','4x8_0x3');
    raw_amp.times = down_t_amplifier_fake;
    %raw_ana                 = pop_importdata('setname','raw_ana','data',[localDir '/' directory '/' exp{1} '/ana.mat'], 'dataformat','matlab','nbchan',8,'xmin',0,'srate',15000,'ref','right_mastoid');
    % add event information
    [raw_amp, eventnumbers] = pop_importevent(raw_amp, 'event',[localDir '/' directory '/' exp{1} '/events.txt'],  'fields',{'type','latency' }, 'append', 'no', 'timeunit', 1 );
    
    clear down_amplifier_data down_t_amplifier_fake down_t_amplifier_real
end


if doMLP | doMBP | doBSP | doRaw5k
    % same for the 5kS high-res data
    raw_amp5k                 = pop_importdata('setname','raw_amp','data',all_amplifier_data, 'dataformat','matlab','nbchan',Namp+Nana,'xmin',0,'srate',real(1/tres_amp),'ref','4x8_0x3');
    raw_amp5k.times           = all_t_amplifier_fake;
    [raw_amp5k, eventnumbers] = pop_importevent(raw_amp5k, 'event',[localDir '/' directory '/' exp{1} '/events.txt'],  'fields',{'type','latency' }, 'append', 'no', 'timeunit', 1 );
    
    clear all_amplifier_data all_t_amplifier_fake all_t_amplifier_real
end

% some old debugging code to see which events were missing
%clear tmp
%for jnd = 1:length(raw_amp.event)
%   tmp(jnd) = raw_amp.urevent(jnd).init_index;
%end

if doRAW
    %% save unfiltered data and un-rereferenced data
    raw = pop_epoch( raw_amp, valSoundStr, [-0.150 0.750]);
    
    raw.setname        = 'Raw';
    
    %raw                = pop_rmbase( raw,      [-50 0] );    % mod TT 20150722
    
    tmp                = raw.data( (Namp+1):(Namp+Nana),:,:); % don't rereference the analog input channels
    raw.data           = reref(raw.data,1, 'keepref', 'on');
    raw.data( (Namp+1):(Namp+Nana),:,:) = tmp;               % restore unreferenced analog input channels
    
    disc               = pop_saveset(raw    , 'filename','raw','filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
end

if doRaw5k

    epoch5k          = pop_epoch( raw_amp5k, valSoundStr ,[-0.100 0.200]);    
    epoch5k.setname  = 'raw5k';
    epoch5k.data     = reref(epoch5k.data,1, 'keepref', 'on');
    disc             = pop_saveset(epoch5k    , 'filename','raw5k','filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile', 'version', '7.3');

end


if doLLP
    %% llp data
    
    % band-pass filter the data
    %locutoff  = 2 ;  % 120 HZ highpass filter
    %hicutoff  = 25;  %
    %filtorder = 826;%1652; %4096;
    
    % this bandpass seems to introduce ringing artefacts. Hence I currently
    % don't use a high-pass at all.
    % actually, the ringing is introduced by the low-pass...
    
    locutoff  = .2 ;  % .1 HZ highpass filter
    %hicutoff  = 50;  %
    %filtorder = 1024; %
    hicutoff  = 70;   % 80
    filtorder = 128;  % 256
       
    ['transition bandwidth = ' num2str(500* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/500*     filtorder) 'ms']
    m = pop_firwsord('blackman', 500, 20)
    %
    
    [raw_amp_filt, com, b] = pop_firws(raw_amp, 'fcutoff',hicutoff, 'forder',filtorder,'ftype','lowpass');
    
    %[raw_amp_filt, com, b] = pop_eegfiltnew(raw_amp, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    % allow filtering of LeftAudio
    raw_amp_filt.data( (Namp+1):(Namp+Nana),:) = raw_amp.data( (Namp+1):(Namp+Nana),:);
    
    llp = pop_epoch( raw_amp_filt, valSoundStr, [-0.150 0.750]);
    clear raw_amp_filt
    
    llp.setname        = 'LongLatencyPotentials';
    %llp                = pop_rmbase( llp,      [-50 0] );  % mod TT 20150722
    
    tmp                = llp.data( (Namp+1):(Namp+Nana),:,:); % don't rereference the analog input channels
    llp.data           = reref(llp.data,1, 'keepref', 'on');
    llp.data( (Namp+1):(Namp+Nana),:,:) = tmp;               % restore unreferenced analog input channels
    
    
    thisFileName = 'llp_2Hz';
    disc     = pop_saveset(llp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');


    
    
    clear llp
    % run short filter settings
    
    filtorder = 28;
    hicutoff  = 20;
    
    ['transition bandwidth = ' num2str(500* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/500*     filtorder) 'ms']
    
    % 
    m = pop_firwsord('blackman', 500, 20)
    %
    
    [raw_amp_filt, com, b] = pop_firws(raw_amp, 'fcutoff',hicutoff, 'forder',filtorder,'ftype','lowpass');
    
    %[raw_amp_filt, com, b] = pop_eegfiltnew(raw_amp, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    % allow filtering of LeftAudio
    raw_amp_filt.data( (Namp+1):(Namp+Nana),:) = raw_amp.data( (Namp+1):(Namp+Nana),:);
    
    llp = pop_epoch( raw_amp_filt, valSoundStr, [-0.150 0.750]);
    clear raw_amp_filt
    
    llp.setname        = 'LongLatencyPotentials';
    %llp                = pop_rmbase( llp,      [-50 0] );  % mod TT 20150722
    
    tmp                = llp.data( (Namp+1):(Namp+Nana),:,:); % don't rereference the analog input channels
    llp.data           = reref(llp.data,1, 'keepref', 'on');
    llp.data( (Namp+1):(Namp+Nana),:,:) = tmp;               % restore unreferenced analog input channels
    
    
    thisFileName = 'llp2_short';
    disc     = pop_saveset(llp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
 
    
    %% make sure tone and all have same number of events
    if length(tone.trial) > size(llp.data,3)
        'warning: removing events from tone. Double check if something is off'
        sall         = size(llp.data,3);
        tone.trialNr = tone.trialNr(1:sall);
        tone.time    = tone.time(1:sall);
        tone.seqNr   = tone.seqNr(1:sall);
        tone.pitch   = tone.pitch(1:sall);
        tone.trial   = tone.trial(1:sall);
        tone.type    = tone.type(1:sall);
        tone.ISI     = tone.ISI(1:sall);
        tone.ampInd  = tone.ampInd(1:sall);
        tone.deltaPitch   = tone.deltaPitch(1:sall);
        tone.relPosition  = tone.relPosition(1:sall);
    end
    
    save( [localDir '/' directory '/' exp{1} '/tone.mat'], 'tone', 'trial');
    clear raw_amp llp raw
end

if doMLP
    %% mid-latency components
    % will be based on the 5kHz data set 'all_amplifier_data'
    
    % setting for the main analysis
    hicutoff  = 500;  %
    filtorder = 512;%% use 280 for 5kS % use 56 for 1kS
    locutoff  = 30; % 50
    
    ['transition bandwidth = ' num2str(5000* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/5000*     filtorder) 'ms']
    
    [raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', [locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');
    %[raw_amp_filt, com, b] = pop_eegfiltnew(raw_amp5k, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    %raw_amp_filt = raw_amp5k;
    mlp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.010 0.110]);
    
    mlp.setname        = 'mlp';
    
    %mlp      = pop_rmbase( mlp,      [-20 0] );
    mlp.data = reref(mlp.data,1, 'keepref', 'on');
    
    thisFileName = 'mlp';
    if shortFilter
        thisFileName = 'mlp_short';
    end
    disc     = pop_saveset(mlp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
    
    
    
    clear mlp
    % run short filter settings
    locutoff = 20;
    filtorder = 128;
    hicutoff  = 500;
    
    ['transition bandwidth = ' num2str(5000* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/5000*     filtorder) 'ms']
    
    [raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', [locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');
    %[raw_amp_filt, com, b] = pop_eegfiltnew(raw_amp5k, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    %raw_amp_filt = raw_amp5k;
    mlp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.010 0.110]);
    
    mlp.setname        = 'mlp';
    
    %mlp      = pop_rmbase( mlp,      [-20 0] );
    mlp.data = reref(mlp.data,1, 'keepref', 'on');
    
    thisFileName = 'mlp_short';
    disc     = pop_saveset(mlp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
    
    %disc     = pop_saveset(mlp,  'filename','mlp','filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
end

if doMBP
    %% mid-brain potentials?
    % will be based on the 5kHz data set 'all_amplifier_data'
    
    % setting for main analysis 20150903; bandpass
    locutoff  = 60;
    hicutoff  = 1000;
    filtorder = 512;
    
    % bandpass not usable, new filter settings: 20150915, only mild
    % low-pass filtering 
    hicutoff  = 500;  % 125 -> 0.5  
    locutoff  = -1;
    filtorder = 30;
    
    ['transition bandwidth = ' num2str(5000* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/5000*     filtorder) 'ms']

    
    %[raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', [locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');
    [raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', hicutoff, 'forder',filtorder,'ftype','lowpass');
    
    %[raw_amp_filt, com, b] = pop_eegfiltnew(raw_amp5k, locutoff, hicutoff, filtorder,    0, 0, 0, 0);
    mbp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.005 0.020]);
    
    bsp.setname        = 'mbp';
    
    %mlp      = pop_rmbase( mlp,      [-20 0] );
    mbp.data = reref(mbp.data,1, 'keepref', 'on');
    
    thisFileName = 'mbp';
    if shortFilter
        thisFileName = 'mbp_short';
    end
    disc     = pop_saveset(mbp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
    
    
    % short filter setting
    clear mbp
    
    hicutoff  = 500;  % 125 -> 0.5
    locutoff  = -1;
    filtorder = 30;
    
    ['transition bandwidth = ' num2str(5000* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/5000*     filtorder) 'ms']
    [raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', hicutoff, 'forder',filtorder,'ftype','lowpass');
   
    mbp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.005 0.020]);
    
    bsp.setname        = 'mbp';
    
    %mlp      = pop_rmbase( mlp,      [-20 0] );
    mbp.data = reref(mbp.data,1, 'keepref', 'on');
    thisFileName = 'mbp_short';
    
    disc     = pop_saveset(mbp,  'filename',thisFileName, 'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
    
end


if doBSP
    %% brain-stem potentials?
    % will be based on the 5kHz data set 'all_amplifier_data'
    locutoff  = 80;  % 75 0 86
    hicutoff  = 1000;
    filtorder = 512;
    
    [raw_amp_filt, com, b] = pop_firws(raw_amp5k, 'fcutoff', [locutoff hicutoff], 'forder',filtorder,'ftype','bandpass');
    
    ['transition bandwidth = ' num2str(5000* 5.5/filtorder) 'Hz']
    ['filter length        = ' num2str(1000/5000*     filtorder) 'ms']
        
    bsp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.005 0.012]);
    
    bsp.setname        = 'bsp';
   
    bsp.data = reref(bsp.data,1, 'keepref', 'on');
    thisFileName = 'bsp';
    disc     = pop_saveset(bsp,  'filename',thisFileName,'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile');
    
    
    
   % short filter setting (in this case, no filter)
   clear bsp 
   raw_amp_filt = raw_amp5k;      
    bsp = pop_epoch( raw_amp_filt, valSoundStr ,[-0.005 0.012]);
    
    bsp.setname        = 'bsp_short';
   
    bsp.data = reref(bsp.data,1, 'keepref', 'on');
    thisFileName = 'bsp_short';
    disc     = pop_saveset(bsp,  'filename',thisFileName,'filepath', [localDir '/' directory '/' exp{1} '/'], 'savemode','onefile'); 
    
    
end

%%
disc = ampOddClick2R(exp{1}, 0); % export data
'done!'

end


