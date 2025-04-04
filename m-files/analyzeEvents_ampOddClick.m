function [ tone, trial, valSoundStr ] = analyzeEvents_ampOddClick( eventCodeCnt, ana, all_t_amplifier_real, all_t_amplifier_fake, args, impedance )
% analyzeEvents_RSkernel
%   part of the analyzeEvents series of functions
%   takes continuous eventCode as input and reports back two structures
%   tone and trial that describe the paradigm
%keyboard
SR = mean(diff(all_t_amplifier_real));
valCodes = [1:48 65:75 99 116]; % only accept these triggers, everything else is caused by loose connection

% Find all events and event
vsyncIndTmp       = find(abs(diff([eventCodeCnt(1) eventCodeCnt]))>0.1);
eventCodeTmp      = eventCodeCnt(vsyncIndTmp);

% remove events that are not on the list
vsyncInd         = vsyncIndTmp( find(ismember(eventCodeTmp,valCodes)) );
eventCode        = eventCodeCnt(vsyncInd);
eventCode(find(eventCode==48)) = 32;

% remove events that are repeated after removing non-members
eventCodeChange = diff([-100 eventCode]);
vsyncInd        = vsyncInd( eventCodeChange~=0 );
eventCode       = eventCode( eventCodeChange~=0 );

eventTime_real   = all_t_amplifier_real(vsyncInd);
eventTime_fake   = all_t_amplifier_fake(vsyncInd);


% find which tones are being used in this paradigm
for i = 1:23
    tbl(i) = length(find(eventCode==i));
end
valSounds = find(tbl>1);
['Valid Sounds: ', int2str(valSounds) ]
soundStr    = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
valSoundStr = soundStr(valSounds);

channel16Broken = 0;
if channel16Broken
    
    toneEndInd    = [1 find(eventCode==32)]; 
    delta32       = diff(toneEndInd);
    toneBefore    = zeros(1,length(toneEndInd)-1);       % only one candidate tone
    toneBeforeMod = zeros(1,length(toneEndInd)-1);       % only one candidat tone modulo 16
    timeBefore    = zeros(1,length(toneEndInd)-1 );      % time of that one tone
    timeBeforeMod = zeros(1,length(toneEndInd)-1 );      % time of the first candiate modulo16
    indexBefore    = zeros(1,length(toneEndInd)-1 );     % index of the one candidate
    indexBeforeMod = zeros(1,length(toneEndInd)-1 );     % index of the first candidate modulo16
    
    
    for t = 1:(length(toneEndInd)-1)
       tmp            = eventCode(toneEndInd(t):toneEndInd(t+1));
       cds            = [0 tmp(1:(end-1))]==99;% codes
       cand           = ismember(tmp,valSounds); %tone candidates
       potentialIndeces = find(~cds & cand);     % index of all viable candidates
       potentialTones = tmp( find(~cds & cand) );
       
       if length(potentialTones)==1
            indexBefore(t)    = toneEndInd(t) + find(~cds & cand) - 1;
            indexBeforeMod(t) = toneEndInd(t) + find(~cds & cand) - 1;
            toneBefore(t)     = potentialTones;
            toneBeforeMod(t)  = mod(potentialTones-1,16)+1;
            timeBefore(t)     = eventTime_real( indexBefore(t) );
            timeBeforeMod(t)  = eventTime_real( indexBefore(t) );
       end
           
       if length(potentialTones)>1
          indexBeforeMod(t)  = toneEndInd(t) + min(find(~cds & cand)) - 1;
          timeBeforeMod(t)   = eventTime_real( indexBeforeMod(t) );
          
          utmp = unique(mod(potentialTones-1,16)+1);
          if(length(utmp)==1)
             toneBeforeMod(t) = potentialTones(end);%utmp; 
             eventCode( toneEndInd(t)-1+potentialIndeces) = potentialTones(end);
          end
       end
    end
    
    eventCodeChange = diff([-100 eventCode]);
    vsyncInd        = vsyncInd( eventCodeChange~=0 );
    eventCode       = eventCode( eventCodeChange~=0 );

    eventTime_real   = all_t_amplifier_real(vsyncInd);
    eventTime_fake   = all_t_amplifier_fake(vsyncInd);
end

% soa = diff([-10, timeBeforeMod]); 
% %novel = -1*ones(1,length(toneEndInd)-1); 
% novelTone = abs(diff( [-1 toneBeforeMod]))>0 ;
% novelTime = abs(diff([50 soa]))>.020;
% 
% novel = abs(diff([50 soa]))>.020 | abs(diff( [-1 toneBeforeMod]))>0 ;
% for t = 2:(length(toneEndInd)-1)
%     
%     
% end


Ind45               = find(eventCode == 116);
%deltaTime           = eventTime_real(Ind45(2:end),);
trialOnsetInd       = Ind45( find(  (diff([eventTime_real(Ind45) max(eventTime_real)]) > 0.1) &  (diff([Ind45 length(eventCode)]) > 10)) );
%trialOnsetInd       = Ind45( find(  (diff([0 eventTime_real(Ind45)]) > 0.1) &  (diff([-250 Ind45]) > 200)) );
if(trialOnsetInd(1)>5) % first trial Onset trigger was missed
    trialOnsetInd = [1; trialOnsetInd];
end

trialOnsetTime_real = all_t_amplifier_real(vsyncInd(trialOnsetInd));
trialOnsetTime_fake = all_t_amplifier_fake(vsyncInd(trialOnsetInd));

% calculate trial number for each event index
trialInd                = zeros(1,length(eventCode));
trialInd(trialOnsetInd) = 1;
trialNrAll              = cumsum(trialInd);


% find all tone onsets
toneOnsetInd = find(ismember(eventCode, valSounds));

% remove the trial info triggers from the list of potential toneOnsets
tmpOnsetInd      = toneOnsetInd(find(toneOnsetInd>1)); % toneOnsetInd should never be the first index, but it can happen if ML does not get turned off after intan crash
trialInfoTrigger = find(eventCode(tmpOnsetInd-1)==99);
toneOnsetInd     = toneOnsetInd(  setdiff(1:length(toneOnsetInd), trialInfoTrigger )  );

% there are some random events that happen between trials, but are
% not associated with an actual tone-presentation. These can be
% recognized by a lack of a 32 marker immediately following it
keep                             = eventCode(toneOnsetInd+1)==32;
eventCode(toneOnsetInd(keep==0)) = 999;
toneOnsetInd                     = toneOnsetInd(keep);
['Found ', int2str(length(find(keep==0))),' non-terminated sounds']


% trial number for tone onsets only
tone.trialNr                 = trialNrAll(toneOnsetInd);
if(tone.trialNr(1)==0)
    % in case the trial onset trigger was not recorded on the first trial
    tone.trialNr = tone.trialNr + 1;
end

tone.time              = eventTime_real(toneOnsetInd);
tone.ISI               = [-1 diff(tone.time)];

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
tone.pitch      = eventCode(toneOnsetInd);
tone.deltaPitch = [-99 diff(tone.pitch)];


%  not being calcuated here. just a placeholder
% calcluate trial type of each tone: 1: oddball; 2: controll mean; 3: controll abs 1; 4: controll abs 2;
tone.trial   = ones(size(tone.pitch)); %trialTypeMap(tmp);

% calcluate the tone type: 0: repeat; 1:novel
tone.type   = ones(size(tone.pitch));

% ---------------------- trial
devInd                                      = find(tone.type==1);
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
tmp             = eventCode(trialOnsetIndVal+3);  % these values should all be 99
if( length(find(tmp~=99))>0 )
    %keyboard
    tmp
    trialOnsetInd
    'no trialType information available in the digital codes'
end

if length(find(tmp==99))>0
    trial.type(find(tmp==99))      = eventCode(trialOnsetIndVal(find(tmp==99))+4) - 64;
end

% calculate the sequence number of the deviant for each trial; use 14 as
% default
tone.deviantPosition = trial.deviantPosition(tone.trialNr);
tone.relPosition     = tone.seqNr - tone.deviantPosition;
tone.deviantPitch    = trial.deviantPitch(tone.trialNr);


%% =================== correct timing based on audio-trigger data

eventTimeCorrected_real = eventTime_real;
eventTimeCorrected_fake = eventTime_fake;

% ana has only one channel, select the correct one in call to read_intan_function

%thrsh = [0.02 0.015];  % for mic inside booth
%thrshMic = [0.025 0.25];  % for mic inside booth
thrshLA  = 2.5*[.1 .05];   % for left audio out
show     = 0; % by default don't show the detection routine

if length(ana)==length(eventCodeCnt)
    
    % identify if left audio (ch 2 or ch 5), photo-diode (ch 3) and mic (ch 4) are
    % present or not. If present determine best threshold
    % approach: more than half the time, there should be no signal on the
    % channels in question. Use the inter-quartile range to estimate
    % noise-floor of the channel.
    
    stimOnInd = toneOnsetInd;
    vsi       = vsyncInd; %find(vsyncInd);
    
    for(i = 1:length(stimOnInd))
        
        thisInd   = (vsi(stimOnInd(i))- (.040/SR)):(vsi(stimOnInd(i))+ (.040/SR));
        thisTime  = all_t_amplifier_real( thisInd  ) - eventTime_real(         stimOnInd(i));   % plus minus 10 ms
        thisMic   = ana(1, thisInd );                                                            % plus minus 10 ms
        
        
        if 1 == 0
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
        if 1 == 1
            thisMicOn = min(find( abs(thisMic- mean(thisMic(1:floor(0.010/SR))) )>thrshLA(1) ));
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
        if micCorrection(i) == -999
            show = 1;
        end
        
        if(show)
            i
            keyboard
            hold off
            plot(thisTime, abs(diff([thisMic(1),thisMic])))
            hold on
            plot(thisTime, thisMic - thisMic(1),'r')
            
            %if chanInd == 4
            %    plot(thisTime, tst);
            %end
            
            plot( 0*[1 1], [-.3,.3] , 'k' )
            plot(  micCorrection(i) * [1 1], [-.5,.5] , 'r' )
            xlim([-0.02 0.02]);
            pause(.5)
            %keyboard
        end
        show=0;
    end
    
    % ============ summarize the correction procedure
    failInd = find(micCorrection<-100);
    if length(failInd)>0
        warning('failed to find correction for some tones')
        length(failInd)
        micCorrection(failInd) = 0;
    end
    
    eventTimeCorrected_real(stimOnInd) =    eventTime_real(stimOnInd) + micCorrection;
    eventTimeCorrected_fake(stimOnInd) =    eventTime_fake(stimOnInd) + micCorrection;
    
    tone.time  = eventTime_fake(stimOnInd) + micCorrection;
    tone.ISI   = [-1 diff(tone.time)];
    
    plot( micCorrection )
    
    mean(micCorrection)
    std(micCorrection)
    % ===========================
end

%% ======= write out 
if ~exist( [args.localDir '/' args.directory '/' args.exp{1}], 'dir' )
        mkdir( [args.localDir '/' args.directory '/' args.exp{1}] );
end

outfile = [args.localDir '/' args.directory '/' args.exp{1} '/events.txt'];
fid     = fopen(outfile, 'w');
bl      = ' ';

for eind = 1:(length(toneOnsetInd))
    cur_line = [ num2str(eventCode(toneOnsetInd(eind))) bl num2str(tone.time(eind) + args.soundDelay ) '\n'];
    fprintf(fid, cur_line);
end
fclose(fid);

%save( [args.localDir '/' args.directory '/' args.exp{1} '/tone.mat'], 'tone', 'trial');
SLsaveObj( [args.localDir '/' args.directory '/' args.exp{1} '/'], 'tone', tone, trial  )

end

