function [ disc ] = RSflex2R( exp, Econfig, SR)

% ampOddClick2R exports summary data from a aMMN experiment to R

% writing out of .txt file does not work on tinyone. Some annoying issue
% with the samba drive

'this is ampOddClick2R'

doRaw   = 1;
doMUA   = 1;
doFFR   = 0;

% preliminaries
%baseDir         = '';
baseDir         = '/Volumes/Drobo5D3/EEG/';

%if exist('/Volumes/rawData/EEG/') 
%    baseDir         = '/Volumes/rawData/EEG/';
%end
%baseDir         = '~/Dropbox/Nibbler/';
IntanComp       = 'PlanetExpress';
monkeyLogicComp = 'Leela';
directory       = 'ampOddClick'; 
localDir        = [baseDir 'EEGLab/']; %'/Volumes/rawData/EEG/EEGLab/';

%localDir        = '/Volumes/rawData/EEG/EEGLab/';

disc            = 1;
RDir            = [baseDir directory '/rda/' exp '/'];
RDirtxt         = ['~/Dropbox/ampOdd_click/Rdata/' exp '/'];

if ~exist(RDir)
    mkdir(RDir)
end

if ~exist(RDirtxt)
    mkdir(RDirtxt)
end


outfile         = [RDir exp '.txt'];
outfile2        = [RDirtxt exp '.txt'];

outfileraw      = [RDir exp '_raw.txt'];
outfile2raw     = [RDirtxt exp '_raw.txt'];


%cmpList         = {'p1','nx','px','n1','p2','n2','frontoCentral','p3','p3b'};
cmpList         = {'p1x','p1a','p1b','p1','nx','px','n1','p2','n2'};%,'frontoCentral','p3','p3b'};
cmpListmlp      = {'p1x','p1a','p1Ta','p1b','p1Tb','p1c'};
cmpListmbp      = {'pSV','n8','n8p'};
fileextn        = '';


saveAll         = 1;   % overwrite existing data?  
saveChannels    = 1;   % save if not existant
saveComponents  = 1;   % save if not existant

% load data
load( [localDir '/' directory '/' exp '/tone.mat']);
trial.type
tone.ampInd = mod(tone.pitch,11);

% load the data
raw     = pop_loadset(['raw_' num2str(SR) '.set'],[localDir '/' directory '/' exp '/']);
llp     = pop_loadset(['LLP_' num2str(SR) '.set'],[localDir '/' directory '/' exp '/']);
mua     = pop_loadset(['MUA_' num2str(SR) '.set'],[localDir '/' directory '/' exp '/']);
%ffr     = pop_loadset(['FFR_' num2str(SR) '.set'],[localDir '/' directory '/' exp '/']);

%raw = lineNoiseNotchEEG( raw, [60 120], 34:(size(llp.data,1)-2),0 );
%llp = lineNoiseNotchEEG( llp, [60 120], 34:(size(llp.data,1)-2),0 );


%% 
tmp          = strsplit(exp,'_');
animal       = tmp(1);
valChan      = [20 21];

llpsr  = SR;          %% sampling rate of llp2
mxLLP  = 500.0000001; % maximum time for llp

all        = llp;
valChan    = [1 2 3 5 6 7    8 9 10 12 13    14 15 16 17 18 19];
if strcmp(animal,'Jesse')
    valChan      = [1 2 3 5 6    8 9 10 12    14 15 16 17 18];
end
if strcmp(animal,'Sam')
    valChan      = [1 2 3 5 6 7    9 10 12    14 15 16 17 18];
    
    if strcmp(exp,'Sam_20161020_1230_zm6350') | strcmp(exp,'Sam_20161020_1420_zm6350')
        valChan = [1 2 3     7    9 10 12    14 15 16 17 18];
    end
    
end

absEEG     = squeeze( max(max(abs(all.data(valChan,:,:)),[],2), [],1) );
minEEG     = squeeze( min(min(all.data(valChan,:,:),[],2), [],1) );
maxEEG     = squeeze( max(max(all.data(valChan,:,:),[],2), [],1) );
p2p        = maxEEG - minEEG;


% get peak-to-peak for vprobe channels only
all.data = reref(all.data,34,'keepref','on');
valVPchan = 34:( size(llp.data,1)-2 );

absVP     = squeeze( max(max(abs(all.data(valVPchan,:,:)),[],2), [],1) );
minVP     = squeeze( min(min(all.data(valVPchan,:,:),[],2), [],1) );
maxVP     = squeeze( max(max(all.data(valVPchan,:,:),[],2), [],1) );
p2pVP     = maxVP - minVP;


clear all

if strcmp(animal,'Jesse')
    valChan      = 20:33;
end

if strcmp(animal,'Walter')
    valChan      = 20:31;
end

if strcmp(animal,'Sam')
    valChan      = 20:31;
end

powEOG     = squeeze( sum( var(mua.data(valChan,:,:),[],2),1) );


%% ==================================================
% new correction routine:
% six different SOA-ranges: 0.2-0.4; 0.4-0.8; 0.8-1.6; 1.6-3.20;3.20-6.40; 6.40-25.00 

tmpllp      = pop_rmbase( llp,      [-50 0] );    % mod TT 20150722
tmpraw      = pop_rmbase( raw,      [-50 0] );    % mod TT 20150722

taper        = ones(1,length(tmpllp.times) );
tpind        = find(tmpllp.times > max(tmpllp.times)-50);
taper(tpind) = 1:(-1/(length(tpind)-1)):0;

tone.futureISI = [tone.ISI(2:length(tone.ISI)) 25];
%brks = [0.0190 0.400 0.800 1.600 3.200 6.400 25];
brks = [0.250 0.500 1.000 2.000 4.000 8.000 25];

regularOnly = sum( tone.pitch>11 )/length(tone.pitch) < .25;
tooFew      = sum( tone.pitch>11 ) < 750;

if ~regularOnly || tooFew
    
    for (snd = 1:6)
        for (and = 1:5)
            vlnd = find( p2pVP' < quantile(p2pVP,.9) & tone.ISI > brks(snd) & tone.ISI < brks(snd+1) & tone.pitch==and+11 & tone.futureISI>.700); % & tone.trial==1 & tone.type==1);
            length(vlnd)
            crrct{ snd, and}  = mean( tmpllp.data( :, :, vlnd   ), 3);
            crrctR{snd, and}  = mean( tmpraw.data( :, :, vlnd   ), 3);
            
            for i = 1:size(crrct{snd,and},1)
                crrct{ snd, and}(i,:)  = crrct{ snd, and}(i,: ).*taper;
                crrctR{snd, and}(i,:)  = crrctR{snd, and}(i,: ).*taper;
            end
        end
    end
    
    brks(7)   = 10000; % make sure that snd is never larger than 6.
    tmpSOA    = tone.ISI;
    tmpSOA(1) = 25;
    tmpSOA( find(tmpSOA>25) ) = 25;
    
    for i = 2:length(tone.ISI)
        if tone.ISI(i)<.700
            %tone.ISI(i)
            fromInd = find( round(tmpllp.times) == 1000*1/llpsr*round(llpsr*(tone.ISI(i)-.150)) ):length(tmpllp.times);
            toInd   = 1:length(fromInd);
            
            snd = sum( tmpSOA(i-1) > brks );
            and = mod(tone.pitch(i-1),11);
            
            llp.data(:,toInd,i)  = llp.data(:,toInd,i) -  crrct{snd,and}(:,fromInd);
            %raw.data(:,toInd,i)  = raw.data(:,toInd,i) - crrctR{snd,and}(:,fromInd);
        end
    end
    
end


%taper        = ones(1,length(tmpllp.times) );
%tpind        = find(tmpllp.times > max(tmpllp.times)-50);
%taper(tpind) = 1:(-1/(length(tpind)-1)):0;

%for (snd = 1:6)
%    for (and = 1:5)
%        for i = 1:size(crrct{snd,and},1)
%            crrct{ snd, and}(i,:)  = crrct{ snd, and}(i,: ).*taper;
%            crrctR{snd, and}(i,:)  = crrctR{snd, and}(i,: ).*taper;
%        end
%    end
%end

%% entire channel data (P3 corrected but not baseline corrected)
if saveChannels == 1
    
    taxis = llp.times;%(tind);
    %save( [RDir    exp '_taxis.mat'],'taxis' )
    %save( [RDirtxt exp '_taxis.mat'],'taxis' )
    SLsave6( RDir,    [exp '_taxis'], taxis );
    
     for ch = 1:size(llp.data,1)
        % long-latency potentials
        
        fileextn = '';
        if ~exist( [RDir exp '_ch', int2str(ch) fileextn '.mat'],'file'  ) || saveAll == 1
            chdat = squeeze( llp.data(ch,:,:));
            %save( [RDir exp '_ch', int2str(ch) fileextn '.mat'],'chdat' )
            SLsave6( RDir, [exp '_ch', int2str(ch) fileextn ], chdat )
        end
        
    end

  if doRaw
        %tind  = find( raw.times>= -150.0000001 & raw.times<=mxLLP );
        taxis = raw.times; %(tind);
        %save( [RDir    exp '_taxis_raw.mat'],'taxis' )
        %save( [RDirtxt exp '_taxis_raw.mat'],'taxis' )
        SLsave6( RDir, [exp '_taxis_raw'], taxis )
        
        for ch = 1:size(raw.data,1)
            % raw data
            if ~exist( [RDir exp '_ch', int2str(ch) '_raw.mat'],'file'  ) || saveAll == 1
                chdat = squeeze( raw.data(ch,:,:));
                %save( [RDir exp '_ch', int2str(ch) '_raw.mat'],'chdat' )
                SLsave6( RDir, [exp '_ch', int2str(ch) '_raw'],chdat )
            end
        end
    end
    
    
    if doMUA
        taxis = mua.times; %(tind);
        SLsave6( RDir, [exp '_taxis_mua'], taxis )
        
        for ch = 1:size(mua.data,1)
            % raw data 5k
            if ~exist( [RDir exp '_ch', int2str(ch) '_mua.mat'],'file'  ) || saveAll == 1
                chdat = squeeze( mua.data(ch,:,:));
                SLsave6( RDir, [exp '_ch', int2str(ch) '_mua'],chdat )
            end
        end
    end

    if doFFR
        taxis = ffr.times; %(tind);
        SLsave6( RDir, [exp '_taxis_ffr'], taxis )
        
        for ch = 1:size(mua.data,1)
            % raw data 5k
            if ~exist( [RDir exp '_ch', int2str(ch) '_ffr.mat'],'file'  ) || saveAll == 1
                chdat = squeeze( ffr.data(ch,:,:));
                SLsave6( RDir, [exp '_ch', int2str(ch) '_ffr'],chdat )
            end
        end
    end

    
    
end


%% group channels by component (P3 corrected but not baseline corrected)
cmpListExpanded                    = cmpList;
cmpListExpanded{length(cmpList)+1} = 'frontoCentral';

if saveComponents==1
    for cmpInd = 1:(length(cmpListExpanded))
        
        
        [ trng ch signalStr ] = getComponentParameters(cmpListExpanded{cmpInd}, animal);
        
        %tind  = find( llp.times>= -500.000001 & llp.times<=mxLLP );
        taxis = llp.times;%(tind);
        
        % LLP2
        fileextn = '';
        if ~exist( [RDir exp '_' cmpListExpanded{cmpInd} fileextn '.mat'],'file' ) || saveAll == 1
            chdat = squeeze( mean( llp.data(ch,:,:)));
            %save( [RDir exp '_' cmpListExpanded{cmpInd} fileextn '.mat'],'chdat' )
            %save( [RDir exp '_taxis.mat'],'taxis' )
            
            SLsave6( RDir, [exp '_' cmpListExpanded{cmpInd} fileextn],chdat )
            SLsave6( RDir, [exp '_taxis'],taxis )
        end
        if strcmp(cmpListExpanded{cmpInd},'frontoCentral')
            %chdat = squeeze( mean( llp.data(ch,:,:)));
            %save( [RDirtxt exp '_' cmpListExpanded{cmpInd} fileextn '.mat'],'chdat' )
            %save( [RDirtxt exp '_taxis.mat'],'taxis' )
        end
        
        
        
        if doRaw
            %tind  = find( raw.times>= -500.000001 & raw.times<=mxLLP );
            taxis = raw.times;%(tind);
            
            % RAW
            if ~exist( [RDir exp '_' cmpListExpanded{cmpInd} '_raw.mat'],'file'  ) || saveAll == 1
                chdat = squeeze( mean( raw.data(ch,:,:)));
                %save( [RDir exp '_' cmpListExpanded{cmpInd} '_raw.mat'],'chdat' )
                %save( [RDir exp '_taxis_raw.mat'],'taxis' )
                
                SLsave6( RDir, [exp '_' cmpListExpanded{cmpInd} '_raw'],chdat )
                SLsave6( RDir, [exp '_taxis_raw'],taxis )
            end
            if strcmp(cmpListExpanded{cmpInd},'frontoCentral')
                %chdat = squeeze( mean( raw.data(ch,:,:)));
                %save( [RDirtxt exp '_' cmpListExpanded{cmpInd} '_raw.mat'],'chdat' )
                %save( [RDirtxt exp '_taxis_raw.mat'],'taxis' )
            end
        end
          
    end
end 

%% =============== baseline correction for all three data types
llp                = pop_rmbase( llp,     [-50 0] );    % mod TT 20150722

if doRaw
    raw                = pop_rmbase( raw,      [-50 0] );    % mod TT 20150722
end

if doMUA
    mua                = pop_rmbase( mua,      [-50 0] );    % mod TT 20151001
end


%% key values for llp (P3 and baseline corrected)
mn         = zeros( length(tone.trialNr), length(cmpList) );
%mns        = zeros( length(tone.trialNr), length(cmpList) );

for i = 1:length(cmpList)
    
    [ trng, ch, signalStr ]  = getComponentParameters(cmpList{i}, animal);
    tind                     = find( llp.times>=trng(1) & llp.times<trng(2) );
    mn(:,i)                  = squeeze( mean(mean(llp.data(ch,tind,:),1 ),2) );
    %mns(:,i)                 = squeeze( mean(mean(llps.data(ch,tind,:),1 ),2) );
    
    
end

%% key values for raw (Only baseline corrected)
if doRaw
    mnraw        = zeros( length(tone.trialNr), length(cmpList) );
    
    for i = 1:length(cmpList)
        
        [ trng, ch, signalStr ] = getComponentParameters(cmpList{i}, animal);
        
        tind                    = find( raw.times>=trng(1) & raw.times<trng(2) );
        mnraw(:,i)              = squeeze( mean(mean(raw.data(ch,tind,:),1 ),2) );
            
    end
end


%% some quick plots
figure
lwd = 2;
ch = [5 6 12 13 18 19];
valInd = find(p2p<450 & powEOG<2000);

plot( llp.times,  squeeze( mean(mean( llp.data( ch, :, valInd   ), 3),1))  ,'k' , 'linewidth',lwd);
hold on
%plot( llp.times,  squeeze( mean(mean( tmpllp.data( ch, :, valInd   ), 3),1))  ,'b' , 'linewidth',lwd);
%plot( llp2.times,  squeeze( mean(mean( raw.data( ch, :, valInd   ), 3),1))  ,'r' , 'linewidth',lwd);
plot( llp.times,  zeros( 1, length(llp.times) ) )


figure
for i = 1:size(llp.data,1)
    plot( llp.times,  10*i + squeeze( mean(mean( llp.data( i, :, valInd  ), 3),1))  ,'k' , 'linewidth',lwd);
    hold on
end
hold off



% all components by ISI
valInd1 = find(p2p<450 & powEOG<2000 & tone.ampInd'==1 & tone.ISI'>1);
valInd2 = find(p2p<450 & powEOG<2000 & tone.ampInd'==2 & tone.ISI'>1);
valInd3 = find(p2p<450 & powEOG<2000 & tone.ampInd'==3 & tone.ISI'>1);
valInd4 = find(p2p<450 & powEOG<2000 & tone.ampInd'==4 & tone.ISI'>1);
valInd5 = find(p2p<450 & powEOG<2000 & tone.ampInd'==5 & tone.ISI'>1);


figure
lwd = 2;
plot( llp.times,  squeeze( mean( raw.data( 5, :, valInd1   ), 3))  ,'k' , 'linewidth',lwd);
hold on
plot( llp.times,  squeeze( mean( raw.data( 5, :, valInd2   ), 3))  ,'y' , 'linewidth',lwd);
plot( llp.times,  squeeze( mean( raw.data( 5, :, valInd3   ), 3))  ,'g' , 'linewidth',lwd);
plot( llp.times,  squeeze( mean( raw.data( 5, :, valInd4   ), 3))  ,'b' , 'linewidth',lwd);
plot( llp.times,  squeeze( mean( raw.data( 5, :, valInd5   ), 3))  ,'r' , 'linewidth',lwd);

%% export data  (P3 and baseline-corrected, based on llp2 data)

fid  = fopen(outfile, 'w');
fid2 = fopen(outfile2, 'w');

%cmpListmlp  = strcat({'m.'},cmpListmlp);  % add mlp. prefix to avoid getting the same names for llp and mlp
%cmpListmlps = strcat(cmpListmlp,{'.s'});
%cmpListmlpr = strcat(cmpListmlp,{'.r'});
%cmpLists    = strcat(cmpList,{'.s'});
cmpListr    = strcat(cmpList,{'.r'});


R_VAR_TONE      = 'trialNr seqNr ampInd trialType toneType ISI deltaPitch time ';
%R_VAR_MBP       =  sprintf('%s ', cmpListmbp{:});
%R_VAR_MLP       =  sprintf('%s ', cmpListmlp{:});
%R_VAR_MLPs      =  sprintf('%s ', cmpListmlps{:});
%R_VAR_MLPr      =  sprintf('%s ', cmpListmlpr{:});
R_VAR_ERP       =  sprintf('%s ', cmpList{:} );
%R_VAR_ERPs      =  sprintf('%s ', cmpLists{:} );
R_VAR_ERPr      =  sprintf('%s ', cmpListr{:} );
R_VAR_NOISE     = 'p2p powEOG \n';


fprintf(fid,R_VAR_TONE);
%fprintf(fid,R_VAR_MBP);
%fprintf(fid,R_VAR_MLP);
%fprintf(fid,R_VAR_MLPs);
%fprintf(fid,R_VAR_MLPr);
fprintf(fid,R_VAR_ERP);
%fprintf(fid,R_VAR_ERPs);
fprintf(fid,R_VAR_ERPr);
fprintf(fid,R_VAR_NOISE)


fprintf(fid2,R_VAR_TONE);
%fprintf(fid2,R_VAR_MBP);
%fprintf(fid2,R_VAR_MLP);
%fprintf(fid2,R_VAR_MLPs);
%fprintf(fid2,R_VAR_MLPr);
fprintf(fid2,R_VAR_ERP);
%fprintf(fid2,R_VAR_ERPs);
fprintf(fid2,R_VAR_ERPr);
fprintf(fid2,R_VAR_NOISE)


bl = ' ';
for tr=1:length(tone.trialNr)
    
    cur_line = [                                                             ...
        num2str(tone.trialNr(tr))                               bl   ...
        num2str(tone.seqNr(tr))                                 bl   ...
        num2str(tone.ampInd(tr))                                 bl   ...
        num2str(trial.type(tone.trialNr(tr)))                    bl   ...
        num2str(tone.type(tr))                                  bl   ...
        num2str(tone.ISI(tr))                                   bl   ...
        num2str(tone.deltaPitch(tr))                            bl   ...
        num2str(tone.time(tr))                                  bl   ...
        num2str(mn(tr,1))                                       bl   ...
        num2str(mn(tr,2))                                       bl   ...
        num2str(mn(tr,3))                                       bl   ...
        num2str(mn(tr,4))                                       bl   ...
        num2str(mn(tr,5))                                       bl   ...
        num2str(mn(tr,6))                                       bl   ...
        num2str(mn(tr,7))                                       bl   ...
        num2str(mn(tr,8))                                       bl   ...
        num2str(mn(tr,9))                                       bl   ...
        num2str(mnraw(tr,1))                                       bl   ...
        num2str(mnraw(tr,2))                                       bl   ...
        num2str(mnraw(tr,3))                                       bl   ...
        num2str(mnraw(tr,4))                                       bl   ...
        num2str(mnraw(tr,5))                                       bl   ...
        num2str(mnraw(tr,6))                                       bl   ...
        num2str(mnraw(tr,7))                                       bl   ...
        num2str(mnraw(tr,8))                                       bl   ...
        num2str(mnraw(tr,9))                                       bl   ...
        num2str(p2p(tr))                                       bl   ...
        num2str(powEOG(tr) )                                   '\n'];
     
    fprintf(fid, cur_line);
    fprintf(fid2, cur_line);
    
end
fclose(fid)
fclose(fid2)


'done'
end
