function [disc] = call_AOCFcn_tinyone(subjINP,expINP,skipINP)

onsetArtefact = 1;
stimArtefact  = 0;
doSpikes      = 0;
SRamp         = 1000;

%
switch subjINP
    case 'Jesse'
        subjNr = 1;
    case 'Rockey'
        subjNr = 2;
    case 'Walter'
        subjNr = 3;
    case 'Sam'
        subjNr = 4;
    case 'Ben'
        subjNr = 5;
end
subjNr
%subjINP  = str2num(subjINP)
subjINP
expINP
skipINP = str2num(skipINP)

% =================== fix the annoying octave library bug
try
    disc = resample( 1:10,1,2);
end

directory         = 'ampOddClick';
dir = ['/Users/teichert/Dropbox/' directory '/'];

addpath(genpath('~/Dropbox/toolbox/'));
addpath(genpath('~/Dropbox/toolbox/testIntan/'));
addpath(genpath('~/Dropbox/RSkernel/m-files/'));
addpath(genpath('~/Dropbox/ampOddClick/m-files/'));

[expList b c d e f g h j k l permission Econfig preproc] = textread([dir,'ampOddclickdual_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 400);

Nexp           = length(expList)-1;
subjList       = [];
dateList       = [];
subjVector     = zeros(length(expList),1);
preprocFlag(1) = -1;

subjVector = zeros(length(expList),1);



for i = 2:(Nexp+1)
    tmp = strsplit(expList{i},'_');
    subjList{i} = tmp(1);
    dateList{i} = tmp(2);
    
    
    if strcmp( tmp{1}, 'Jesse')
        subjVector(i) = 1;
    end
    
    if strcmp( tmp{1}, 'Rockey')
        subjVector(i) = 2;
    end
    
    if strcmp( tmp{1}, 'Walter')
        subjVector(i) = 3;
    end
    
    if strcmp( tmp{1}, 'Sam')
        subjVector(i) = 4;
    end
    
    preprocFlag(i) = str2num(preproc{i});
end


subjVector
ind = find(  strcmp(permission,expINP) & subjVector == subjNr);
expList(ind)

%%
for i = (length(ind)+1-skipINP):(length(ind)+1-skipINP)
    i
    close all
    if ~(permission{ind(i)}=='x')
        
        rawDataFolder = [ '/Volumes/rawData/EEG/', subjList{ind(i)}{1} '/' subjList{ind(i)}{1} '_' dateList{ind(i)}{1}  '/' expList{ind(i)} '/']
        spellCheck(i) = exist( rawDataFolder , 'dir');
        args          = prep_ampOddClick_function( expList(ind(i)), Econfig{ind(i)}, SRamp, 1 );
        
        tmp  = expList(  ind(i) )
        disc = ampOddClick2R( tmp{1}, Econfig(ind(i)) ,SRamp);
        
        
        %args = prep_FFR_function( expList{ind(i)}, Econfig{ind(i)}, SRamp, doSpikes, toneArtefact(ind(i)), onsetArtefact(ind(i)), NFreq(ind(i)), triggerTime(ind(i)) );
        %load( ['/Volumes/Drobo5D3/EEG/EEGLab/' directory '/' expList{ind(i)} '/args_' num2str(SRamp) '.mat'] );args = object;clear object
        %disc = FFR2R( args );
    end
end
disc = -1;
'all done'
