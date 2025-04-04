
%% --------
clear all
clear all

try
    disc = resample( 1:10,1,2);
catch
end

dir = '/Users/teichert/Dropbox/ampOdd_click/';


addpath(genpath('~/Dropbox/toolbox/'));
addpath(genpath('~/Dropbox/toolbox/testIntan/'));
addpath(genpath('~/Dropbox/RSkernel/m-files/'));

[expList b c d e f g h j k l permission Econfig preproc] = textread([dir,'ampOddclickdual_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 400);
%[expList b c d e f g h j k l permission Econfig preproc] = textread([dir,'ampOddclickEEG_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 400);

Nexp = length(expList)-1;
subjList = [];
dateList = [];
subjVector = zeros(length(expList),1);
preprocFlag(1) = -1;


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

%ind = [16 32];
%ind = [32];
% check out expList(4) seems weird on ch 6 or 7
%ind = [21 22 23 25 29 32];
ind = find(  (strcmp(permission,'g') & subjVector == 3));%| (strcmp(permission,'l'))))%% & subjVector>0) )    ; %|  strcmp(permission,'l') | strcmp(permission,'s'));
%ind =  find (strcmp(permission,'m') );
%ind = find(  preprocFlag == 0 )

%ind = [290];
%ind = ind( ind>69)
%ind = find( subjVector == 4)
%ind = [152 148]
%ind = find(subjVector==2 & strcmp(permission,'p') & (strcmp(k,'TB4')|strcmp(k,'manger') ) ); %%( strcmp(k,'ER1') | strcmp(k,'ER2')) );
%ind = 24
%ind = ind(63);
expList(ind)

% trigger problems with Walt_20150813_0825

%%

numWorkers = min(1,length(ind));
poolobj = gcp('nocreate');
if numWorkers > 1
    poolobj = parpool('local',numWorkers, 'SpmdEnabled', false)
end
%for i = (length(ind)-27+1):(length(ind))
for i = (length(ind)):(length(ind))
%parfor (i = 1:length(expList(ind)), numWorkers)
    i
    close all
    %try
        if ~(permission{ind(i)}=='t')
            
            rawDataFolder = [ '/Volumes/rawData/EEG/', subjList{ind(i)}{1} '/' subjList{ind(i)}{1} '_' dateList{ind(i)}{1}  '/' expList{ind(i)} '/']
            spellCheck(i) = exist( rawDataFolder , 'dir');
            
            disc = prep_ampOddClick_function( expList(ind(i)), Econfig{ind(i)}, 1000, 0 );
            %%%disc = prep_RS_flex_causal_function( expList(ind(i)) );
        end
        
    %catch
    %    disc = -1;
    %end
    %end
    
    if (permission{ind(i)}=='t')
        'not allowed to run hypothesis testing data set at this point in time!!!'
    end
    
    %catch
    % disc = -1;
    %end
    %try
        tmp  = expList(  ind(i) )
        disc = ampOddClick2R( tmp{1}, Econfig(ind(i)) ,1000);
    %catch
        disc = -1;
    %end
    
end

'all done'
delete(poolobj)

spellCheck
