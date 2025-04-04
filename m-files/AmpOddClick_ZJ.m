%% --------
clear all
clear all

disc = resample( 1:10,1,2);
dir = '/Users/ZJ/Dropbox/ampOdd_click/';


addpath(genpath('~/Dropbox/toolbox/'));
addpath(genpath('~/Dropbox/toolbox/testIntan/'));
addpath(genpath('~/Dropbox/RSkernel/m-files/'));

[expList b c d e f g h j k l permission Econfig preproc] = textread([dir,'ampOddclickdual_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', 400);

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
ind = find(  (strcmp(permission,'g') & subjVector >= 4));%| (strcmp(permission,'l'))))%% & subjVector>0) )    ; %|  strcmp(permission,'l') | strcmp(permission,'s'));
%ind =  find (strcmp(permission,'m') );
ind = find(  preprocFlag == 0 )

%ind = [290];
%ind = ind( ind>69)
%ind = find( subjVector == 4)
%ind = [152 148]
%ind = find(subjVector==2 & strcmp(permission,'p') & (strcmp(k,'TB4')|strcmp(k,'manger') ) ); %%( strcmp(k,'ER1') | strcmp(k,'ER2')) );
%ind = 24
%ind = ind([2,4]);
expList(ind)