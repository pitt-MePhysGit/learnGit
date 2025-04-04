
%% --------
clear all
clear all

disc = resample( 1:10,1,2);
dir = '/Users/teichert/Dropbox/ampOdd_click/';
[expList b c d e f g h j k permission] = textread([dir,'ampOddclick_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s', 400);

Nexp = length(expList)-1;
subjList = [];
subjVector = zeros(length(expList),1);
for i = 2:(Nexp+1)
    tmp = strsplit(expList{i},'_');
    subjList{i} = tmp(1);
    
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
    
end
    
%ind = [16 32];
%ind = [32];
% check out expList(4) seems weird on ch 6 or 7
%ind = [21 22 23 25 29 32];
ind = find(  (strcmp(permission,'k') & subjVector == 3));%| (strcmp(permission,'l'))))%% & subjVector>0) )    ; %|  strcmp(permission,'l') | strcmp(permission,'s'));
%ind =  find (strcmp(permission,'m') );

%ind = [290];
%ind = ind( ind>69)
%ind = find( subjVector == 4)
%ind = [152 148]
%ind = find(subjVector==2 & strcmp(permission,'p') & (strcmp(k,'TB4')|strcmp(k,'manger') ) ); %%( strcmp(k,'ER1') | strcmp(k,'ER2')) );


expList(ind)

% trigger problems with Walt_20150813_0825
%%
for i = (length(expList(ind))):length(expList(ind))
    i
    close all
    try
    if ~(permission{ind(i)}=='t')
        %disc = prep_ampOddClick_function_old( expList(ind(i)) );
        %%disc = prep_RS_flex_causal_function( expList(ind(i)) );
    end
    end
    
    if (permission{ind(i)}=='t')
        'not allowed to run hypothesis testing data set at this point in time!!!'
    end
    
    %catch
       % disc = -1;
    %end
    
    tmp  = expList(  ind(i) )
    disc = ampOddClick2R( tmp{1}, 0);
    
end

'all done'
