%% --------
clc
clear all
clear df args cst chNr


addpath(genpath('~/Dropbox/RSkernel/m-files/'));

%directory         = 'ampOdd_click/';
directory         = 'ampOddClick/';
SR                = 1000;

%EEGLabDir         = '/Volumes/DATA test/'; %find_localDir();
%EEGLabDir         = '/Volumes/Drobo5D3/EEG/EEGLab/';
EEGLabDir         = '/Volumes/EEGLab/EEGLab/';
%baseDir           = find_baseDir();
baseDir            = '/Volumes/Drobo5D3/EEG/';


toDir             = 'spikeTriggered';

%load( [baseDir directory exp{1} '/args_' num2str(SR) '.mat']);

args.method    = 'nev';  % select either nev or ws (for BOSS or Wavesorter, respectively)
args.baseDir   = baseDir;
args.localDir  = EEGLabDir;
args.directory = directory;
args.toDir     = toDir;
%args.spkDir    = '~/Dropbox/RSkernel/spikeSorting/';
%args.spkDir    = '/Volumes/Drobo5D3/EEG/sortedZJ_ampOddClick/';
args.spkDir    = '/Volumes/Drobo5D3/EEG/sortedZJ_ampOddClick/';
args.RDir      = [args.baseDir args.toDir '/rda/'];
%

disc = resample( 1:10,1,2);
%dir = '/Volumes/Drobo5D3/EEG/sortedZJ_ampOddClick/';

[expList b c d e f g h j k l m n o Econfig channel cell q r fxtn t] = textread([args.spkDir,'ampOddClick_SUA_ZJ_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s');

stimArtefact = 0; % if 1, try to remove stimulation artefact in read_Intan_function.m

for (j = 1:length(o))
    valInd(j) = strcmp(o{j},'g') ;
end

ind = find(valInd);
expList(ind)

getSpikeShape = 1;
%%
clear j
for j = 1:length(ind)
    i = ind(j);
    close all
    
    args.exp      = expList{ i };
    
    chNr          = str2num( channel{ i });
    args.chNr     = chNr;
    args.cellExtn = cell{ i };
    
    tmp           = strsplit( cell{i}, '_' );
    args.thrsh    = tmp{2};
    args.Namp     = 33+24;
    args.Nana     = 1;
    args.Ndig     = 1;
    args.SRamp    = SR;
    args.fxtn     = fxtn{ i };
    args.getSpikeShape = getSpikeShape;
    
    
    try
        disc               = spikes2Rbase( args );
    catch
    end
end



