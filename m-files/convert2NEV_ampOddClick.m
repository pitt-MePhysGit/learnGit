%% --------
clc
clear all
clear all

disc = resample( 1:10,1,2);
%dir = '/Users/teichert/Dropbox/ampOdd_click/';
ddir = '~/Dropbox/ampOddClick/';

%[expList b c d e f g h j k l permission Econfig preproc] = textread([dir,'ampOddclickdual_link.txt'], '%s %s %s %s %s %s %s %s %s %s %s %s %s %s', -1);

%fid = fopen([dir, 'ampOddclickdual_link.txt']);
%dat = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s');
%expList1 = dat{1};
%permission1 = dat{12};
%Econfig1 = dat{13};

fid = fopen([ddir, 'ampOddclick_link.txt']);
dat = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s');
expList1 = dat{1};
permission1 = dat{11};
Econfig1 = dat{12};

valsess = cellfun(@(x) numel(regexp(x, 'z'))>0, expList1);
expList1 = expList1(valsess, 1);
permission1 = permission1(valsess, 1);
Econfig1 = Econfig1(valsess, 1);


fid = fopen([ddir, 'ampOddclickVProbe_link.txt']);
dat = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s');
expList2 = dat{1};
permission2 = dat{11};
Econfig2 = dat{12};

fclose('all');

expList = cat(1, expList1, expList2{2:end, 1});
permission = cat(1, permission1, permission2{2:end, 1});
Econfig = cat(1, Econfig1, Econfig2{2:end, 1});

Nexp = length(expList)-1;
%extn   = 'm25';
%WSdir  = '/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/';
WSdir  = '/Volumes/EEGLab/EEGLab/ampOddClick/';
%NEVdir = '/Volumes/Drobo5D3/EEG/EEGLab/ampOddClick/';
NEVdir = '/Volumes/EEGLab/EEGLab/ampOddClick/';

newNEVdir = '/Volumes/EEGLab/EEGLab/ampOddClick/';
newprocdir = '/Volumes/Drobo5D3/EEG/ampOddClick/rda/';

fileset = cell(size(expList, 1), 4);
fileset{1, 1} = 'session';

for i = 2:(Nexp+1)
    %check if preprocessed
    fileset{i, 1} = expList{i, 1};
    tmpproc = dir([newprocdir, expList{i, 1}, '/', '*', '.mat']);    
    fileset{i, 2} = 0;
    if size(tmpproc, 1)>30
        fileset{i, 2} = 1;
    end
        
    %then check if all_m30-type text files exist
    tmptxt = dir([newNEVdir, expList{i, 1}, '/ch_1/all_', '*', '.txt']);
    fileset{i, 3} = tmptxt;
    
    %then if nevs exist
    tmpnev = dir([newNEVdir, expList{i, 1}, '/ch_1/', '*', '.nev']);
    fileset{i, 4} = tmpnev;
end


overwrite = 1;

%NEVdir = '/Volumes/DATA test/EEG/EEGLab/RS/';
%NEVdir = '~/teichert/Dropbox/RSkernel/spikeSorting/';

SR = 1000;

subjList = [];
subjVector = zeros(length(expList),1);
%preprocFlag(1) = 1;

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
    
    if strcmp( tmp{1}, 'Cirque')
    subjVector(i) = 5;
    end

%    preprocFlag(i) = str2num(preproc{i});
end
stimArtefact = 0; % if 1, try to remove stimulation artefact in read_Intan_function.m

% check out expList(4) seems weird on ch 6 or 7
%ind = [21 22 23 25 29 32];
%ind = find(  (strcmp(permission,'g') & subjVector > 2 & ~strcmp(Econfig,'ed')));%| (strcmp(permission,'l'))))%% & subjVector>0) )    ; %|  strcmp(permission,'l') | strcmp(permission,'s'));
%ind = find(  subjVector >= 3 & strcmp(permission,'n')& (strcmp(Econfig,'ed') | strcmp(Econfig,'edd')) );
%ind = find(  preprocFlag == 0 );

%%data prior to 20180103 for Sam and Walter isn't on the current set of
%%disks in the drobos
dataExist = zeros(size(fileset, 1), 1);
for i = 2:size(fileset, 1)
    tmp = strsplit(fileset{i, 1}, '_');
    dataExist(i, 1) = str2num(tmp{2})>20180103;
end


%different from overwrite?
%just checking for any all_*.txt files that haven't been converted
ind = find(cellfun(@(x) size(x, 1)>0, fileset(:, 3)) & cellfun(@(x) size(x, 1)==0, fileset(:, 4)) & dataExist);
%this doesn't mean that all all_*.txt files for a session have been, only at least one

[~, ia, ~] = unique(expList(ind));
ind = ind(ia);

%ind = ind(16);
expList(ind)

%% ===================================== experiment loop
 %smpFile = '~/Dropbox/RSkernel/m-files/sampleData.nev';
smpFile = '~/Dropbox/toolbox/example.nev';
NEV     = openNEV(smpFile);
       
ind 
 
for xnd = 1:length(ind)
    
    % load args file
    argsFile = [WSdir expList{ind(xnd)} '/args_' num2str(SR) '.mat'];
    load(argsFile);
    %%not sure what 'object' is, doesn't seem to be in argsFile which has
    %%only the object args, passed on unchanged instead of below, actually
    %%sometimes it is called object
    if numel(who('object'))==1
        args = object;
    end
    
    tmp = split(args.electrodeLayout, '');
    eCount = 0;
    for i = 2:(numel(tmp)-1)
        einc = 0;
        switch tmp{i}
            case 'v'
                einc = 24;
            case 't'
                einc = 32;
            case 'd'
                einc = 24;
        end
        eCount = eCount + einc;
    end
    
    for chnd = 1:eCount
%    for (chnd = 1:length(args.spkChan))
        %%
        clear NEW a wav tm mxWav mnWav p2p valInd
        
        allfNames = dir([NEVdir expList{ind(xnd)} '/ch_' num2str(chnd) '/all_*.txt']);
        allextn = {};
        for i = 1:numel(allfNames)
            tmp1 = strsplit(allfNames(i).name, '_');
            tmp2 = strsplit(tmp1{2}, '.');
            allextn{i} = tmp2{1};
        end
        
        for threshi = 1:numel(allextn)
            extn = allextn{threshi};            

            NEVfile = [NEVdir expList{ind(xnd)} '/ch_' num2str(chnd) '/' extn '.nev'];
            NEVfile

            if (~exist(NEVfile,'file') || overwrite)
                %% load example wave-forms
                spkFile = [WSdir expList{ind(xnd)} '/ch_' num2str(chnd) '/all_' extn '.txt']
                %spkFile = '~/Dropbox/RSkernel/spikeSorting/Sam_20160425_1220_zm9000/ch_15/all_m30.txt';

                [Nspk Ntst]    = textread(spkFile, '%n %n', 1);
                a              = textread(spkFile, '%n', -1, 'headerlines', 1);
                a              = reshape( a, [Ntst+1, Nspk]);
                wav            = a( 2:(Ntst+1),:)';

                timeFile = [WSdir expList{ind(xnd)} '/ch_' num2str(chnd) '/time_' extn '.txt'];

                tm       = textread(timeFile, '%n', -1);


                %% ================================== cut the outliers

                mxWav = max(wav,[],2);
                mnWav = min(wav,[],2);
                p2p   = mxWav - mnWav;

                valInd = find(p2p<400);

                % to work with BOSS, some spikes got pruned in convert2NEV
                % hence the spike time file is not accurate and needs to be 
                % re-written.

                spikeTimeFile = [NEVdir expList{ind(xnd)} '/ch_' num2str(chnd) '/time_' extn '_nev.txt'];

                fid     = fopen(spikeTimeFile, 'w');
                bl      = ' ';

                for eind = 1:(length(valInd))
                    cur_line = [ num2str(tm(eind)) '\n'];
                    fprintf(fid, cur_line);
                end
                fclose(fid);

                %%
                NEW = NEV;

                NEW.MetaTags.FilePath = [NEVdir expList{ind(xnd)} '/ch_' num2str(chnd)];
                NEW.MetaTags.Filename  = extn;

                %NEV.Data.Spikes
                %NEW.ElectrodesInfo = NEV.ElectrodesInfo(1);
                %NEW.MetaTags.ChannelID = NEV.MetaTags.ChannelID(1);

                NEW.Data.Spikes.Electrode  = ones(1,length(tm(valInd)) );
                NEW.Data.Spikes.TimeStamp  = tm(valInd)'*30000;
                NEW.Data.Spikes.Unit       = ones(1,length(tm(valInd)) );
                NEW.Data.Spikes.Waveform   = wav(valInd,8:55)';

    %             NEW.Data.Tracking = [];
    %             
    %             %NEW.Data.Comments.TimeStamp = [];
    %             NEW = rmfield( NEW, 'ObjTrackInfo');
    %             NEW = rmfield( NEW, 'VideoSyncInfo');
    %             NEW = rmfield( NEW, 'IOLabels');
    %             
    %             NEW.Data.SerialDigitalIO.TimeStamp = [];
    %             NEW.Data.Comments.TimeStamp        = [];
    %             NEW.Data.VideoSync.TimeStamp       = [];
    %             NEW.Data.Tracking                  = [];
    %             NEW.Data.PatientTrigger.TimeStamp  = [];
    %             
    %             NEW.Data.Comments.Color = ones(1,length(NEW.Data.Comments.CharSet) );
    %             
                saveNEV(NEW, [NEVdir expList{ind(xnd)} '/ch_' num2str(chnd) '/' extn '.nev'], 'noreport')
                %keyboard
                %tst = openNEV( [NEVdir expList{ind(xnd)} '/ch_' num2str(chnd) '/' extn '.nev']);
            end
        end
    end
end



%NEV     = openNEV('~/Dropbox/RSkernel/m-files/handSorted.nev'); 
%NEV     = openNEV('/Volumes/DATA/EEG/EEGLab/RS/Sam_20160425_1220_zm9000/ch_4/m30.nev'); 
%NEV     = openNEV('~/Dropbox/toolbox/example.nev'); 
%NEV     = oNEV('~/Dropbox/toolbox/example.nev');

%NEV     = openNEV('/Volumes/EEGLab/EEGLab/RSkernel/Sam_20160504_1445_zm11300/ch_4/m30.nev'); 


