function dataStructForGPFA_free2D_v2(eyeinfo,varargin)

%% updates
% 8/2/2021 now using get_saccade_index2D

%%
params.speed_thr = 500;
params.mag = 3;
params.fr = 90;

condVS = {'free'};

% Bin size for PSTH in ms
params.binsize_hist = 1; 
before = -500; % ADDed to saccade time onset, subtract kernSD/2 such that a bin falls around 0
after = 500; % ADDed to saccade time onset
kernSD = 20;

assignopts(who, varargin);
params.before = before;
params.after = after;
params.kernSD = kernSD;
params.condNames = condVS;

before = before-kernSD/2;
after = after+kernSD/2;

% edges for PSTH
params.edges = before:params.binsize_hist:after; 
%% Get file path SU
while 1
    [tetnames,tetpath]=uigetfile('times_polytrodeX.mat','Choose times_polytrodeX.mat files','MultiSelect','off');
    if tetpath == 0
        error('Execution cancelled');
    end
    if ~iscell(tetnames)
        tetnames=cellstr(tetnames);
    end
    nameMatch = strfind(tetnames{1},'times_polytrode');

    if ~isempty(nameMatch)
        break;
    end
end

%% get saccade timing
answer = get_saccade_index2D(eyeinfo,params);
indexes = answer.indexes;
        
Saccade_timing = cell(1,7);

if ~isempty(indexes)
    
    allindexes = [];
    eyeTraceE = [];
    eyeTraceA = [];
    eyeTraceT = [];
    for j=1:size(indexes,1)
        timing_i = [eyeinfo(4,indexes(j,1)) eyeinfo(4,indexes(j,2))];
        allindexes = [allindexes;timing_i];
        eyeTraceA(j,:) = eyeinfo(1,indexes(j,1)-46:indexes(j,1)+46);
        eyeTraceE(j,:) = eyeinfo(2,indexes(j,1)-46:indexes(j,1)+46);
        eyeTraceT(j,:) = eyeinfo(4,indexes(j,1)-46:indexes(j,1)+46);
    end
    % get rid of indexes that are too close to each other (the second
    % is usually not a real saccade)
    iSac = [false; diff(allindexes(:,1))<100];
    Saccade_timing{1,1} = allindexes(~iSac,:); % place allindexes in the cell
    Saccade_timing{1,2} = nan;
    Saccade_timing{1,3} = answer.amp(~iSac);
    Saccade_timing{1,4} = answer.pos(~iSac); 
    Saccade_timing{1,5} = eyeTraceA(~iSac,:); 
    Saccade_timing{1,6} = eyeTraceE(~iSac,:); 
    Saccade_timing{1,7} = eyeTraceT(~iSac,:); 
end


%% For every spike cluster, restructure and place in dat.spikes
dat=[];

fname = tetnames{1};
tetread = load(fullfile(tetpath,fname),'-mat');
clusClass = tetread.cluster_class; % index is a row vector of all spike times
% fastvsreg = tetread.fastvsreg;
% depthum = tetread.depthAllum;
clear tetread

lastClass = max(clusClass(:,1)); % the largest class ID
recEnd = max(clusClass(:,2)); % approximate end time in ms of the recording

goodUnits = [];
% figure out good spiking units
for i=1:lastClass
    testCol = ones(size(clusClass,1),1)*i;

    % Figure out the rows that include the spike cluster
    % Take out the spikes times
    spiketimes = clusClass(clusClass(:,1)==testCol,2)';
    estFR = length(spiketimes)/recEnd*1e3;
    if estFR>=0 %&& strcmp(fastvsreg(i+1),'regular')
        goodUnits = [goodUnits i];
    end
end

% differentiate between conditions
template = zeros(length(goodUnits),length(params.edges)-1);

numEvents = 1;
allSaccades=[];
sacVS=cell(0);
sacTime = cell(0);
sacAmplitudes = cell(0);
sacPositions = cell(0);
eyeTraces = cell(0);
for i=1:size(Saccade_timing,1)
    saccadeTime = Saccade_timing{i,1};
%     directions = Saccade_timing{i,2};
    amplitudes = Saccade_timing{i,3};
    positions = Saccade_timing{i,4};
    eyeTrace = Saccade_timing{i,5};
    eyeTraceE = Saccade_timing{i,6};
    eyeTraceT = Saccade_timing{i,7};
    if ~isempty(saccadeTime)
        for ii=1:size(saccadeTime,1)
            dat(end+1).spikes = template;
            dat(end).trialId = numEvents;
            % saccade timing onset
            sacTimeOnOff = saccadeTime(ii,:);
            sacTime{end+1} = sacTimeOnOff;
            sacVS(end+1) = condVS(1);
            % amplitudes
            amp = amplitudes(ii);
            sacAmplitudes{end+1} = amp;
            % positions
            pos = positions(ii);
            sacPositions{end+1} = pos;
            % eyetraces
            eyeT = [eyeTrace(ii,:);eyeTraceE(ii,:);eyeTraceT(ii,:)];
            eyeTraces{end+1} = eyeT;

            numEvents = numEvents+1;
        end
        allSaccades = [allSaccades;saccadeTime];
    end
end

datAve(1).spikes = template;
datAve(1).trialId = 1;
datAve(2).spikes = template;
datAve(2).trialId = 2;


for l=1:length(goodUnits)
    testCol = ones(size(clusClass,1),1)*goodUnits(l);

    % Figure out the rows that include the spike cluster
    % Take out the spikes times
    spiketimes = clusClass(clusClass(:,1)==testCol,2)';
    
    numEvents1 = 0;
    numEvents2 = 0;
    spikeSum1 = zeros(1,length(params.edges)-1);
    spikeSum2 = zeros(1,length(params.edges)-1);
    
    for j=1:size(allSaccades,1)
        saccadeTimeOnset=allSaccades(j,1);
        spks = spiketimes(spiketimes>=saccadeTimeOnset+before & spiketimes<=saccadeTimeOnset+after)-saccadeTimeOnset;
        spike01 = histc(spks,params.edges);
        spike01 = spike01(1:end-1);
        dat(j).spikes(l,:) = spike01;
        if strcmp(sacVS{j},condVS{1}) 
            spikeSum1 = spikeSum1+spike01;
            numEvents1 = numEvents1+1;
        else
            spikeSum2 = spikeSum2+spike01;
            numEvents2 = numEvents2+1;
        end
    end
    spikeAve1 = spikeSum1/numEvents1;
    spikeAve2 = spikeSum2/numEvents2;
    
    datAve(1).spikes(l,:) = spikeAve1;
    datAve(2).spikes(l,:) = spikeAve2;
end

for i=1:length(sacVS)
    dat(i).VS = sacVS{i};
    dat(i).time = sacTime{i};
    dat(i).amplitude = sacAmplitudes{i};
    dat(i).position = sacPositions{i};
    dat(i).eyeTrace = eyeTraces{i};
end

% rename for DataHigh
[dat.data]=dat.spikes;
dat=rmfield(dat,'spikes');

result.dat = dat;
result.datAve = datAve;
% unitinfo.depth = depthum(goodUnits+1)';
% save file
answer = questdlg('Do you want to save the new data set?', ...
'Save new data set?', ...
'No','Yes','Yes');
% Handle response
switch answer
case 'Yes'
    savefilename = ['21-X_free_dX' '_2D_DataHighFormat.mat'];
    [savefile,savepath] = uiputfile(savefilename);
    if savepath
%         save(fullfile(savepath,savefile),'dat','datAve','params','unitinfo','-mat');
        save(fullfile(savepath,savefile),'dat','datAve','params','-mat');
    end
end

end


