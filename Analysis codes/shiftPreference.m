function shiftPreference(varargin)

% calculates aroc for the direction of full-field image shift

%%
% kernel SD used for GPFA
kernSD = 20;

before = -500;
after = 500;
condVS = {'condition1'};
sortBy = 'VS'; % colors to separate different events alt: 'direction'

assignopts(who, varargin);

%% New file?
answer = questdlg('Use existing data set or create a new data set?', ...
    'Create new data set?', ...
    'Create new','Use existing','Use existing');
% Handle response
switch answer
case 'Create new'
    newFile = 1;
case 'Use existing'
    newFile = 0;
end

%% If new file, PCA or GPFA
if newFile
%     answer = questdlg('PCA on trial average or GPFA?', ...
%         'Analysis type', ...
%         'PCA','GPFA','GPFA');
%     % Handle response
%     switch answer
%     case 'PCA'
%         pcagpfa = 'pca';
%         % if PCA, option to choose multiple experiments
%         answer2 = inputdlg('How many experiments do you want to concatenate?');
%         numExp = str2double(answer2{1});
%     case 'GPFA'
%         pcagpfa = 'gpfa';
%         numExp = 1;
%     end
    numExp = 1;
    % for the number of experiments, create a data structure, and then concatenate
    for i=1:numExp        
%         result = dataStructForGPFA_saccades(kernSD,'before',before,'after',after,'condVS',condVS,'sortBy',sortBy); % condition is based on blank trials/beh file
%         result = dataStructForGPFA_saccades(kernSD,'before',before,'after',after,'condVS',condVS,'sortBy',sortBy); % condition is based on VS channel
%         result = dataStructForGPFA_shiftLoc(kernSD,'before',before,'after',after,'condVS',condVS,'sortBy',sortBy); % condition is based on VS channel
        result = dataStructForGPFA_shiftVariable(kernSD,'before',before,'after',after,'condVS',condVS,'sortBy',sortBy);
        
        D = result.dat;
    end
    
% If old file, upload
else 
    while 1
        [fname,fpath]=uigetfile('_DataHighFormat.mat','Choose existing data structure with raw spike series for analysis','MultiSelect','off');
        if fpath == 0
            error('Execution cancelled');
        else
            if contains(fname,'DataHighFormat')
                D = load(fullfile(fpath,fname));
                D = D.dat;
                break
            end
        end
    end
end

%% Copy sortBy data to field 'condition', change colors accordingly
[D.condition] = D.(sortBy);
conds = {D.condition};
condTypes = unique(conds);
mag = find(strcmp(conds,condTypes{1}));
for i=1:length(D)
    if ismember(i,mag)
        D(i).epochColors = [0.8 0 0.8];
    else
        D(i).epochColors = [0 0.8 0.8];
    end
end

%% Select relevant subsets from D
c = unique({D(:).VS});
if length(c)>1
    % CHOOSE condition to analyze
    answer = questdlg('Which VS condition do you want to analyze?', ...
    'Choose VS condition', ...
    c{1},c{2},'Both','Both');
    if ~strcmp(answer,'Both')
        D = IndexedStructCopy(D, strcmp({D(:).VS},answer));
    end
end

d = unique({D(:).direction});
if length(d)<1
    error('no directions to discriminate')
end

%%
% Create two matrices for each direction. SU x time x events
numBins = abs(before) + after + kernSD;
numSUs = size([D(1).data],1);
numEvents1 = sum(strcmp({D(:).direction},'NasalMimic'));
numEvents2 = sum(strcmp({D(:).direction},'TemporalMimic'));
Dir1 = zeros(numSUs,numBins,numEvents1);
Dir2 = zeros(numSUs,numBins,numEvents2);
Dir3 = zeros(numSUs,numBins,numEvents1+numEvents2);
cm1 = 0;
cm2 = 0;
cm3 = 0;
for i=1:length(D)
    if strcmp({D(i).direction},'NasalMimic')
        cm1 = cm1+1;
        data = D(i).data;
        Dir1(:,:,cm1) = data;
    else
        cm2 = cm2+1;
        data = D(i).data;
        Dir2(:,:,cm2) = data;
    end
    cm3 = cm3+1;
    Dir3(:,:,cm3) = data;
end

% for each SU, calculate the average number of spikes and the std dev for
% each direction
% winWidth = 200;
% midT = ceil(numBins/2);
% spkNasal = squeeze(sum(Dir1(:,midT-winWidth:midT,:),2));
% baseNasal = squeeze(sum(Dir1(:,midT-201-winWidth:midT-201,:),2));
% 
% spkTemporal = squeeze(sum(Dir2(:,midT-winWidth:midT,:),2));
% baseTemporal = squeeze(sum(Dir2(:,midT-201-winWidth:midT-201,:),2));

winWidth = 100;
midT = ceil(numBins/2);
spkNasal = squeeze(sum(Dir1(:,midT:midT+winWidth,:),2));
baseNasal = squeeze(sum(Dir1(:,midT-201-winWidth:midT-201,:),2));

spkTemporal = squeeze(sum(Dir2(:,midT:midT+winWidth,:),2));
baseTemporal = squeeze(sum(Dir2(:,midT-201-winWidth:midT-201,:),2));


% Ignore units with extremely low baseline
bFRn = mean(baseNasal,2)*1e3/winWidth; % base line spike rate nasal
bFRt = mean(baseTemporal,2)*1e3/winWidth; % base line spike rate temporal

% tU = (bFRn >= 0.2) & (bFRt >= 0.2); % only units with higher than 0.2Hz baseline
tU = (bFRn >= 0) & (bFRt >= 0); % only units with higher than 0.2Hz baseline

numSUs = sum(tU); % update the number of SUs used for analysis

spkNasal = spkNasal(tU,:);
baseNasal = baseNasal(tU,:);
numSpkNasal = mean(spkNasal,2); 
numSpkNasalChange = numSpkNasal-mean(baseNasal,2);
stdSpkNasal = std(spkNasal,0,2);
spkTemporal = spkTemporal(tU,:);
baseTemporal = baseTemporal(tU,:);
numSpkTemporal = mean(spkTemporal,2); 
numSpkTemporalChange = numSpkTemporal-mean(baseTemporal,2);
stdSpkTemporal = std(spkTemporal,0,2);


changeSig = zeros(numSUs,2);
for i=1:numSUs
    p = ranksum(spkNasal(i,:),baseNasal(i,:));
    changeSig(i,1) = p;
    p = ranksum(spkTemporal(i,:),baseTemporal(i,:));
    changeSig(i,2) = p;
end

% baseMean = mean([baseNasal baseTemporal],2);
baseMean = [mean(baseNasal,2) mean(baseTemporal,2)];

sig = false(numSUs,1);
pDir = ones(numSUs,1);
for i=1:numSUs
    p=ranksum(spkNasal(i,:),spkTemporal(i,:));
%     p=ranksum(spkNasal(i,:)-mean(baseNasal(i,:)),spkTemporal(i,:)-mean(baseTemporal(i,:)));
%     p=ranksum(spkNasal(i,:)-baseNasal(i,:),spkTemporal(i,:)-baseTemporal(i,:));
    pDir(i) = p;
    if p<0.05
        sig(i) = 1;        
    end
end

binsize = 0.1;
histbins = -1:binsize:1;

foldChange = [numSpkNasal numSpkTemporal]./baseMean;

df = zeros(numSUs,1);
for i=1:numSUs
    if sum(spkTemporal(i,:))==0 && sum(spkNasal(i,:))==0
        df(i) = 0;
    else
        a = roc_curve(spkTemporal(i,:)',spkNasal(i,:)',0,0);
%         a = roc_curve(spkTemporal(i,:)'-mean(baseTemporal(i,:)),spkNasal(i,:)'-mean(baseNasal(i,:)),0,0);
%         a = roc_curve(spkTemporal(i,:)'-baseTemporal(i,:),spkNasal(i,:)'-baseNasal(i,:),0,0);
        df(i) = 2*a.param.AROC-1;
    end
end

nonsigHist = histcounts(df(~sig),histbins)'/length(df);
sigHist = histcounts(df(sig),histbins)'/length(df);
% plotting
figure;
x=histbins(1:end-1)+binsize/2;
b=bar(x,[sigHist nonsigHist],1,'stacked');
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.3010 0.7450 0.9330];
xlabel('Gini coeff')
ylabel('Frequency of occurence')

figure;
scatter(foldChange(~sig,1),foldChange(~sig,2)); hold on
scatter(foldChange(sig,1),foldChange(sig,2),'filled');
xlabel('Fold change over baseline, nasal saccades')
ylabel('Fold change over baseline, temporal saccades')

includedUnits = tU;
aroc = df;
save(['shiftPref_aroc_',num2str(winWidth)],'numSpkNasal','stdSpkNasal','numSpkNasalChange','numSpkTemporal','stdSpkTemporal','numSpkTemporalChange','aroc','foldChange','pDir','sig','changeSig','includedUnits')


end

function result = dataStructForGPFA_shiftLoc(kernSD,varargin)
% for VS shift experiments to mimic saccades

%% Obtain parameters
R = InputInitParams;

params = R{1};
eliminate_spikes = R{2};
% Bin size for PSTH in ms
params.binsize_hist = 1; 

before = -500; % ADDed to saccade time onset, subtract kernSD/2 such that a bin falls around 0
after = 500; % ADDed to saccade time onset

condVS = {'condition1'};

sortBy = 'direction'; % colors to separate different events alt: 'direction'

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
fullPath = tetpath(1:(strfind(tetpath,'_ephys')+5));
lastchar=strfind(fullPath,'_ephys')-1;
firstchar = strfind(fullPath,'\');
firstchar=firstchar(firstchar<lastchar);
firstchar = firstchar(end)+1;
exp_id = fullPath(firstchar:lastchar);


%% Determine trigger (trial) timing

% Specify file path
trigger_fname = strcat(exp_id,'_DIN00.mat');

% From DIN_00 channel, obtain trigger timings
dg00read = load(fullfile(fullPath,trigger_fname),'-mat');
din00 = dg00read.data;

% Get trial timings in samples
trialtiming = Get_timing_digital(din00,'both'); % first column, trial start, second column, trial end

% Convert to ms
trialtiming = trialtiming * params.sample_duration;

td = trialtiming(:,2)-trialtiming(:,1);
trialtiming = trialtiming(td>5,:);

clear dg00read din00;

%% Determine LED timings

% Specify file path
led_fname = strcat(exp_id,'_ADC00.mat');

if exist(fullfile(fullPath,led_fname),'file')
    % From ADC_00 channel, obtain led voltage
    an00read = load(fullfile(fullPath,led_fname),'-mat');
    an00 = an00read.v;

    % filter data
    an00 = medfilt1(an00);

    % Take the voltage divider into account
    an00 = an00/6.5*9.8; % 3.3kOhm and 6.5kOhm resistors

    thr=0.2;

    % Convert to T/F
    an00Dg = an00>thr;

    % Get led timings
    trgt = params.freq*0.0005;
    ledtiming = Get_timing_digital_loc(an00Dg','both',trgt);

    % Convert to ms
    ledtiming = ledtiming * params.sample_duration;
else
    ledtiming = [];
end

% 
clear din00 an00 an00Dg

%% Determine shift timings

% Specify file path
pulse_fname = strcat(exp_id,'_DIN03.mat');

% From DIN_00 channel, obtain trigger timings
dg03read = load(fullfile(fullPath,pulse_fname),'-mat');
din03 = dg03read.data;

% Get TTL timings in samples
shifttiming = Get_timing_digital(din03,'both'); % Column vector of onset timings

% Convert to ms
shifttiming = shifttiming * params.sample_duration;

clear dg03read din03;


%% determine which vs happened in which trials
% load shift direction file
fread = load('F:\MATLAB3\shift100x300','-mat');
shiftDirMat = fread.shiftDir;
clear fread

directions = [];
sT = [];
for i=1:size(trialtiming,1)
    np = (shifttiming(:,1) - trialtiming(i,1)).*(shifttiming(:,1) - trialtiming(i,2)) < 0;
    thistrial = find(np);
    numvs = length(thistrial);
    directions = [directions;shiftDirMat(i,1:numvs)'];
    sT = [sT;shifttiming(thistrial,:)];
end

% delete vstiming in case it's detected outside of trials... highly
% unlikely
shifttiming(isempty(directions))=[];
directions(isempty(directions))=[];

%% For every spike cluster, restructure and place in dat.spikes
dat=[];

fname = tetnames{1};
tetread = load(fullfile(tetpath,fname),'-mat');
clusClass = tetread.cluster_class; % index is a row vector of all spike times
fastvsreg = tetread.fastvsreg;
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
allVS=[];
sacVS=cell(0);
sacDirections = cell(0);
sacAmplitudes = cell(0);
sacPositions = cell(0);
eyeTraces = cell(0);
for i=1:length(directions)
    vsTime = sT(i,1);

    dat(end+1).spikes = template;
    dat(end).trialId = numEvents;
    sacVS(end+1) = condVS(1); 
    % direction
    direc = directions(i);
    if direc>0.5 % 1 is nasal mimic, 0 is temporal mimic
        sacDirections(end+1) = {'NasalMimic'};
        % amplitudes
        amp = 11;
        sacAmplitudes{end+1} = amp;
    else
        sacDirections(end+1) = {'TemporalMimic'};
        % amplitudes
        amp = 10;
        sacAmplitudes{end+1} = amp;
    end
    
    % positions
    pos = [];
    sacPositions{end+1} = pos;
    % eyetraces
    eyeT = [];
    eyeTraces{end+1} = eyeT;

    numEvents = numEvents+1;

    allVS = [allVS;vsTime];

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
    
    for j=1:size(allVS,1)
        vsTimeOnset=allVS(j,1);
        spks = spiketimes(spiketimes>=vsTimeOnset+before & spiketimes<=vsTimeOnset+after)-vsTimeOnset;
        spike01 = histc(spks,params.edges);
        spike01 = spike01(1:end-1);
        dat(j).spikes(l,:) = spike01;
        if strcmp(sacDirections{j},'NasalMimic') 
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
    dat(i).direction = sacDirections{i};
    dat(i).amplitude = sacAmplitudes{i};
    dat(i).position = sacPositions{i};
    dat(i).eyeTrace = eyeTraces{i};
    if strcmp(sortBy,'VS')
        if strcmp(sacVS{i},condVS{1})
            dat(i).epochColors = [0.8 0 0.8];
        else
            dat(i).epochColors = [0 0.8 0.8];
        end
    else
        if strcmp(sacDirections{i},'TemporalMimic')
            dat(i).epochColors = [0.8 0 0.8];
        else
            dat(i).epochColors = [0 0.8 0.8];
        end
    end
end

% rename for DataHigh
[dat.data]=dat.spikes;
dat=rmfield(dat,'spikes');

result.dat = dat;
result.datAve = datAve;
params.sortBy = sortBy;
% unitinfo.depth = depthum(goodUnits+1)';
% save file
answer = questdlg('Do you want to save the new data set?', ...
'Save new data set?', ...
'No','Yes','Yes');
% Handle response
switch answer
case 'Yes'
    savefilename = [exp_id '_DataHighFormat.mat'];
    currentFolder = pwd;
    cd('F:\MATLAB3\GPFA results');
    [savefile,savepath] = uiputfile(savefilename);
    cd(currentFolder);
    if savepath
%         save(fullfile(savepath,savefile),'dat','datAve','params','unitinfo','-mat');
        save(fullfile(savepath,savefile),'dat','datAve','params','-mat');
    end
end

end


function answer = InputInitParams

prompt = {'Ephys sampling freq:',...
    'Eye tracking frame rate:',...
    'Eye movement threshold (degrees):',...
    'Minimum magnitude of saccades to analyze:',...
    'Saccade speed threshold for detection (deg/s):',...
    'Eliminate artifact spikes (1/0)?',...
    'Magnitude of the artifact spikes to eliminate (deg):',...
    'Max number of artifact spikes to eliminate (ratio to the number of samples):'};
dlg_title = 'Initial parameters';
num_lines = 1;
defaultans = {'30000','200','3','3','150','1','3','0.1'};
R = inputdlg(prompt,dlg_title,num_lines,defaultans);

% Ephys sampling frequency samples/s
params.freq = str2num(R{1});
% ms/sample
params.sample_duration = 1/params.freq*1000;
% Frame rate of the eye tracking camera
params.fr = str2num(R{2});
% Eye movement threshold in degrees
params.mv_threshold = str2num(R{3});
% minimum magnitude of saccades to consider
params.mag = str2num(R{4});
% Saccade speed threshold for detection (degrees/s)
params.speed_thr = str2num(R{5});
% eliminate artifactual spikes in eye movement?
eliminate_spikes = str2num(R{6});
% magnitude of the spikes to eliminate in degrees
params.artifact_mag = str2num(R{7});
% maximum number of spikes to eliminate (ratio to the number of samples)
params.spk_elim_ratio = str2num(R{8});

answer = {params eliminate_spikes};

end

%%
function timing = Get_timing_digital_loc(dg_array,specifier,trgt)
% This function will obtain the onsets and offsets of digital pulses and
% return their indices. specifier specifies whether to return onset,
% offset, or both. In case of both, first column, onset, second column,
% offset.
% LOG
% 3/10/16 change vec2mat to reshape

breaks = [true logical(diff(dg_array))];
runStart = find([breaks true]);
elem = dg_array(breaks);
len = diff(runStart);

elem(len<trgt)=0;

dg_array = runlen(elem,len);

dgDiff = diff(dg_array);

ind = find(dgDiff); % return indices of 1s and -1s.

%temp = vec2mat(ind,2);
temp = reshape(ind,2,[])';
temp(:,1)=temp(:,1)+1;

if strcmp(specifier,'onset')
    timing = temp(:,1);
elseif strcmp(specifier,'offset')
    timing = temp(:,2);
elseif strcmp(specifier,'both')
    timing = temp;
else
    error('invalid specifier');
end
end
