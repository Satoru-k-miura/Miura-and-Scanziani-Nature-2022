function result = dataStructForGPFA_saccades2(kernSD,varargin)
% conditions of saccades are based on VS timing instead of blank trials

%% Obtain parameters
R = InputInitParams;

params = R{1};
eliminate_spikes = R{2};
% Bin size for PSTH in ms
params.binsize_hist = 1; 

before = -500; % ADDed to saccade time onset, subtract kernSD/2 such that a bin falls around 0
after = 500; % ADDed to saccade time onset

condVS = {'condition1','condition2'};

sortBy = 'VS'; % colors to separate different events alt: 'direction'

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
trialtiming = trialtiming(td>1,:);

clear dg00read din00;

%% Select trials
x = size(trialtiming,1);
prompt = {[num2str(x) ' trials found. Specify trial range:']};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
trange = str2num(answer{1});

trialtiming = trialtiming(trange,:);

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

%% Determine TTL (video frame) timing

% Specify file path
pulse_fname = strcat(exp_id,'_DIN02.mat');

% From DIN_02 channel, obtain pulse timings
dg02read = load(fullfile(fullPath,pulse_fname),'-mat');
din02 = dg02read.data;

% Get TTL timings in samples
frametiming = Get_timing_digital(din02,'onset'); % Column vector of onset timings

% Convert to ms
frametiming = frametiming * params.sample_duration;

clear dg02read din02;

%% Generate a cell array of frame timing for each trial

% Initialize cell array. Each cell is a trial. Each cell contains a column
% vector of frame timings for that trial.
frames_per_trial = cell(size(trialtiming,1),1);

for i=1:size(trialtiming,1)-1
    % Obtain all the pulse timings (ms) for the current trial
    current_array = frametiming(frametiming >= trialtiming(i,1) & frametiming < trialtiming(i+1,1));
    % Place the array in the cell as a column vector
    frames_per_trial{i} = current_array';
end

% For the last trial, we take all the frame timings from the trial start to
% the end
current_array = frametiming(frametiming >= trialtiming(end,1));
frames_per_trial{end} = current_array';

clear current_array frametiming;

%% Determine eye positions during every trial. Create an array of azimuth and elevation coordinates with their timing info.

% Specify file path. Have the eye data in the ephys rig.
eye_fname = strcat(exp_id,'.eye');

lastchar = strfind(fullPath,'Ephys data');
eyePath = fullfile(fullPath(1:lastchar-1),'Eye data');
% Obtain information on the file
eye_fileinfo = dir(fullfile(eyePath,eye_fname));
num_samples = eye_fileinfo.bytes/4; % single = 4 bytes

% Open the file for read only
eye_fid = fopen(fullfile(eyePath,eye_fname),'r');

% Read in the data as singles
data_eye = fread(eye_fid, num_samples, 'single')';

% Close the file
fclose(eye_fid);

% The first four elements in data_eye are (1, "number of trials", "number
% of elements for eye position info",0). This is followed by the eye
% position in every trial. The first two elements for each trial are
% ("number of frames for the trial",3). These are followed by (x, y,
% radius) for every frame. The rest is not useful for this purpose.

% Number of trials
trials = data_eye(2);

% Starting element of the first trial
index = 5;

% Save the eye positions for each trial in a cell array
Eye_per_trial = cell(trials,1);

% Extract just the eye position information. Because the number of pulses
% does not match with the number of eye positions, we need to go trial by
% trial.
for i=1:trials
    % Number of caputured frames for this trial
    numframes_current = typecast(single(data_eye(index)),'int32');
    % There are 3 elements from each frame, (x, y, r)
    numelements_current = numframes_current*3;
    % Calculate the starting index for the next trial
    nextindex = index + 1 + numelements_current +1;
    
    % Extract the eye positions from data_eye for the current trial
    eye_pos_i = data_eye(index+2:nextindex-1);
    % Reshape the array into a matrix of three columns, x, y, and r
    eye_pos_i = reshape(eye_pos_i,3,numframes_current)';
    % Discard the first frame
    eye_pos_i = eye_pos_i(2:end,:);
    
    frame_array = frames_per_trial{i};
    
    % Because the number of pulses is more than the captured frames, we
    % discard the extra pulses at the end. The timing is saved in the 4th
    % column
    eye_pos_i(:,4) = frame_array(1:numframes_current-1)+1/params.fr*1000/2;
    
    % Place the matrix in eye_per_trial. Discard the first position, as it
    % is not an accurate frame. Transpose, such that the top row is
    % azimuth, the second row elevation, the third row radius, and the
    % fourth timing.
    Eye_per_trial{i} = eye_pos_i(2:end,:)';
    
    % Update index
    index = nextindex;
end

clear eye_fname eye_fileinfo eye_pos_i frame_array index nextindex numframes_current numelements_current;

%% Eliminate artifactual spikes in eye movement traces

if eliminate_spikes == 1
    
    Eye_per_trial = elim_eyespk2(trials,Eye_per_trial,params);
    
end

%% Obtain saccade timings 

% Create a cell array to save the saccade timings for each trial. Each cell
% contains a 2x? matrix, where the first column is the onset and the second
% column is the offset of each saccade.
Saccade_timing = cell(trials,5);
% the second column in the cell matrix was added for directions, amplitudes, and position 3/9/19

for i=1:trials
    eyeinfo = Eye_per_trial{i}; % take out the entire eye information for the trial
    eyetrace = eyeinfo(1,:); % Take out the azimuth trace

    answer = get_saccade_index(eyetrace,params); % obtain saccade indexes
    indexes = answer.indexes;

    if ~isempty(indexes)
        allindexes = [];
        for j=1:size(indexes,1)
            timing_i = [eyeinfo(4,indexes(j,1)) eyeinfo(4,indexes(j,2))];
            allindexes = [allindexes;timing_i];
        end
        % get rid of indexes that are too close to each other (the second
        % is usually not a real saccade)
        iSac = [false; diff(allindexes(:,1))<500];
        Saccade_timing{i,1} = allindexes(~iSac,:); % place allindexes in the cell
        Saccade_timing{i,2} = answer.directions(~iSac);
        Saccade_timing{i,3} = answer.amp(~iSac);
        Saccade_timing{i,4} = answer.pos(~iSac); 
        Saccade_timing{i,5} = answer.eyeTrace(~iSac,:); 
    end

    eyeinfo(1,:) = answer.eyetrace; % The cleaned eyetrace
    Eye_per_trial{i} = eyeinfo; % replace the old eyeinfo

end

clear eyeinfo eyetrace allindexes timing_i;

%% If there was LED, exclude saccades around LED
afterLED = 3000;
if ~isempty(ledtiming)
    ledtiming(:,2) = ledtiming(:,2)+afterLED;
    for i=1:size(Saccade_timing,1)
        saccadeTime = Saccade_timing{i,1};
        directions = Saccade_timing{i,2};
        amplitudes = Saccade_timing{i,3};
        positions = Saccade_timing{i,4};
        eyeTrace = Saccade_timing{i,5};
        whichRows = true(size(saccadeTime,1),1);
        for ii=1:size(saccadeTime,1)
            checkMat = (ledtiming-saccadeTime(ii,1))>0;
            checkN = sum(checkMat(:,1) - checkMat(:,2));
            if ~checkN==0
                whichRows(ii) = false;
            end
        end
        Saccade_timing{i,1} = saccadeTime(whichRows,:);
        Saccade_timing{i,2} = directions(whichRows);
        Saccade_timing{i,3} = amplitudes(whichRows);
        Saccade_timing{i,4} = positions(whichRows);
        Saccade_timing{i,5} = eyeTrace(whichRows,:);
    end
end

%% Determine visual stimulus timings

% Specify file path
pulse_fname = strcat(exp_id,'_DIN01.mat');

% From DIN_00 channel, obtain trigger timings
dg01read = load(fullfile(fullPath,pulse_fname),'-mat');
din01 = dg01read.data;

% Get TTL timings in samples
vstiming = Get_timing_digital(din01,'both'); % Column vector of onset timings

% Convert to ms
vstiming = vstiming * params.sample_duration;

clear dg01read din01;

% Select for vs timing within the trials
vstiming = vstiming(vstiming(:,1)>=trialtiming(1,1) & vstiming(:,2)<=trialtiming(end,2),:);

%% EXCLUDE saccades too close to VS switching

% for i=1:size(Saccade_timing,1)
%     saccadeTime = Saccade_timing{i,1};
%     directions = Saccade_timing{i,2};
%     amplitudes = Saccade_timing{i,3};
%     positions = Saccade_timing{i,4};
%     eyeTrace = Saccade_timing{i,5};
%     if ~isempty(saccadeTime)
%         validSac = ones(size(saccadeTime,1),1);
%         for ii=1:size(saccadeTime,1)
%             sacTimeOn = saccadeTime(ii,1);
%             % test whether it falls on VS or not
%             npTest = vstiming - ones(size(vstiming))*sacTimeOn;
%             if any(any(abs(npTest)<500))
%                 validSac(ii)=0;
%             end
%         end
%         Saccade_timing{i,1} = saccadeTime(logical(validSac),:);
%         Saccade_timing{i,2} = directions(logical(validSac));
%         Saccade_timing{i,3} = amplitudes(logical(validSac));
%         Saccade_timing{i,4} = positions(logical(validSac));
%         Saccade_timing{i,5} = eyeTrace(logical(validSac),:);
%     end
% end



%% Run a GUI to select saccades

% Saccade_timing = selectSaccades(Saccade_timing, Eye_per_trial, exp_id);
Saccade_timing = selectSaccades2(Saccade_timing, Eye_per_trial, exp_id); % Saccade_timing has more than the timing information, use selectSaccades2


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
sacDirections = cell(0);
sacAmplitudes = cell(0);
sacPositions = cell(0);
eyeTraces = cell(0);
for i=1:size(Saccade_timing,1)
    saccadeTime = Saccade_timing{i,1};
    directions = Saccade_timing{i,2};
    amplitudes = Saccade_timing{i,3};
    positions = Saccade_timing{i,4};
    eyeTrace = Saccade_timing{i,5};
    if ~isempty(saccadeTime)
        for ii=1:size(saccadeTime,1)
            dat(end+1).spikes = template;
            dat(end).trialId = numEvents;
            % saccade timing onset
            sacTimeOn = saccadeTime(ii,1);
            % test whether it falls on VS or not
            npTest = vstiming - ones(size(vstiming))*sacTimeOn;
            npTest = npTest(:,1).*npTest(:,2);
            if any(npTest<0)
                sacVS(end+1) = condVS(1); % falls in VS
            else
                sacVS(end+1) = condVS(2);
            end
            % direction
            direc = directions(ii);
            if direc>=0
                sacDirections(end+1) = {'Temporal'};
            else
                sacDirections(end+1) = {'Nasal'};
            end
            % amplitudes
            amp = amplitudes(ii);
            sacAmplitudes{end+1} = amp;
            % positions
            pos = positions(ii);
            sacPositions{end+1} = pos;
            % eyetraces
            eyeT = eyeTrace(ii,:);
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
        if strcmp(sacDirections{i},'Temporal')
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

