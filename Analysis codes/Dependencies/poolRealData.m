function S=poolRealData

% 

%% choose files to pool (from the same directory)
while 1
    [fnames,fpath]=uigetfile('_DataHighFormat.mat','Choose DataHighFormat files','MultiSelect','on');
    if fpath == 0
        error('Execution cancelled');
    end
    if ~iscell(fnames)
        fnames=cellstr(fnames);
    end
    numFiles = length(fnames);
    nameMatch = false(1,numFiles);
    for i=1:numFiles
        tf = strfind(fnames{i},'_DataHighFormat');
        nameMatch(i) = isempty(tf);
    end
    if ~any(nameMatch)
        break;
    end
end

%% Load files
fnames = sort(fnames);

numEachDir = zeros(numFiles,2); % first column for the number of nasal, second temporal
numSU = zeros(numFiles,1);

for i=1:numFiles
    C=load(fullfile(fpath,fnames{i}));
    D = C.dat;
    D = IndexedStructCopy(D,strcmp({D(:).VS},'statVertical'));
    % Make sure nasal is negative
    DNasal = IndexedStructCopy(D, strcmp({D(:).direction},'Nasal'));
    
    numEachDir(i,1) = length(DNasal);
    numEachDir(i,2) = length(D)-length(DNasal);
    
    ampNasal = mean([DNasal(:).amplitude]);
    if ampNasal>=0
        for ii=1:length(D)
            if strcmp(D(ii).direction,'Nasal')
                D(ii).amplitude = -D(ii).amplitude;
            end
        end
    end
    
    [~,sortOrder] = sort([D(:).amplitude]);
    D = D(sortOrder); % sort by amplitude
    numSU(i) = size(D(1).data,1);
    E{i} = D;
end

finalN = min(numEachDir(:,1)); % minimum number of nasal saccades
finalT = min(numEachDir(:,2));
SUtally = sum(numSU);

% popup of summary
f = msgbox({['Number of nasal saccades: ' num2str(finalN)];...
    ['Number of temporal saccades: ' num2str(finalT)];...
    ['Number of single units: ' num2str(SUtally)]});
uiwait(f)

% trim each data to the same size
for i=1:length(E)
    D = E{i};
    num_N = numEachDir(i,1);
    
    idx1 = num_N-finalN+1;
    idx2 = num_N+finalT;
    
    D = D(idx1:idx2);
    E{i} = D;
end

dat = struct([]);
for i=1:finalN+finalT
    dat(i).VS = 'statVertical';
    dat(i).trialId = i;
    
    data = [];
    cA = zeros(1,length(E));
    for j = 1:length(E)
        datE = E{j};
        cA(j) = datE(i).amplitude;
        data = [data;datE(i).data];
    end
    dat(i).data = data;
    dat(i).amplitude = mean(cA);
    if mean(cA)<0
        dat(i).direction = 'Nasal';
    else
        dat(i).direction = 'Temporal';
    end
    
    dat(i).originalAmplitudes = cA;
    
    S.dat = dat;
    S.params = C.params;
end






