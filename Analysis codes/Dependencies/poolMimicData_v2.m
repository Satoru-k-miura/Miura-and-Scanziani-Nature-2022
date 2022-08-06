function S=poolMimicData_v2

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

for i=1:numFiles
    C=load(fullfile(fpath,fnames{i}));
    D = C.dat;

    % Make sure nasal is negative
    DNasal = IndexedStructCopy(D, strcmp({D(:).direction},'NasalMimic'));
    
    numEachDir(i,1) = length(DNasal);
    numEachDir(i,2) = length(D)-length(DNasal);
    
    ampNasal = mean([DNasal(:).amplitude]);
    if ampNasal>=0
        for ii=1:length(D)
            if strcmp(D(ii).direction,'NasalMimic')
                D(ii).amplitude = -D(ii).amplitude;
            end
        end
    end
    
    [~,sortOrder] = sort([D(:).amplitude]);
    D = D(sortOrder); % sort by amplitude

    E{i} = D;
end

finalN = min(numEachDir(:,1)); % minimum number of nasal saccades
finalT = min(numEachDir(:,2));

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
        dat(i).direction = 'NasalMimic';
    else
        dat(i).direction = 'TemporalMimic';
    end
    
    dat(i).originalAmplitudes = cA;
    
    S.dat = dat;
    S.params = C.params;
end






