function saccadeFeatureRandomization_v4(varargin)
% Calculation of the cross-validated performance on pseudo-saccades after
% randomizing a certain number of top features... randomizing training set!

%% Params
% params for PCA to LDA
thrFR = 1;
cnvSD = 20; % Convolve with Gaussian, SD in ms, cnvSD=0 is no smoothing
spikeShape = 'both';
numLDALoops = 1;
numK = 2;   % K fold cross validation
shuffle = 0;
totalVariance = 75;
numRandSample = 5;
modelClass = 'LDA';
assignopts(who, varargin);

params.thrFR = thrFR;
params.cnvSD = cnvSD;
params.spikeShape = spikeShape;
params.numLDALoops = numLDALoops;
params.numK = numK;
params.shuffle = shuffle;
params.totalVariance = totalVariance; % min total variance accounted for in PCA
params.modelClass = modelClass;


poolobj = parpool('local',4);


allScore = [];
arocPerIter = [];
arocPvalPerIter = [];
elim_perc = 0:0.1:1;%[0 0.01:0.01:0.1 0.12 0.16 0.2 0.25 0.3 0.5 0.75 1];
allAcc = zeros(numRandSample,length(elim_perc));
AccFeaturePerm = zeros(numRandSample,length(elim_perc));
PCAexplainedVariance = zeros(numRandSample,length(elim_perc));

MimicD1 = poolMimicData_v2;
for ii=1:numRandSample
    disp(['Iteration ' num2str(ii) ' of ' num2str(numRandSample)])
    %% Generate balanced pooled dataset
    MimicD = resampleLoc(MimicD1);
    MimicD.dat = MimicD.dat(1:floor(floor(length(MimicD.dat)/numK)/2)*2*numK);
    [aroc_i,gini_p] = calculateAROC(MimicD);

    %% Perform PCA to LDA cross validation with permutation of each unit
    disp(' Calculating feature importance...')
    tic;
    params.unitCondition = [];
    [allacc,allaccPerm,unitsKeep,permPval] = PCAtoXDALoc(MimicD,params,ii);
    T = toc;
    disp([' Done. Elapsed time (s): ' num2str(T)])
    fullset = find(unitsKeep);
    score = nan(length(unitsKeep),1);
    allaccPerm = allacc - allaccPerm; %PFI for every k
    allaccPerm = nanmean(allaccPerm);
    for i=1:sum(unitsKeep)
        score(fullset(i)) = allaccPerm(i);
    end  

    H = false(length(unitsKeep),1);
    p = all(permPval<0.05,1);
    for i=1:sum(unitsKeep)
        H(fullset(i)) = p(i);
    end

    ScoreSort = sort(score,'descend');
    numNan = sum(isnan(score));
    disp(' Calculating accuracy with randomized features...')
    for i=1:length(elim_perc)
        unitCondition = true(length(score),1);
        p = elim_perc(i);
        disp(['  Randomizing top ' num2str(p*100) '% of all units'])
        if ~numNan && ~p
        else
            targetScore = ScoreSort(numNan + round((length(score)-numNan)*p));
            unitCondition(score>=targetScore) = false;
        end
        params.unitCondition = unitCondition; % Stable units, all others randomize

        [acc,accFeaturePerm,varExplained] = PCAtoXDALoc2(MimicD,params,ii);
        allAcc(ii,i) = acc;
        AccFeaturePerm(ii,i) = accFeaturePerm;
        PCAexplainedVariance(ii,i) = varExplained;
    end

    allScore = [allScore score];
    arocPerIter = [arocPerIter aroc_i];
    arocPvalPerIter = [arocPvalPerIter gini_p];
end
delete(poolobj);

%% save file
if strcmp(params.modelClass,'LDA')
    fname = 'LDAFeatureRandomization_v4';
elseif strcmp(params.modelClass,'QDA')
    fname = 'QDAFeatureRandomization_v4';
elseif strcmp(params.modelClass,'NN')
    fname = 'NNFeatureRandomization_v4';
end
save(fname,'allScore','allAcc','AccFeaturePerm','elim_perc','arocPerIter','arocPvalPerIter','PCAexplainedVariance','params','MimicD1');

end

function [allacc,allaccPerm,unitsKeep,permPval] = PCAtoXDALoc(MimicD,params,~)

thrFR = params.thrFR;
cnvSD = params.cnvSD;
spikeShape = params.spikeShape;

params.kernSD = MimicD.params.kernSD;
params.before = MimicD.params.before;
params.after = MimicD.params.after;

%%

D = MimicD.dat;

%% Select relevant conditions

% d = unique({D(:).direction});
d{1} = 'NasalMimic';
d{2} = 'TemporalMimic';


%% Copy sortBy data to field 'condition', change colors accordingly

for i=1:length(D)
    if strcmp(D(i).direction,d{1})
        D(i).condition = [d{1}];
        D(i).epochColors = [127/255 0 1];
    elseif strcmp(D(i).direction,d{2})
        D(i).condition = [d{2}];
        D(i).epochColors = [0 127/255 1];
    end
end

%% Remove low FR neurons and bin time points
allSpk = [D(:).data];
fr = mean(allSpk,2)*1e3;
unitsKeep = fr>=thrFR;
% 
% if ~isempty(unitCondition)
%     unitsKeep = unitsKeep & unitCondition;
% end

allSpk=cat(1,D(:).data)';
aS = reshape(allSpk,1020,[],length(D));
aS2 = squeeze(sum(aS(21:220,:,:))); % 200ms of baseline
fr = mean(aS2,2)/2;

aS3 = squeeze(sum(aS(511:610,:,:))); % 100ms post onset
% stddev = std(aS3,[],2);
% 
% unitsKeep = fr*10>=thrFR; %fr*10>=thrFR & stddev > 0.1;

nU = sum(unitsKeep);

for i=1:length(D)
%     t = aS3(unitsKeep,i);
%     t = t - fr(unitsKeep);
%     D(i).data = t;

    t = D(i).data(unitsKeep,:);
    t = t';
    t = reshape(t,20,[]);
    t = sum(t,1);
    t = reshape(t,[],nU)';
    if cnvSD > 0
        t = smoother(t,cnvSD,20);
    end
    D(i).data = t(:,29);

%     t = D(i).data(unitsKeep,:);
%     if cnvSD > 0
%         t = smoother(t,cnvSD,1);
%     end
%     D(i).data = mean(t(:,511:711),2);
end

%% Cross validate with permutation

[allacc,allaccPerm,permPval] = XDAonPCAcrossvalPerm(D,params);

end

%% Classification after elimination
function [acc,accFeaturePerm,varExplained] = PCAtoXDALoc2(MimicD,params,~)

thrFR = params.thrFR;
cnvSD = params.cnvSD;
unitCondition = params.unitCondition;
spikeShape = params.spikeShape;

if size(unitCondition,1)==1
    unitCondition = unitCondition'; % column vector
end

params.kernSD = MimicD.params.kernSD;
params.before = MimicD.params.before;
params.after = MimicD.params.after;

%%

D = MimicD.dat;

%% Select relevant conditions

% d = unique({D(:).direction});
d{1} = 'NasalMimic';
d{2} = 'TemporalMimic';


%% Copy sortBy data to field 'condition', change colors accordingly

for i=1:length(D)
    if strcmp(D(i).direction,d{1})
        D(i).condition = [d{1}];
        D(i).epochColors = [127/255 0 1];
    elseif strcmp(D(i).direction,d{2})
        D(i).condition = [d{2}];
        D(i).epochColors = [0 127/255 1];
    end
end

%% Remove low FR neurons, choose units in appropriate depths, and bin time points, smoothen if necessary
allSpk = [D(:).data];
fr = mean(allSpk,2)*1e3;
unitsKeep = fr>=thrFR;

allSpk=cat(1,D(:).data)';
aS = reshape(allSpk,1020,[],length(D));
aS2 = squeeze(sum(aS(21:220,:,:)));
fr = mean(aS2,2)/2;

aS3 = squeeze(sum(aS(511:610,:,:))); % 100ms post onset
% stddev = std(aS3,[],2);

% unitsKeep = fr*10>=thrFR; %fr*10>=thrFR & stddev > 0.1;
% 
% if ~isempty(unitCondition)
%     unitCondition = unitCondition(unitsKeep);
% end

params.unitCondition = unitCondition(unitsKeep)'; % row vector

nU = sum(unitsKeep);

for i=1:length(D)
%     t = aS3(unitsKeep,i);
%     t = t - fr(unitsKeep);
%     D(i).data = t;

    t = D(i).data(unitsKeep,:);
    t = t';
    t = reshape(t,20,[]);
    t = sum(t,1);
    t = reshape(t,[],nU)';
    if cnvSD > 0
        t = smoother(t,cnvSD,20);
    end
    D(i).data = t(:,29);

%     t = D(i).data(unitsKeep,:);
%     if cnvSD > 0
%         t = smoother(t,cnvSD,1);
%     end
%     D(i).data = mean(t(:,511:711),2);
end

%% Cross validate
[acc,accFeaturePerm,varExplained] = XDAonPCAcrossvalPerm2(D,params);

end

%% Cross validation with permutation of features
function [allacc,allaccPerm,permPval] = XDAonPCAcrossvalPerm(D,params)

numK = params.numK;
totalVariance = params.totalVariance;
modelClass = params.modelClass;

data=cat(1,D(:).data);
nU = size(D(1).data,1);

d = data; % concatenated data for the condition
d = reshape(d,nU,[]); % reshape such that nU x events
d = d'; % events x nU

% partition
cp = cvpartition(length(D),'KFold',numK);

cumN = 0; %track partition
allevents = 1:length(D);
allacc = zeros(numK,1);
allaccPerm = zeros(numK,nU);
permPval = zeros(numK,nU);
permNum = 20;

for k=1:numK
    testevents = cumN+1:cumN+cp.TestSize(k); % events to exclude//test
    cumN = cumN+cp.TestSize(k);
    trainevents = allevents(~ismember(allevents,testevents));
    
    d_train = d(trainevents,:);
    meand_train = mean(d_train);
    stdd_train = std(d_train);
    d_train = (d_train - meand_train)./stdd_train;
    
    d_train(isnan(d_train)) = 0;
    
    [coeff,score,~,~,explained,~] = pca(d_train);
    
    expTot = 0;
    dims = 0;
    while sum(expTot)<totalVariance && dims < length(trainevents)/5
        dims = dims+1;
        expTot = sum(explained(1:dims));
    end
    varExplained(1,1) = expTot;
    finalDims(1,1) = dims;

    % reduce
    dataPCA_train = score(:,1:dims);
    
    d_test = d(testevents,:);
    
    d_test = (d_test - meand_train)./stdd_train;
    dataPCA_test = d_test*coeff;
    dataPCA_test = dataPCA_test(:,1:dims);

    trainlabel = {D(trainevents).direction};
    
    if strcmp(modelClass,'LDA')
        Mdl = fitcdiscr(dataPCA_train,trainlabel,'DiscrimType','diaglinear','SaveMemory','off','FillCoeffs','off');
    elseif strcmp(modelClass,'QDA')
        Mdl = fitcdiscr(dataPCA_train,trainlabel,'DiscrimType','quadratic','SaveMemory','off','FillCoeffs','off'); 
    elseif strcmp(modelClass,'NN')
        Mdl = fitcnet(dataPCA_train,trainlabel); 
    else
        error('Specify model name')
    end

    truelabel = {D(testevents).direction};

    cerror = loss(Mdl,dataPCA_test,truelabel,'LossFun','classiferror');
    allacc(k) = 1-cerror;

    A = zeros(permNum,nU);
    for ii=1:permNum
        parfor i=1:nU
            d_perm = d_train;
            d_perm(:,i) = d_perm(randperm(size(d_train,1)),i);
            
            [coeff,score,~,~,explained,~] = pca(d_perm);
    
            expTot = 0;
            dims = 0;
            while sum(expTot)<totalVariance && dims < length(trainevents)/5
                dims = dims+1;
                expTot = sum(explained(1:dims));
            end
        
            % reduce
            dataPCA_perm = score(:,1:dims);
            
            if strcmp(modelClass,'LDA')
                Mdl = fitcdiscr(dataPCA_perm,trainlabel,'DiscrimType','diaglinear','SaveMemory','off','FillCoeffs','off');
            elseif strcmp(modelClass,'QDA')
                Mdl = fitcdiscr(dataPCA_perm,trainlabel,'DiscrimType','quadratic','SaveMemory','off','FillCoeffs','off'); 
            elseif strcmp(modelClass,'NN')
                Mdl = fitcnet(dataPCA_perm,trainlabel); 
            else
                error('Specify model name')
            end

            dataPCA_test = d_test*coeff;
            dataPCA_test = dataPCA_test(:,1:dims);

            cerror = loss(Mdl,dataPCA_test,truelabel,'LossFun','classiferror');

            A(ii,i) = 1-cerror;
        end
    end
    for i=1:nU
        [~,p] = ttest(A(:,i));
        permPval(k,i) = p;
    end

    A = mean(A);
    allaccPerm(k,:) = A;
end

end

%% Cross validation
function [acc,accFeaturePerm,varExplained] = XDAonPCAcrossvalPerm2(D,params)

modelClass = params.modelClass;
numK = params.numK;
totalVariance = params.totalVariance;
unitCondition = params.unitCondition; % stable units
unitsToRandomize = find(~unitCondition); 

data=cat(1,D(:).data);
nU = size(D(1).data,1);

d = data; % concatenated data for the condition
d = reshape(d,nU,[]); % reshape such that nU x events
d = d'; % events x nU

% partition
cp = cvpartition(length(D),'KFold',numK);

cumN = 0; %track partition
allevents = 1:length(D);
allacc = zeros(numK,1);
allaccPerm = zeros(numK,1);
permNum = 10;

for k=1:numK
    testevents = cumN+1:cumN+cp.TestSize(k); % events to exclude//test
    cumN = cumN+cp.TestSize(k);
    trainevents = allevents(~ismember(allevents,testevents));
    
    d_train = d(trainevents,:);
    meand_train = mean(d_train);
    stdd_train = std(d_train);
    d_train = (d_train - meand_train)./stdd_train;
    
    d_train(isnan(d_train)) = 0;
    
    [coeff,score,~,~,explained,~] = pca(d_train);
    
    expTot = 0;
    dims = 0;
    while sum(expTot)<totalVariance && dims < length(trainevents)/5
        dims = dims+1;
        expTot = sum(explained(1:dims));
    end
    varExplained(1,1) = expTot;
    finalDims(1,1) = dims;

    % reduce
    dataPCA_train = score(:,1:dims);
    trainlabel = {D(trainevents).direction};
    
    if strcmp(modelClass,'LDA')
        Mdl = fitcdiscr(dataPCA_train,trainlabel,'DiscrimType','diaglinear','SaveMemory','off','FillCoeffs','off');
    elseif strcmp(modelClass,'QDA')
        Mdl = fitcdiscr(dataPCA_train,trainlabel,'DiscrimType','quadratic','SaveMemory','off','FillCoeffs','off'); 
    elseif strcmp(modelClass,'NN')
        Mdl = fitcnet(dataPCA_train,trainlabel); 
    else
        error('Specify model name')
    end

    % Test
    d_test = d(testevents,:);
    d_test = (d_test - meand_train)./stdd_train;
    dataPCA_test = d_test*coeff;
    dataPCA_test = dataPCA_test(:,1:dims);
    
    truelabel = {D(testevents).direction};

    cerror = loss(Mdl,dataPCA_test,truelabel,'LossFun','classiferror');
    allacc(k) = 1-cerror;
    
    A = zeros(permNum,1);
    parfor ii=1:permNum
        d_perm = d_train;
        for i=1:length(unitsToRandomize)
            d_perm(:,unitsToRandomize(i)) = d_perm(randperm(size(d_train,1)),unitsToRandomize(i));
        end

        [coeff,score,~,~,explained,~] = pca(d_perm);
    
        expTot = 0;
        dims = 0;
        while sum(expTot)<totalVariance && dims < length(trainevents)/5
            dims = dims+1;
            expTot = sum(explained(1:dims));
        end

        % reduce
        dataPCA_perm = score(:,1:dims);

        if strcmp(modelClass,'LDA')
            Mdl = fitcdiscr(dataPCA_perm,trainlabel,'DiscrimType','diaglinear','SaveMemory','off','FillCoeffs','off');
        elseif strcmp(modelClass,'QDA')
            Mdl = fitcdiscr(dataPCA_perm,trainlabel,'DiscrimType','quadratic','SaveMemory','off','FillCoeffs','off'); 
        elseif strcmp(modelClass,'NN')
            Mdl = fitcnet(dataPCA_perm,trainlabel); 
        else
            error('Specify model name')
        end

        dataPCA_test = d_test*coeff;
        dataPCA_test = dataPCA_test(:,1:dims);

        cerror = loss(Mdl,dataPCA_test,truelabel,'LossFun','classiferror');

        A(ii) = 1-cerror;
    end
    allaccPerm(k) = mean(A);
end
    
acc = mean(allacc);
accFeaturePerm = mean(allaccPerm);
end

%% Mimic data resampling
function MimicD = resampleLoc(MD)

D = MD.dat;

% Check is nasal amplitude is negative, if not convert

for i=1:length(D)
    if strcmp(D(i).direction,'NasalMimic') && D(i).amplitude>0
        D(i).amplitude = -D(i).amplitude;
    end
end

DNasal = IndexedStructCopy(D, strcmp({D(:).direction},'NasalMimic'));

% ampNasal = mean([DNasal(:).amplitude]);
% if ampNasal>=0
%     for i=1:length(D)
%         if strcmp(D(i).direction,'NasalMimic')
%             D(i).amplitude = -D(i).amplitude;
%         end
%     end
% end

ampNasal = [DNasal(:).amplitude];

DTemporal = IndexedStructCopy(D, strcmp({D(:).direction},'TemporalMimic'));
R = randperm(length(DTemporal),length(DTemporal));
E = DTemporal(R);
ampTemp = [E(:).amplitude];
delta_amplitude = [];

ngflag = false(length(ampTemp),1);
for i=1:length(ampTemp)
    if ~isempty(DNasal)
        amp = ampTemp(i);
        [~,I] = min(abs(abs(ampNasal)-amp));
        delta = abs(ampNasal(I))-amp;
        if abs(delta)<4
            ngflag(i) = true;
            E = [E DNasal(I)];
            delta_amplitude(end+1) = delta;
            ngflag = [ngflag;true];
        else
        end
        %update DNasal to exclude I
        tr = ~((1:length(DNasal))==I);
        DNasal = IndexedStructCopy(DNasal,tr);
        ampNasal = [DNasal(:).amplitude];
    else
    end
        
end

E = E(ngflag);

l = 1:length(E);
so = randsample(l,length(l),'false');

E = E(so);

n = find(strcmp({E.direction},'NasalMimic'));
t = find(strcmp({E.direction},'TemporalMimic'));
nt = reshape([n;t],1,[]);

E = E(nt);

MimicD.dat=E;
MimicD.params = MD.params;

end

%% Calculate discriminability for each neuron
function [ginicoeff,gini_p] = calculateAROC(MimicD)

D = MimicD.dat;
hf = 1020/2+1;
allSpk=cat(1,D(:).data)';
aS2 = reshape(allSpk,1020,[],length(D));
baseSpkCount = squeeze(sum(aS2(21:220,:,:)))';
fr = mean(baseSpkCount)/2;
responseSpkCount = squeeze(sum(aS2(hf:hf+99,:,:)))';
responseSpkCount = responseSpkCount - repmat(fr,size(responseSpkCount,1),1);

nU = size(responseSpkCount,2);

nEvents = strcmp({D(:).direction},'NasalMimic');
tEvents = strcmp({D(:).direction},'TemporalMimic');

ginicoeff = zeros(nU,1);
gini_p = ones(nU,1);
for i=1:nU
    spk_i = responseSpkCount(:,i);
    spkNasal = spk_i(nEvents);
    spkTemporal = spk_i(tEvents);
    if all(sort(spkNasal) == sort(spkTemporal)) % If they are identical
        ginicoeff(i) = 0;
    else
        a = roc_curve(spkTemporal,spkNasal,0,0);
        ginicoeff(i) = 2*a.param.AROC-1;
        gini_p(i) = ranksum(spkTemporal,spkNasal);
    end
end

end