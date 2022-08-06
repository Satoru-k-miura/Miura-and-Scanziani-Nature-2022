function saccadeClassifier(varargin)

%% Params
% params for PCA to LDA
thrFR = 1;
cnvSD = 20; % Convolve with Gaussian, SD in ms, cnvSD=0 is no smoothing
spikeShape = 'both';
numLDALoops = 1;
numSUperm = 10;
nSU = 10;
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
params.numSUperm = numSUperm;
params.nSU = nSU;
params.numK = numK;
params.shuffle = shuffle;
params.totalVariance = totalVariance; % min total variance accounted for in PCA
params.modelClass = modelClass;
%% pool data
MimicD1 = poolMimicData_v2;
RealD = poolRealData;

resultsMimicSummary = [];
resultsRealSummary = [];
poolobj = parpool('local',4);
for i=1:numRandSample
    disp(['Iteration ' num2str(i) ' of ' num2str(numRandSample)])
    MimicD = resampleLoc(MimicD1);
    MimicD.dat = MimicD.dat(1:floor(floor(length(MimicD.dat)/numK)/2)*2*numK);
    %% Perform PCA to QDA
    params.unitCondition = [];
    [resultsMimic,resultsReal,params] = PCAtoXDALoc(MimicD,RealD,params,i);
    resultsMimicSummary = [resultsMimicSummary;resultsMimic];
    resultsRealSummary = [resultsRealSummary;resultsReal];

end
delete(poolobj);

%% Save summary
if strcmp(params.modelClass,'LDA')
    fname = 'simpleClassifierLDASummary';
elseif strcmp(params.modelClass,'QDA')
    fname = 'simpleClassifierQDASummary';
else
    fname = 'simpleClassifierNNSummary';
end
save(fname,'resultsMimicSummary','resultsRealSummary','nSU','params');


%% Plot
len = length(nSU)+1;
dk = [0 107 179]/255;
lt = [165 222 255]/255;
colors_p = [linspace(lt(1),dk(1),len)', linspace(lt(2),dk(2),len)', linspace(lt(3),dk(3),len)'];

figure;hold on
plot(0,0.5,'o','Color',colors_p(1,:));
pks = [];
for i=1:length(nSU)
    acc = cell2mat(resultsMimicSummary(:,i));
    accMean = mean(acc);
    pks = [pks accMean];
    plot(nSU(i),accMean,'o','Color',colors_p(i+1,:));
end
xlim auto
x0 = 0;
y0 = 0.5;
% fit exp
g = @(A,p,x)(A-(A-y0)*exp(-p*(x-x0)));
f=fit([0 nSU]',[0.5 pks]',g);
plot(f)
hold off
title('K-fold classification, fake saccades')

figure;hold on
plot(0,0.5,'o','Color',colors_p(1,:));
pks = [];
for i=1:length(nSU)
    acc = cell2mat(resultsRealSummary(:,i));
    accMean = mean(acc);
    pks = [pks accMean]; 
    plot(nSU(i),accMean,'o','Color',colors_p(i+1,:));
end
xlim auto
x0 = 0;
y0 = 0.5;
% fit exp
g = @(A,p,x)(A-(A-y0)*exp(-p*(x-x0)));
f=fit([0 nSU]',[0.5 pks]',g);
plot(f)
hold off
title('Classification of real saccades')
end

function [resultsMimic,resultsReal,params] = PCAtoXDALoc(MimicD,RealD,params,randNum)

thrFR = params.thrFR;
cnvSD = params.cnvSD;
unitCondition = params.unitCondition;
spikeShape = params.spikeShape;
numLDALoops = params.numLDALoops;
numSUperm = params.numSUperm;
nSU = params.nSU;
shuffle = params.shuffle;
totalVariance = params.totalVariance;

if size(unitCondition,1)==1
    unitCondition = unitCondition';
end

params.kernSD = MimicD.params.kernSD;
params.before = MimicD.params.before;
params.after = MimicD.params.after;

%%
D = MimicD.dat;
kernSD = MimicD.params.kernSD;
before = MimicD.params.before;
after = MimicD.params.after;

%%
resultsMimic = cell(numSUperm,length(nSU));
varExplained = zeros(numSUperm,length(nSU));
finalDims = zeros(numSUperm,length(nSU));


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

if ~isempty(unitCondition)
    unitsKeep = unitsKeep & unitCondition;
end

allSpk=cat(1,D(:).data)';
aS = reshape(allSpk,1020,[],length(D));
aS2 = squeeze(sum(aS(21:220,:,:)));
fr = mean(aS2,2)/2;

aS3 = squeeze(sum(aS(511:610,:,:))); % 100ms post onset

% unitsKeep = fr*10>=thrFR;
% 
% if ~isempty(unitCondition)
%     unitsKeep = unitsKeep & unitCondition;
% end

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

R = RealD.dat;
if shuffle % If shuffle, shuffle condition
    so = randsample(1:length(R),length(R),'false'); % shuffling order
    [R.direction] = R(so).direction;
end
allSpk=cat(1,R(:).data)';
aS = reshape(allSpk,1020,[],length(R));
aS2 = squeeze(sum(aS(21:220,:,:)));
fr = mean(aS2,2)/2;
aS3 = squeeze(sum(aS(511:610,:,:))); % 100ms post onset
for i=1:length(R) % select the same units from real saccade data
%     t = aS3(unitsKeep,i);
%     t = t - fr(unitsKeep);
%     R(i).data = t;

    t = R(i).data(unitsKeep,:);
    t = t';
    t = reshape(t,20,[]);
    t = sum(t,1);
    t = reshape(t,[],nU)';
    if cnvSD > 0
        t = smoother(t,cnvSD,20);
    end
    R(i).data = t(:,29);

%     t = R(i).data(unitsKeep,:);
%     if cnvSD > 0
%         t = smoother(t,cnvSD,1);
%     end
%     R(i).data = mean(t(:,511:711),2);
end

resultsReal = cell(numSUperm,length(nSU));

for l = 1:length(nSU)
    nS = nSU(l);
    
    parfor i=1:numSUperm
        
        % Determine singularity
        sing = true;
        while sing

            tu = randsample(nU,nS);
            data=cat(1,D(:).data);
            d = data; % concatenated data for the condition
            d = reshape(d,nU,[]); % reshape such that nU x events
            d = d(tu,:)';

            %% Perform XDA on PCA
            meand = mean(d);
            stdd = std(d);

            d = (d-meand)./stdd;
            
            d(isnan(d)) = 0;

            [coeff,score,~,~,explained,~] = pca(d);
            w = warning('query','last');
            warning('off',w.identifier);

            expTot = 0;
            dims = 0;
            while sum(expTot)<totalVariance && dims < length(D)/5 %dims <= nS-1 && dims < length(D)/5 
                dims = dims+1;
                expTot = sum(explained(1:dims));
            end
%             varExplained(i,l) = expTot;
%             finalDims(i,l) = dims;

            % reduce
            dataPCA = score(:,1:dims);
            nrows = strcmp({D(:).direction},'NasalMimic');
            trows = strcmp({D(:).direction},'TemporalMimic');
            
            A = cov(dataPCA(nrows,:));
            B = cov(dataPCA(trows,:));
            if rcond(A)>10e-10 && rcond(B)>10e-10
                sing = false;
            end
        end
        if strcmp(params.modelClass,'LDA')
            Mdl = fitcdiscr(dataPCA,{D(:).direction},'SaveMemory','on','FillCoeffs','on');
        elseif strcmp(params.modelClass,'QDA')
            Mdl = fitcdiscr(dataPCA,{D(:).direction},'DiscrimType','quadratic','SaveMemory','on','FillCoeffs','on'); 
        elseif strcmp(params.modelClass,'NN')
            Mdl = fitcnet(dataPCA,{D(:).direction}); 
        else
            error('Specify model name')
        end
    
        dataReal = cat(1,R(:).data);

        dR = dataReal;
        dR = reshape(dR,nU,[]); % reshape such that nU x events
        dR = dR(tu,:)'; % events x nU

        centerdR = dR-meand;
        centerdR = centerdR./stdd;

        score2 = centerdR*coeff;
        dataPCAReal = score2(:,1:dims);

        label = predict(Mdl,dataPCAReal);
        for ii=1:length(label)
            if strcmp(label{ii},'NasalMimic')
                label{ii} = 'Nasal';
            elseif strcmp(label{ii},'TemporalMimic')
                label{ii} = 'Temporal';
            end
        end
        CM = confusionmat({R(:).direction},label);
        acc = (CM(1,1)+CM(2,2))/length(label);
        resultsReal{i,l} = acc;

        %% Cross validate
        Dcv = D;
        for ii = 1:length(Dcv)
            t = Dcv(ii).data(tu,:);
            Dcv(ii).data = t;
        end 
        acc = XDAonPCAcrossval(Dcv,params);
        resultsMimic{i,l} = acc;
    end
end
if strcmp(params.modelClass,'LDA')
    fname = ['saccadeClassifier_LDA_rep' num2str(randNum)];
elseif strcmp(params.modelClass,'QDA')
    fname = ['saccadeClassifier_QDA_rep' num2str(randNum)];
elseif strcmp(params.modelClass,'NN')
    fname = ['saccadeClassifier_NN_rep' num2str(randNum)];
end
save(fname, 'kernSD','before','after','resultsMimic','resultsReal','finalDims','varExplained','params');

end

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

%% Cross validation
function acc = XDAonPCAcrossval(D,params)

numK = params.numK;
totalVariance = params.totalVariance;

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
    while sum(expTot)<totalVariance && dims < length(trainevents)/5 %dims <= nU-1 && dims < length(trainevents)/5 %sum(expTot)<totalVariance && dims < length(trainevents)/5
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
    
    nrows = strcmp({D(trainevents).direction},'NasalMimic');
    trows = strcmp({D(trainevents).direction},'TemporalMimic');

    A = cov(dataPCA_train(nrows,:));
    B = cov(dataPCA_train(trows,:));
    if rcond(A)<10e-10 || rcond(B)<10e-10
        allacc(k) = nan;
    else
        if strcmp(params.modelClass,'LDA')
            Mdl = fitcdiscr(dataPCA_train,{D(trainevents).direction},'SaveMemory','on','FillCoeffs','on'); %DA
        elseif strcmp(params.modelClass,'QDA')
            Mdl = fitcdiscr(dataPCA_train,{D(trainevents).direction},'DiscrimType','quadratic','SaveMemory','on','FillCoeffs','on'); %DA
        elseif strcmp(params.modelClass,'NN')
            Mdl = fitcnet(dataPCA_train,{D(trainevents).direction});
        end
        label_test = predict(Mdl,dataPCA_test);

        CM = confusionmat({D(testevents).direction},label_test);
        acc = (CM(1,1)+CM(2,2))/length(label_test);
        allacc(k) = acc;
    end
end
    
acc = nanmean(allacc);
end