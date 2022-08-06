function DirdiscSU_DataHighFormat(varargin)
% parameters for DataHigh

% kernel SD used for GPFA
kernSD = 20;

before = -500;
after = 500;
condVS = {'condition1','condition2'}; % only if new file
sortBy = 'VS'; % colors to separate different events alt: 'direction'
fnameAppend = 'X';

assignopts(who, varargin);
% varargin:
%   before: time before saccade onset, default -500
%   after: time after saccade onset, default 500
%   condNames: condition names, provided as 1x2 cell array, first is
%   'blank' trials, second is 'non-blank' trials
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
        result = dataStructForGPFA_saccades2(kernSD,'before',before,'after',after,'condVS',condVS,'sortBy',sortBy); % condition is based on VS channel
%         result = dataStructForGPFA_saccades_trialAve(kernSD,'before',before,'after',after,'condNames',condNames); % condition is based on VS channel
        
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
numEvents1 = sum(strcmp({D(:).direction},'Nasal'));
numEvents2 = sum(strcmp({D(:).direction},'Temporal'));
Dir1 = zeros(numSUs,numBins,numEvents1);
Dir2 = zeros(numSUs,numBins,numEvents2);
Dir3 = zeros(numSUs,numBins,numEvents1+numEvents2);
cm1 = 0;
cm2 = 0;
cm3 = 0;
for i=1:length(D)
    if strcmp({D(i).direction},'Nasal')
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
% spkNasal = squeeze(sum(Dir1(:,midT+100:midT+100+winWidth,:),2));
baseNasal = squeeze(sum(Dir1(:,midT-201-winWidth:midT-201,:),2));

spkTemporal = squeeze(sum(Dir2(:,midT:midT+winWidth,:),2));
% spkTemporal = squeeze(sum(Dir2(:,midT+100:midT+100+winWidth,:),2));
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

baseNasal2 = squeeze(sum(Dir1(:,midT-201-200:midT-201,:),2));
baseTemporal2 = squeeze(sum(Dir2(:,midT-201-200:midT-201,:),2));
baseSig = false(numSUs,1);
pBase = ones(numSUs,1);
for i=1:numSUs
    p=ranksum(baseNasal2(i,:),baseTemporal2(i,:));
%     p=ranksum(spkNasal(i,:)-mean(baseNasal(i,:)),spkTemporal(i,:)-mean(baseTemporal(i,:)));
%     p=ranksum(spkNasal(i,:)-baseNasal(i,:),spkTemporal(i,:)-baseTemporal(i,:));
    pBase(i) = p;
    if p<0.05
        baseSig(i) = 1;        
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
% figure;
% x=histbins(1:end-1)+binsize/2;
% b=bar(x,[sigHist nonsigHist],1,'stacked');
% b(1).EdgeColor = 'none';
% b(2).EdgeColor = 'none';
% 
% b(1).FaceColor = [0 0.4470 0.7410];
% b(2).FaceColor = [0.3010 0.7450 0.9330];
% xlabel('Gini coeff')
% ylabel('Frequency of occurence')
% 
% figure;
% scatter(foldChange(~sig,1),foldChange(~sig,2)); hold on
% scatter(foldChange(sig,1),foldChange(sig,2),'filled');
% xlabel('Fold change over baseline, nasal saccades')
% ylabel('Fold change over baseline, temporal saccades')

includedUnits = tU;
aroc = df;
save(['Dirdisc_SU_aroc_',num2str(winWidth),'_',fnameAppend],'numSpkNasal','stdSpkNasal','numSpkNasalChange','numSpkTemporal','stdSpkTemporal','numSpkTemporalChange','aroc','foldChange','pDir','sig','changeSig','includedUnits','baseSig','pBase')
%% Plot SUs
% Parameters for the raster plots


% LineFormat = struct();
% LineFormat.Color = [0 0 0];
% LineFormat.Linewidth = 1.0;
% LineFormat.LineStyle = '-';
% wdw = 15; % window size for psth
% pd = floor(mod(numBins,wdw)/2);
% numNewBins = floor(numBins/wdw);
% for i=1:numSUs
%     figure;
%     vs = logical(squeeze(Dir1(i,:,:))');
%     subplot(4,1,1);
%     plotSpikeRaster(vs,'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',0.8,'AutoLabel',true,'XForLogical',[before-kernSD/2 after+kernSD/2]);
%     set(gca,'XLim',[before after],'YLim',[0 size(vs,1)],'box','off','Visible','off');
%     subplot(4,1,2);
%     vs2 = convLogSpk(vs,15);
%     vs2 = mean(vs2,1);
%     vsN = zeros(1,numNewBins);
%     for ii=1:numNewBins
%         strt = (ii-1)*wdw+pd+1;
%         vsN(ii) = sum(vs2(strt:strt+wdw-1))/wdw*1e3;
%     end
%     x = (1+pd+floor(wdw/2):wdw:1+pd+floor(wdw/2)+wdw*(numNewBins-1))-midT;
%     plot(x,mean(vsN,1),'k');
%     xlim([before after]);
%     set(gca,'XTick',before:100:after,'TickDir','Out','box','off');
%     
%     vs = logical(squeeze(Dir2(i,:,:))');
%     subplot(4,1,3);
%     plotSpikeRaster(vs,'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',0.8,'AutoLabel',true,'XForLogical',[before-kernSD/2 after+kernSD/2]);
%     set(gca,'XLim',[before after],'YLim',[0 size(vs,1)],'box','off','Visible','off');
%     subplot(4,1,4);
%     vs = convLogSpk(vs,15);
%     vs = mean(vs,1);
%     vsN = zeros(1,numNewBins);
%     for ii=1:numNewBins
%         strt = (ii-1)*wdw+pd+1;
%         vsN(ii) = sum(vs(strt:strt+wdw-1))/wdw*1e3;
%     end
%     x = (1+pd+floor(wdw/2):wdw:1+pd+floor(wdw/2)+wdw*(numNewBins-1))-midT;
%     plot(x,mean(vsN,1),'k');
%     xlim([before after]);
%     set(gca,'XTick',before:100:after,'TickDir','Out','box','off');
%     
%     suptitle(['SU ' num2str(i)])
%     
%     figure;
%     vs = logical(squeeze(Dir3(i,:,:))');
%     subplot(4,1,1);
%     plotSpikeRaster(vs,'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',0.8,'AutoLabel',true,'XForLogical',[before-kernSD/2 after+kernSD/2]);
%     set(gca,'XLim',[before after],'YLim',[0 size(vs,1)],'box','off','Visible','off');
%     subplot(4,1,2);
%     vs2 = convLogSpk(vs,15);
%     vs2 = mean(vs2,1);
%     vsN = zeros(1,numNewBins);
%     for ii=1:numNewBins
%         strt = (ii-1)*wdw+pd+1;
%         vsN(ii) = sum(vs2(strt:strt+wdw-1))/wdw*1e3;
%     end
%     x = (1+pd+floor(wdw/2):wdw:1+pd+floor(wdw/2)+wdw*(numNewBins-1))-midT;
%     plot(x,mean(vsN,1),'k');
%     xlim([before after]);
%     set(gca,'XTick',before:100:after,'TickDir','Out','box','off');
% end
end

function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S)';
end 
FieldList{2,1} = {};
T = struct(FieldList{:});
numDat = sum(Condition);
T(1:numDat) = S(Condition);
end

function s = convLogSpk(logSpk,sigma)
resolution = 1;
%Time ranges form -3*st. dev. to 3*st. dev.
edges = (-3*sigma:resolution:3*sigma);
%Evaluate the Gaussian kernel
kernel = normpdf(edges,0,sigma);
%Multiply by bin width so the probabilities sum to 1
kernel = kernel.*resolution; 
s=[];
% for each event, convolve spike data with the kernel
for i=1:size(logSpk,1)
    c = conv(logSpk(i,:),kernel);
    s = [s;c];
end
%Find the index of the kernel center
center = ceil(length(edges)/2); 
%Trim out the relevant portion of the spike density
s = s(:,center:end-center+1); 
end