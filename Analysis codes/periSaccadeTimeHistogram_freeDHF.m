function periSaccadeTimeHistogram_freeDHF(varargin)

% create average PSTH for each neuron around saccades

% kernel SD used for GPFA
kernSD = 20;
before = -500;
after = 500;
condVS = {'condition1','condition2'}; % only relevant when making new file
sortBy = 'VS'; % colors to separate different events alt: 'direction'
thrFR = 1;
cnvSD = 15; % Convolve with Gaussian, SD in ms, cnvSD=0 is no smoothing
unitCondition = [];
spikeShape = 'both';
numLoops = 1;
numK = 5;
shuffle = 0;
SUs = [];
assignopts(who, varargin);
if size(unitCondition,1)==1
    unitCondition = unitCondition';
end

%% Read in Data High format file
while 1
    [fname,fpath]=uigetfile('_DataHighFormat.mat','Choose existing data structure with raw spike series for analysis','MultiSelect','off');
    if fpath == 0
        error('Execution cancelled');
    else
        if contains(fname,'DataHighFormat')
            C = load(fullfile(fpath,fname));
            D = C.dat;
            break
        end
    end
end

%% Select relevant conditions

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

%% Parameters
% Bin size for PSTH in ms
params.binsize_hist = 20; 
% edges for PSTH
params.edges = -2000:params.binsize_hist:2000+params.binsize_hist; 
% params.edges = -1500:params.binsize_hist:1500+params.binsize_hist;
% Parameters for the raster plots
LineFormat = struct();
LineFormat.Color = [0 0 0];
LineFormat.Linewidth = 1.0;
LineFormat.LineStyle = '-';

%% pool all data for all neurons
spikesN = [];
spikesT = [];
for i=1:length(D)
    data=D(i).data;
    
    % Create D.direction (doesn't exist in DHF from freely moving exp)
    ang = D(i).theta;
%     if (ang>=0 && ang<45) || (ang>=315)
%         D(i).direction = 'Nasal';
%         spikesN = cat(3,spikesN,data);
%     elseif ang>=135 && ang<225
%         D(i).direction = 'Temporal';
%         spikesT = cat(3,spikesT,data);
%     else
%         D(i).direction = 'NA';
%     end
    if (ang>=0 && ang<22.5) || (ang>=337.5)
        D(i).direction = 'Nasal';
        spikesN = cat(3,spikesN,data);
    elseif ang>=157.5 && ang<205.2
        D(i).direction = 'Temporal';
        spikesT = cat(3,spikesT,data);
    else
        D(i).direction = 'NA';
    end
end

if isempty(SUs)
    r = 1:size(spikesN,1);
else
    r = SUs;
end

for i=r
    C = squeeze(spikesN(i,:,:))';
    D = squeeze(spikesT(i,:,:))';
    
    figure;
    set(gcf,'color','white');
    subplot(3,1,1);
    plotSpikeRaster(logical(C),'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',0.8,'AutoLabel',true,'rasterWindowOffset',-0.510);
    line([0 0],[0 size(C,1)+1],'LineWidth',1,'Color',[1 0 0]);
    set(gca,'XLim',[-500 500],'YLim',[0 size(C,1)+1],'box','off','Visible','off','XTick',-500:500:500); 
    subplot(3,1,2);
    plotSpikeRaster(logical(D),'PlotType','vertline','LineFormat',LineFormat,'VertSpikeHeight',0.8,'AutoLabel',true,'rasterWindowOffset',-0.510);
    line([0 0],[0 size(D,1)+1],'LineWidth',1,'Color',[1 0 0]);
    set(gca,'XLim',[-500 500],'YLim',[0 size(D,1)+1],'box','off','Visible','off','XTick',-500:500:500); 
    
    N = sum(C,1);
    psthN = [];
    st = 1;
    for j=1:(length(N)/params.binsize_hist)
        ed = st+params.binsize_hist-1;
        psthN = [psthN sum(N(st:ed))];
        st = ed+1;
    end
    psthN = psthN/size(C,1); % per trial
    psthN = psthN/params.binsize_hist*1e3; %spks/sec
    
    T = sum(D,1);
    psthT = [];
    st = 1;
    for j=1:(length(T)/params.binsize_hist)
        ed = st+params.binsize_hist-1;
        psthT = [psthT sum(T(st:ed))];
        st = ed+1;
    end
    psthT = psthT/size(D,1); % per trial
    psthT = psthT/params.binsize_hist*1e3; %spks/sec
    
    subplot(3,1,3);
    x = before:params.binsize_hist:after;
    plot(x,psthN,'b','LineWidth',1.5);
    hold on
    plot(x,psthT,'r','LineWidth',1.5);
    line([0 0],[0 200],'LineWidth',1,'Color',[1 0 0]);
    set(gca,'XLim',[-500 500],'YLim',[0 30],'XTick',-500:500:500,'TickDir','Out','box','off');
    ylabel('Firing rate (spikes/s)');xlabel('Time from saccade onset (ms)');
    
    tet = fname(1:strfind(fname,'_DataHighFormat')-1);
    str = strcat(tet,' Cluster',num2str(i));
    h=suptitle(str);
    set(h,'interpreter','none');
    
end
