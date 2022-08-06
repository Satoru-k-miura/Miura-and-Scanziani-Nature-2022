function sacResp = saccadeDirectionResponse2D
%% Read in Data High format file
while 1
    [fname,fpath]=uigetfile('_DataHighFormat.mat','Choose existing data structure with raw spike series for analysis','MultiSelect','off');
    if fpath == 0
        error('Execution cancelled');
    else
        if contains(fname,'DataHighFormat')
            C = load(fullfile(fpath,fname));
            dat = C.dat;
            break
        end
    end
end

%% Calculate the saccade direction
theta=[];
numspk = [];
bl_numspk = [];
nr=[];
for i=1:length(dat)
    et = dat(i).eyeTrace;
    sacT = dat(2).time;
    idx = round((sacT(2)-sacT(1))/11);
    y = round(et(2,47+idx)-et(2,47));
    x = round(et(1,47+idx)-et(1,47));
    if isreal(y)&&isreal(x)
        t = atan2d(y,x);
        theta(i)=t;
        data = dat(i).data;
        bl_data = sum(data(:,211:311),2); %-300 to -200 ms
        data = sum(data(:,511:611),2);
        numspk(i,:) = data';
        bl_numspk(i,:) = bl_data';
    else
        theta(i)=nan;
    end
end
theta=theta+180; % nasal is 0, dorsal is 90
theta = theta';

theta(:,2)=0;
for i=1:7
    theta(theta(:,1)>(i-1)*45+22.5 & theta(:,1)<=i*45+22.5,2) = i*45;
end

theta(isnan(theta(:,1)),2)=nan;

Theta = num2cell(theta);
[dat.theta] = deal(Theta{:,1});
[dat.thetaGroup] = deal(Theta{:,2});

datAve=C.datAve;
params = C.params;

save(fullfile(fpath,fname),'dat','datAve','params','-mat'); %save new dat

for j=1:size(numspk,2)
    numEvents = zeros(1,8);
    meanSpk = zeros(1,8);
    stdSpk = zeros(1,8);
    spks = numspk(:,j);
    bl_spks = bl_numspk(:,j);
    kw=[];
    for i=0:7
        tt = i*45;
        sac = theta(:,2)==tt;
        numEvents(i+1) = sum(sac);
        s = spks(sac)';
        bl_s = bl_spks(sac)';
        meanSpk(i+1) = mean(s)*10; %Hz
        stdSpk(i+1) = std(s)*10;
        eval(['sacResp(j).spikeCount' num2str(tt) '= s;']);
        eval(['sacResp(j).baselinespikeCount' num2str(tt) '= bl_s;']);
        bl_s = bl_s';
        bl_s(:,2) = i;
        kw = [kw;bl_s];
    end
    sacResp(j).avSpikeCount = meanSpk;
    sacResp(j).stdSpikeCount = stdSpk;
    sacResp(j).order = 0:45:315;
    sacResp(j).binSize = 100; % ms
    
    p = kruskalwallis(kw(:,1),kw(:,2),'off');
    sacResp(j).pBase = p;
%     if p<0.05
%         sacResp(j).baseSig = true;
%     else
%         sacResp(j).baseSig = false;
%     end

    [~,idx] = max(meanSpk);
    m = (idx-1)*45;
    sacResp(j).prefDir = m;
    
    eval(['prefSpkCount = sacResp(j).spikeCount' num2str(m) ';']);
    if idx<=4
        np = m+180;
    else
        np = m+180-360;
    end
    eval(['nonprefSpkCount = sacResp(j).spikeCount' num2str(np) ';']);
    
    p = ranksum(prefSpkCount,nonprefSpkCount);
    sacResp(j).pSig = p;
%     if p<0.05
%         sacResp(j).sig = true;
%     else
%         sacResp(j).sig = false;
%     end
    a = roc_curve(nonprefSpkCount',prefSpkCount',0,0);
    sacResp(j).GiniCoeff = 2*a.param.AROC-1;
    
    meanSpkShift = circshift(meanSpk,3-idx);
    stdSpkShift = circshift(stdSpk,3-idx);
    
    sacResp(j).meanSpikeShift = meanSpkShift;
    sacResp(j).stdSpikeShift = stdSpkShift;
    
    figure;
    polarplot(deg2rad(0:45:360),[meanSpk meanSpk(1)])
    thetaticks(0:45:315)
    thetaticklabels({'Nasal','','Dorsal','','Temporal','','Ventral',''})
    title(['Unit ' num2str(j)])

    prefEvents = theta(:,2)==m;
    nonprefEvents = theta(:,2)==np;
    
    datPref = dat(prefEvents);
    datNonpref = dat(nonprefEvents);
    
    dataPref = [datPref.data];
    dataPref = reshape(dataPref(j,:),1020,[]);
    dataPref = dataPref'; % events by time
    dataNonpref = [datNonpref.data];
    dataNonpref = reshape(dataNonpref(j,:),1020,[]);
    dataNonpref = dataNonpref'; % events by time
    
    sacResp(j).preferredResponse = dataPref;
    sacResp(j).nonpreferredResponse = dataNonpref;
    
end

