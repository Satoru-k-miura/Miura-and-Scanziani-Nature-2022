function fastvsreg = plotclusWFs

while 1
    [tetname,tetpath]=uigetfile('times_polytrodeAll.mat','Choose times_polytrodeAll.mat files','MultiSelect','off');
    if tetpath == 0
        error('Execution cancelled');
    end

    nameMatch = strfind(tetname,'times_polytrodeAll');

    if ~isempty(nameMatch)
        break;
    end
end

polyAll = load(fullfile(tetpath,tetname),'-mat');

normShapes = [];
for i=1:size(polyAll.clusWFs,1)
    clusWF = squeeze(polyAll.clusWFs(i,:,:))';
    ampWF = squeeze(squeeze(abs(min(clusWF,[],2)-mean(clusWF(:,1:8),2))));
    [~,maxChan] = max(ampWF);
    medShape = clusWF(maxChan,:)-mean(clusWF(maxChan,1:8));
    normShapes(i,:) = medShape/abs(min(medShape));
%     figure;plot(medShape)
end

x=((1:64)-20)/30;

% find min
[~,Ii] = min(normShapes(:,1:20),[],2);
% find max between the trough and the end
[~,I] = max(normShapes(:,20:end),[],2);
t2p = (I+19-Ii)/30;

w = [];
thr=0.5; % relative to min trough
for i=1:size(normShapes,1) % for every spike
    idx = find((normShapes(i,1:20)-(-thr))<=0,1);
    ini = abs(thr-normShapes(i,idx-1))/abs(normShapes(i,idx-1)-normShapes(i,idx))/30 + x(idx-1);
    idx = find((normShapes(i,21:end)-(-thr))>=0,1);
    sec = abs(thr-normShapes(i,19+idx))/abs(normShapes(i,19+idx)-normShapes(i,20+idx))/30 + x(19+idx);
    if isempty(sec) || isempty(ini)
        w(end+1) = NaN;
    else
        w(end+1) = abs(sec-ini);
    end
end

figure;scatter(w,t2p);
fastvsreg=cell(length(t2p),1);
for i=1:length(t2p)
    if t2p(i)<0.5
        fastvsreg{i} = 'fast';
    else
        fastvsreg{i} = 'regular';
    end
end
% fastvsreg(t2p<0.5) = 1; % 1 for fast
% fastvsreg(t2p>=0.5) = 2; % 2 for regular
save(fullfile(tetpath,tetname),'fastvsreg','-append');