function [pref,nonpref,nonz]=psthPrefNonpref_batch_freeDHF(aroc,unitCriteria)

%% Read in Data High format file

while 1
    [fname,fpath]=uigetfile('_DataHighFormat.mat','Choose existing data structure with raw spike series for analysis','MultiSelect','on');
    if fpath == 0
        error('Execution cancelled');
    else
        if all(contains(fname,'DataHighFormat'))
            break
        end
    end
end
fname = sort(fname);

VS = 'free';
eval(['psthNT.psth_N_' VS '= [];'])
eval(['psthNT.psth_T_' VS '= [];'])
for i=1:length(fname)
    fname_i = fname{i};
    R = generatepsthforeach(fpath,fname_i,VS);
    eval(['psthNT.psth_N_' VS '= [psthNT.psth_N_' VS ';R.psth_N_' VS '];'])
    eval(['psthNT.psth_T_' VS '= [psthNT.psth_T_' VS ';R.psth_T_' VS '];'])
end

%Figure
[pref,nonpref,nonz]=psth_prefNonpref_VS(psthNT,aroc,unitCriteria,VS);

% significance
for i=1:size(pref,2)
    p(i) = signrank(pref(:,i),nonpref(:,i),'tail','right');
end

sig = benjamini_hochberg(p,0.1);
p(~sig) = inf;

figure;imagesc(p);
pbaspect([1 1 1])
end

function R = generatepsthforeach(fpath,fname,VS)
%%
kernSD = 20;

C = load(fullfile(fpath,fname));
D = C.dat;

%% check for number of conditions

D1 = IndexedStructCopy(D, strcmp({D(:).VS},VS));

for i=1:length(D1)
    data=D1(i).data;
    
    % Create D.direction (doesn't exist in DHF from freely moving exp)
    ang = D1(i).theta;
    if (ang>=0 && ang<45) || (ang>=315)
        D1(i).direction = 'Nasal';
    elseif ang>=135 && ang<225
        D1(i).direction = 'Temporal';
    else
        D1(i).direction = 'NA';
    end
end
%% Calculate PSTH
resp = strcmp({D1(:).direction},'Nasal');
Dn = D1(resp);
Dt = D1(~resp);

spkCounts = zeros(size(Dn(1).data));
for i=1:length(Dn)
    spkCounts = spkCounts+Dn(i).data;
end
PSTH_N=[];
st = 1;
for i=1:size(spkCounts,2)/kernSD
    ed = st+kernSD-1;
    PSTH_N = [PSTH_N sum(spkCounts(:,st:ed),2)];
    st = ed+1;
end
PSTH_N = PSTH_N/length(Dn);
eval(['R.psth_N_' VS '= PSTH_N/kernSD*1e3;']); %spks/sec

spkCounts = zeros(size(Dt(1).data));
for i=1:length(Dt)
    spkCounts = spkCounts+Dt(i).data;
end

PSTH_T=[];
st = 1;
for i=1:size(spkCounts,2)/kernSD
    ed = st+kernSD-1;
    PSTH_T = [PSTH_T sum(spkCounts(:,st:ed),2)];
    st = ed+1;
end
PSTH_T = PSTH_T/length(Dt);
eval(['R.psth_T_' VS '= PSTH_T/kernSD*1e3;']); %spks/sec

end

function [prefGrating,nonprefGrating,nonz]=psth_prefNonpref_VS(psthNT,aroc,R,VS)

if isempty(R)
    R = true(size(aroc));
end

eval(['nGrating = psthNT.psth_N_' VS '(R,:);'])
eval(['tGrating = psthNT.psth_T_' VS '(R,:);'])

aroc = aroc(R);

nt = aroc>0; % false temporal, true nasal

%% Standardize

% bl_nGrating = mean(nGrating(:,1:(300/20)),2);
% bl_tGrating = mean(tGrating(:,1:(300/20)),2);
% 
% nonz = bl_nGrating>0.5 & bl_tGrating>0.5;
% 
% bl_nGrating = bl_nGrating(nonz);
% nGrating = nGrating(nonz,:);
% nGrating = nGrating ./ bl_nGrating;
% 
% % % 
% bl_tGrating = bl_tGrating(nonz);
% tGrating = tGrating(nonz,:);
% tGrating = tGrating ./ bl_tGrating;
% 
% nt = nt(nonz);
% 
% assignin('base','nonz',nonz)

%%

prefGrating = zeros(size(nGrating));
nonprefGrating = zeros(size(nGrating));

for i=1:length(nt)
    if nt(i) % if nasal is the preferred direction
        prefGrating(i,:) = nGrating(i,:);
        nonprefGrating(i,:) = tGrating(i,:);
    else
        prefGrating(i,:) = tGrating(i,:);
        nonprefGrating(i,:) = nGrating(i,:);
    end
end

numU = length(nt); % number of units
x = -500:20:500;
X = [x,fliplr(x)];

% bootstatPref = bootstrp(100,@(x)[mean(x)],prefGrating);
% bootstatNonpref = bootstrp(100,@(x)[mean(x)],nonprefGrating);
% 
% prefGrating = bootstatPref;
% nonprefGrating = bootstatNonpref;

figure; % grating
set(gcf,'color','white');
hold on

avPSTH = nanmean(prefGrating,1);
stdPSTH = nanstd(prefGrating,0,1);
stdEPSTH = stdPSTH / sqrt(numU);
Y = [avPSTH-stdEPSTH,fliplr(avPSTH+stdEPSTH)];
fill(X,Y,[0.9 0.6 0.6],'EdgeColor','none');
plot(x,avPSTH,'r');

avPSTH = nanmean(nonprefGrating,1);
stdPSTH = nanstd(nonprefGrating,0,1);
stdEPSTH = stdPSTH / sqrt(numU);
Y = [avPSTH-stdEPSTH,fliplr(avPSTH+stdEPSTH)];
fill(X,Y,[0.6 0.6 0.9],'EdgeColor','none');
plot(x,avPSTH,'b');

xlabel('Time from saccade onset (ms)');
ylabel('Average spike rate');
title('Pref vs Non-pref on gratings')
pbaspect([1 1 1])
set(gca,'TickDir','out')

%%
% figure; % grating
% set(gcf,'color','white');
% hold on
% 
% avPSTH = nanmean(tGrating,1);
% stdPSTH = nanstd(tGrating,0,1);
% stdEPSTH = stdPSTH / sqrt(numU);
% Y = [avPSTH-stdEPSTH,fliplr(avPSTH+stdEPSTH)];
% fill(X,Y,[0.9 0.6 0.6],'EdgeColor','none');
% plot(x,avPSTH,'r');
% 
% avPSTH = nanmean(nGrating,1);
% stdPSTH = nanstd(nGrating,0,1);
% stdEPSTH = stdPSTH / sqrt(numU);
% Y = [avPSTH-stdEPSTH,fliplr(avPSTH+stdEPSTH)];
% fill(X,Y,[0.6 0.6 0.9],'EdgeColor','none');
% plot(x,avPSTH,'b');
% 
% xlabel('Time from saccade onset (ms)');
% ylabel('Average spike rate');
% title('Nasal vs Temporal on gratings')
% pbaspect([1 1 1])

end
