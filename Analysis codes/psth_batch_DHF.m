function [PSTH,nonz]=psth_batch_DHF(unitCriteria)

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

VS = 'Blind';
eval(['psthNT.psth_' VS '= [];'])
for i=1:length(fname)
    fname_i = fname{i};
    R = generatepsthforeach(fpath,fname_i,VS);
    eval(['psthNT.psth_' VS '= [psthNT.psth_' VS ';R.psth_' VS '];'])
end

%Figure
[PSTH,nonz]=psth_VS(psthNT,unitCriteria,VS);

end

function R = generatepsthforeach(fpath,fname,VS)
%%
kernSD = 20;

C = load(fullfile(fpath,fname));
D = C.dat;

%% check for number of conditions

D1 = IndexedStructCopy(D, strcmp({D(:).VS},VS));

%% Calculate PSTH

spkCounts = zeros(size(D1(1).data));
for i=1:length(D1)
    spkCounts = spkCounts+D1(i).data;
end

PSTH=[];
st = 1;
for i=1:size(spkCounts,2)/kernSD
    ed = st+kernSD-1;
    PSTH = [PSTH sum(spkCounts(:,st:ed),2)];
    st = ed+1;
end
PSTH = PSTH/length(D1);
eval(['R.psth_' VS '= PSTH/kernSD*1e3;']); %spks/sec


end

function [allGrating,nonz]=psth_VS(psthNT,R,VS)

if isempty(R)
    eval(['R = true(size(psthNT.psth_' VS ',1),1);'])
end

eval(['allGrating = psthNT.psth_' VS '(R,:);'])

%% Standardize

bl_allGrating = mean(allGrating(:,1:(300/20)),2);

nonz = bl_allGrating>0.5;

bl_allGrating = bl_allGrating(nonz);
allGrating = allGrating(nonz,:);
allGrating = allGrating ./ bl_allGrating;

%%

numU = size(allGrating,1); % number of units
x = -500:20:500;
X = [x,fliplr(x)];

figure; % grating
set(gcf,'color','white');
hold on

avPSTH = nanmean(allGrating,1);
stdPSTH = nanstd(allGrating,0,1);
stdEPSTH = stdPSTH / sqrt(numU);
Y = [avPSTH-stdEPSTH,fliplr(avPSTH+stdEPSTH)];
fill(X,Y,[0.7 0.7 0.7],'EdgeColor','none');
plot(x,avPSTH,'k');

xlabel('Time from saccade onset (ms)');
ylabel('Average spike rate');
pbaspect([1 1 1])

end
