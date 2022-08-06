function resampleShift3

while 1
    [fnameS,fpathS]=uigetfile('_DataHighFormat.mat','Choose DHF from shift experiment for resampling','MultiSelect','off');
    if fpathS == 0
        error('Execution cancelled');
    else
        if contains(fnameS,'DataHighFormat')
            C = load(fullfile(fpathS,fnameS));
            D = C.dat;
            break
        end
    end
end

% Check is nasal amplitude is negative, if not convert
DNasal = IndexedStructCopy(D, strcmp({D(:).direction},'NasalMimic'));
ampNasal = mean([DNasal(:).amplitude]);
if ampNasal>=0
    for i=1:length(D)
        if strcmp(D(i).direction,'NasalMimic')
            D(i).amplitude = -D(i).amplitude;
        end
    end
end
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
        if abs(delta)<2
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
dat=E;
params = C.params;

[~,name,ext] = fileparts(fnameS);

save(fullfile(fpathS,[name '_resamp3' ext]),'dat','params','delta_amplitude');
    
    
    
