function resampleShift(varargin)

twice = false;
matchingVS = 'statVertical';

assignopts(who,varargin);

currentFolder = pwd;

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
ampShift = [D(:).amplitude];
    

while 1
    [fname,fpath]=uigetfile('_DataHighFormat.mat','Choose DHF to which saccades will be matched','MultiSelect','off');
    if fpath == 0
        error('Execution cancelled');
    else
        if contains(fname,'DataHighFormat')
            F = load(fullfile(fpath,fname));
            E = F.dat;
            break
        end
    end
end

E = IndexedStructCopy(E, strcmp({E(:).VS},matchingVS));

amp = [E(:).amplitude];

delta_amplitude = []; % for reporting the difference between real and fake saccades
idx=[];

dat = D(1); %just for initialization. delete later.

for i=1:length(amp)
    ai = amp(i);
    [~,I] = min(abs(ampShift - ai));
    dat = [dat D(I)];
    delta = ampShift(I)-ai;
    delta_amplitude(end+1) = delta;
    
    %update D to exclude I
    tr = ~((1:length(D))==I);
    D = IndexedStructCopy(D,tr);
    ampShift = [D(:).amplitude];
end

if twice
    for i=1:length(amp)
        ai = amp(i);
        [~,I] = min(abs(ampShift - ai));
        dat = [dat D(I)];
        delta = ampShift(I)-ai;
        delta_amplitude(end+1) = delta;

        %update D to exclude I
        tr = ~((1:length(D))==I);
        D = IndexedStructCopy(D,tr);
        ampShift = [D(:).amplitude];
    end
    
    dat = dat(2:end);

    params = C.params;

    [~,name,ext] = fileparts(fnameS);

    save(fullfile(currentFolder,[name '_resampTWICE' ext]),'dat','params','delta_amplitude');
    
else
    dat = dat(2:end);

    params = C.params;

    [~,name,ext] = fileparts(fnameS);

    save(fullfile(currentFolder,[name '_resamp' ext]),'dat','params','delta_amplitude');

end
        