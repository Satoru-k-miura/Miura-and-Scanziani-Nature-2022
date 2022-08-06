function eyeinfo = freeEyeTracking2(exp_id)
%%
% 3 files needed:
%       exp_id_DLCeyetrack.csv - eye tracking result from DLC
%       exp_id_DIN04 - eye tracking trigger on Intan
%       exp_id_timestamps.csv - Raspberry pi timestamps

% includes correction for CR tracking

%%
% params
sample_duration = 1/30000*1000; % Intan sample duration in ms
conv_x = 0.0192; % mm/pxl
conv_y = 0.0192; % mm/pxl
slope_Rp = -0.15036; % -0.1687; % slope of Rp to pupil diameter relationship 21-3: -0.15036 21-5: -0.05211
int_Rp = 1.2423; %1.2087; % intercept 21-3: 1.2423 21-5: 1.09702

% Read DeepLabCut CSV
dlc_fname = [exp_id '_DLCeyetrack.csv'];
delimiter = ',';
startRow = 4;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(dlc_fname,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

dataArray = dataArray([2 3 5 6 8 9 11 12 14 15 17 18 20 21 23 24 26 27]);
freeEye = cell2mat(dataArray);

clearvars delimiter startRow formatSpec fileID dataArray;

% Extract x and y coordinates for the pupil
freeEye_x = freeEye(:,[1 3 5 7 9 11 13 15]);
freeEye_y = freeEye(:,[2 4 6 8 10 12 14 16]);
CR = freeEye(:,[17 18]);

clear freeEye

% x1bound = [];
% x2bound = [];
% x3bound = [];
% x4bound = [];
% x5bound = [];
% x6bound = [];
% x7bound = [];
% x8bound = [];
% xbound = [-100 100];
% 
% y1bound = [];
% y2bound = [];
% y3bound = [];
% y4bound = [];
% y5bound = [];
% y6bound = [];
% y7bound = [];
% y8bound = [];
% ybound = [-100 100];

figure;
for i=1:8
    subplot(8,1,i)
    histogram(freeEye_x(:,i),'binWidth',1)
    ylim([0 50])
end
suptitle('X coordinates pupil')

figure;
for i=1:8
    subplot(8,1,i)
    histogram(freeEye_y(:,i),'binWidth',1)
    ylim([0 50])
end
suptitle('Y coordinates pupil')

figure;
histogram(CR(:,1),'binWidth',1)
ylim([0 25])
title('X coordinates CR')

figure;
histogram(CR(:,2),'binWidth',1)
ylim([0 25])
title('Y coordinates CR')

prompt = {'Enter x-coordinate bounds for pupil:','Enter y-coordinate bounds for pupil:',...
    'Enter x-coordinate bounds for CR:','Enter y-coordinate bounds for CR:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'[50 225]','[25 200]','[135 168]','[25 200]'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

roi_x = str2num(answer{1}); %[50 225]; % region of interest within the image
roi_y = str2num(answer{2}); %[25 200];
roi_CR_x = str2num(answer{3}); %[135 168];
roi_CR_y = str2num(answer{4}); %[25 200];


% fit ellipse to each frame
flat = [];
d_pxl = [];
x_pxl = [];
y_pxl = [];
phi = []; % tilt angle

for i=1:size(freeEye_x,1)
    dlc_x = freeEye_x(i,:)';
    dlc_y = freeEye_y(i,:)'; 
    inlim_x = (dlc_x<roi_x(2)&dlc_x>roi_x(1));
    inlim_y = (dlc_y<roi_y(2)&dlc_y>roi_y(1));
    inlim = inlim_x&inlim_y;
    if sum(inlim)>=5
        dlc_x = dlc_x(inlim);
        dlc_y = dlc_y(inlim);
        el = fit_ellipse(dlc_x,dlc_y,'n');
        if isempty(el.status)
            if el.X0_in>roi_x(1) && el.X0_in<roi_x(2) && el.Y0_in>roi_y(1) && el.Y0_in<roi_y(2) 
                f = (el.long_axis/2 - el.short_axis/2) / (el.long_axis/2);
                flat = [flat f];
                d_pxl = [d_pxl el.long_axis];
                x_pxl = [x_pxl el.X0_in];
                y_pxl = [y_pxl el.Y0_in];
                if (el.b>el.a)&&el.phi<0
                    phi = [phi 90+rad2deg(el.phi)];
                elseif (el.b>el.a)&&el.phi>0
                    phi = [phi 90-rad2deg(el.phi)];
                elseif (el.a>el.b)&&el.phi<0
                    phi = [phi abs(rad2deg(el.phi))];
                else
                    phi = [phi rad2deg(el.phi)];
                end
            else
                flat = [flat nan];
                d_pxl = [d_pxl nan];
                x_pxl = [x_pxl nan];
                y_pxl = [y_pxl nan];
                phi = [phi nan];
            end
        else
            flat = [flat nan];
            d_pxl = [d_pxl nan];
            x_pxl = [x_pxl nan];
            y_pxl = [y_pxl nan];
            phi = [phi nan];
        end
    else
        flat = [flat nan];
        d_pxl = [d_pxl nan];
        x_pxl = [x_pxl nan];
        y_pxl = [y_pxl nan];
        phi = [phi nan];
    end
end

clear freeEye_x freeEye_y

%% Correct for the movement of CR on the eye surface

% replace bad CR tracking with nan
tr = CR(:,1)>roi_CR_x(2) | CR(:,1)<roi_CR_x(1);
tr = tr | (CR(:,2)>roi_CR_y(2) | CR(:,2)<roi_CR_y(1));
CR(tr,1) = nan;
CR(tr,2) = nan;

clear tr

% filter
CR = filterCR(CR);

% determine the systematic relationship between the pupil center and the CR
% along the x-axis
mdl = fitlm(x_pxl,CR(:,1),'linear');

aX = mdl.Coefficients.Estimate(2); % slope

bX = mdl.Coefficients.Estimate(1); % intercept

mdl = fitlm(y_pxl,CR(:,2),'linear');

aY = mdl.Coefficients.Estimate(2); % slope

bY = mdl.Coefficients.Estimate(1); % intercept

clear mdl

% determine the 50%ile of CR position to use as the landmark
cr50_x = nanmedian(CR(:,1));
cr50_y = nanmedian(CR(:,2));

% corrected CR x-axis
CR_x = CR(:,1)' + (cr50_x - (aX * x_pxl + bX));

% uncorrected CR y-axis
CR_y = CR(:,2)' + (cr50_y - (aY * y_pxl + bY));

clear CR cr50_x cr50_y aX bX aY bY

% correct for movement using CR
x_pxl = x_pxl - CR_x;
y_pxl = y_pxl - CR_y;

clear CR_x CR_y

%% figure out the center
len = length(0:0.01:1);
dk = [0 107 179]/255;
lt = [165 222 255]/255;
colors_p = [linspace(lt(1),dk(1),len)', linspace(lt(2),dk(2),len)', linspace(lt(3),dk(3),len)'];

x_range = [-60 40];
y_range = [-80 20];

F = zeros(101,101);
P = zeros(101,101);
count = zeros(101,101);
for i=1:length(flat)
    x = round(x_pxl(i));
    y = round(y_pxl(i));
    if x >= x_range(1) && x <= x_range(2) && y >= y_range(1) && y <= y_range(2)
        x = x + abs(x_range(1))+1;
        y = y + abs(y_range(1))+1;
        F(y,x) = F(y,x)+flat(i);
        P(y,x) = P(y,x)+phi(i);
        count(y,x) = count(y,x)+1;
    end
end
avFlat = F./count;
avPhi = P./count;

figure;s=surf(x_range(1):x_range(2),y_range(1):y_range(2),avFlat);
s.EdgeColor = 'none';
s.FaceAlpha = 1;
xlabel('pupil center x (pxl)')
ylabel('pupil center y (pxl)')
title('Flattening')
set(gca,'YDir','reverse')
pbaspect([1 1 1])
view(2)
caxis([0 0.2])
figure;s=surf(x_range(1):x_range(2),y_range(1):y_range(2),avPhi);
s.EdgeColor = 'none';
s.FaceAlpha = 1;
xlabel('pupil center x (pxl)')
ylabel('pupil center y (pxl)')
title('Orientation')
set(gca,'YDir','reverse')
pbaspect([1 1 1])
view(2)
caxis([0 90])
freq = count/sum(sum(count));
freq(isnan(avFlat)) = nan;
figure;s=surf(x_range(1):x_range(2),y_range(1):y_range(2),freq);
s.EdgeColor = 'none';
s.FaceAlpha = 1;
xlabel('pupil center x (pxl)')
ylabel('pupil center y (pxl)')
title('Frequency')
set(gca,'YDir','reverse')
pbaspect([1 1 1])
view(2)

save('eyeProfile.mat','avFlat','avPhi','freq','x_range','y_range','-mat');

% correct the coordinate using center
prompt = {'Enter x-coordinate of the center:','Enter y-coordinate of the center:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','0'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

x_center = str2num(answer{1});
y_center = str2num(answer{2});

x_pxl = x_pxl - x_center;
y_pxl = y_pxl - y_center;

% convert d_pxl into d_mm
d_mm = d_pxl * mean([conv_x conv_y]);

% determine Rp (distance from the center of the eyeball to the pupil) from d_mm
Rp = d_mm * slope_Rp + int_Rp;

% convert x_pxl to x_mm
x_mm = x_pxl * conv_x;

% determine azimuth in degrees
azim = asind(x_mm./Rp);

% convert y_pxl to y_mm
y_mm = y_pxl * conv_y;

% determine elevation in degrees
elev = asind(y_mm./Rp);

%% Determine trigger timing

% Specify file path
trigger_fname = strcat(exp_id,'_DIN04.mat');

% From DIN_04 channel, obtain trigger timings
dg04read = load(trigger_fname,'-mat');
din04 = dg04read.data;

% Get trial timings in samples
triggertiming = Get_timing_digital(din04,'both'); % first column, trial start, second column, trial end

% Convert to ms
triggertiming = triggertiming * sample_duration;

clear dg04read din04;

% Read timestamps from Raspberry pi
rasp_fname = strcat(exp_id,'_timestamps.csv');
Rasp_time = readtable(rasp_fname);
rasp_time = table2array(Rasp_time)/1e3; % ms

% % time difference between Intan and Raspberry pi
% time_diff = triggertiming(1,1)-rasp_time(1);
% 
% % Convert acquisition timing to Intan timing (ms)
% camtiming = round(rasp_time(:,1)'+time_diff);

% time difference between frame grabbing and TTL on RPi (ms)
time_diff = rasp_time(:,2)-rasp_time(:,1);

% Convert acquisition timing to Intan timing (ms)
if length(time_diff) == size(triggertiming,1)-1 % sometimes the last trigger is missed on raspberry pi
    camtiming2 = round(triggertiming(1:end-1,1)-time_diff)';
    azim = azim(1:end-1);
    elev = elev(1:end-1);
    d_mm = d_mm(1:end-1);
else
    camtiming2 = round(triggertiming(:,1)-time_diff)';
end

% 1st row: azimuth, 2nd: elevation, 3rd: radius, 4th: timing
azim = -azim; % convert nasal to negative
eyeinfo = [azim;elev;d_mm;camtiming2];

end

function A = filterCR(A)

% A is a two column vector, 1st column is CR_x, 2nd is CR_y

[elem,len] = runlen(isnan(A(:,1)'));
len_nonnan = len(1:2:end);

for i=1:length(len_nonnan)-1 % look for two consecutive counts >thr
    if len_nonnan(i)>30
        if len_nonnan(i+1)>30
            % Check if the length of NaN flanked by i and i+1 is 1
            ii = i*2; % index of the Nan flanked by i and i+1
            if len(ii)==1
                % if there is only 1 NaN in between, then filter
                % figure out the index of the NaN in A
                idx = sum(len(1:ii));
                A(idx,:) = mean([A(idx-1,:);A(idx+1,:)],1);
            end
        end
    end
end

clear elem len len_nonnan idx
end
