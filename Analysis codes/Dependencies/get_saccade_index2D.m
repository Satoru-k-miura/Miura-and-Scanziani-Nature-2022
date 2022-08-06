function returned_values = get_saccade_index2D(eyetrace, params)
% This function takes a 1-D array of eye positions and will detect all saccade
% events that fullfil the speed and magnitude criteria provided. Artefact
% eye movements will also be deleted if detected. Returns a 2-D array of
% indexes, where the 1st column is the onset and the 2nd column is the
% offset of saccades and rows correspond to saccade events
% speed_thr: degrees per second, fr: frame rate, Hz, mag: degrees

% editing note:
% 11/19/2015 Added eyetrace_diff_shift2
% 1/30/2016 Added variables direction and speed_thr_max
% 10/31/2016 Added directions
% 11/1/2016 Added amp
% 8/2/2021 converted to 2D eyetrace with time index
%   eyetrace is 4 x N where 1st row: azimuth, 2nd row: elevation, 3rd row:
%   diameter, 4th row: time (ms)

%% Parameters

stationary = 1.5; % eye movement allowed during stationary period
initial_mv = 200; % the speed (deg/s) that the first frame must exceed 
trj_thr = 45; % limit for the change in trajectory over the course of one saccade (deg)
speed_thr = params.speed_thr;
speed_thr_max = 2000; % saccades very unlikely to exceed 2000 deg/sec
fr=params.fr;
mag = params.mag; % minimum magnitude of saccades to consider
stp = 300/(1e3/fr); % search step, every 300ms
if isfield(params,'edges')
    sBefore = abs(round(params.edges(1)/(1e3/fr)));
    sAfter = abs(round(params.edges(end)/(1e3/fr)));
else
    sBefore = 3;
    sAfter = 2;
end

%% Figure out the indexes in the 2-D array, where there is a change in position

% calculate movement between each frame (deg)
eyetrace_sp = [];
for i=1:size(eyetrace,2)-1
    [arclen,~] = distance(eyetrace(2,i),eyetrace(1,i),eyetrace(2,i+1),eyetrace(1,i+1));
    eyetrace_sp(i) = arclen/(eyetrace(4,i+1)-eyetrace(4,i))*1e3; %deg/s
end

eyetrace_sp_shift = circshift(eyetrace_sp,[0 -1]); % Shift the difference
eyetrace_sp_shift2 = circshift(eyetrace_sp_shift,[0 -1]); % Shift the difference again

% find the first index that exceeds the speed threshold
idx = find(eyetrace_sp>=speed_thr & eyetrace_sp<speed_thr_max,1);

%% Determine the saccade onset and offset indexes

indexes = []; % Initialize the return value
amp=[];
pos=[];
eyeX = [];
eyeY = [];

while ~isempty(idx) % keep searching for the saccade as long as idx is found
    if idx>sBefore && idx<length(eyetrace)-sAfter % Make sure the index is not too close to the edge of the recording
        % To avoid noise, make sure the points around the speed change are
        % NOT close to each other
        if distance(eyetrace(2,idx),eyetrace(1,idx),eyetrace(2,idx+2),eyetrace(1,idx+2)) > stationary && ~any(isnan(eyetrace(1,idx-sBefore:idx+sAfter))) && ~any(isnan(eyetrace(2,idx-sBefore:idx+sAfter))) && isreal(eyetrace(1,idx-sBefore:idx+sAfter)) && isreal(eyetrace(2,idx-sBefore:idx+sAfter))
            % by comparing the vector direction and the speed to the 5
            % frames before idx, determine the first frame of the saccade
            y = diff(eyetrace(2,idx-5:idx+1));
            x = diff(eyetrace(1,idx-5:idx+1));
            theta = atan2d(y,x); % vector angle for each frame including one following idx
            % the trajectory should be consistent with the idx frame and
            % the speed also must exceed initial_mv in order to be the
            % first frame of saccade
            [~,len] = runlen((abs(theta(end)-theta)<trj_thr) & (eyetrace_sp(idx-5:idx)>initial_mv));
            if len(end)<=3 % we allow up to 2 frames before
                idx_st = idx-(len(end)-1); % this is the beginning frame of the saccade
                
                % Look for the end of the saccade
                % 1. the speed slows down for one frames OR
                % 2. change in the course
                y = diff(eyetrace(2,idx_st:idx_st+10));
                x = diff(eyetrace(1,idx_st:idx_st+10));
                theta = atan2d(y,x);
                idx_end = idx_st+find(eyetrace_sp(idx_st+1:idx_st+10)<100 | abs(theta(1)-theta)>trj_thr,1);
                
                if isempty(idx_end) || idx_end>size(eyetrace,2)-10 % if not found, assume ~40ms for now
                    idx_end = idx_st+round(40/(1e3/fr));
                end
                % Finally, make sure that the movement exceeds the
                % minimum magnitude and that idx_st and idx_end are no more than
                % 100ms away to avoid artefacts
                if (idx_end-idx_st+1)*(1e3/fr)>100
                    % if too far, assume ~40ms for now
                    idx_end = idx_st+round(40/(1e3/fr));
                end

                sacAmp = distance(eyetrace(2,idx_st),eyetrace(1,idx_st),eyetrace(2,idx_end),eyetrace(1,idx_end));
                y = diff(eyetrace(2,idx_st:idx_st+2));
                x = diff(eyetrace(1,idx_st:idx_st+2));
                theta = atan2d(y,x);
                if sacAmp > mag && idx_end-idx_st>1 && (abs(max(theta)-min(theta))<trj_thr || abs(max(theta)-min(theta))>360-trj_thr) ...
                        && distance(eyetrace(2,idx_st),eyetrace(1,idx_st),eyetrace(2,idx_st+2),eyetrace(1,idx_st+2)) > stationary ...
                        && distance(eyetrace(2,idx_st+1),eyetrace(1,idx_st+1),eyetrace(2,idx_st+2),eyetrace(1,idx_st+2)) > stationary...
                        && distance(eyetrace(2,idx_st-1),eyetrace(1,idx_st-1),eyetrace(2,idx_st+1),eyetrace(1,idx_st+1)) > stationary...
                        && distance(eyetrace(2,idx_st-1),eyetrace(1,idx_st-1),eyetrace(2,idx_st),eyetrace(1,idx_st)) < 5
                    indexes = [indexes; [idx_st idx_end]];
                    amp = [amp; sacAmp];
                    pos = [pos; eyetrace([1 2],idx_st)'];
                    eyeX = [eyeX; eyetrace(1,idx_st-sBefore:idx_st+sAfter)];
                    eyeY = [eyeY; eyetrace(2,idx_st-sBefore:idx_st+sAfter)];

                    idx = find(eyetrace_sp(idx+stp:end)>=speed_thr & eyetrace_sp(idx+stp:end)<speed_thr_max,1)+idx+stp-1;   
                else % move on to the next one
                    idx = find(eyetrace_sp(idx+10:end)>=speed_thr & eyetrace_sp(idx+10:end)<speed_thr_max,1)+idx+9;       
                end 
            else % the saccade beginning cannot be defined clearly, move on to the next one
                idx = find(eyetrace_sp(idx+1:end)>=speed_thr & eyetrace_sp(idx+1:end)<speed_thr_max,1)+idx;
            end           
        else % the change could be an artefact/noise or may include nans around the time
            idx = find(eyetrace_sp(idx+2:end)>=speed_thr & eyetrace_sp(idx+2:end)<speed_thr_max,1)+idx+1; % update idx
        end
    elseif idx<=sBefore % if idx is too close to the start, look for the next one
        idx = find(eyetrace_sp(idx+1:end)>=speed_thr & eyetrace_sp(idx+1:end)<speed_thr_max,1)+idx; % If idx is at the beginning or the end of eyetrace, look for the next one
        % idx will return an empty value if no more idx is found
    else % no more, stop searching
        break;
    end
end
    
returned_values.indexes = indexes;
returned_values.eyetrace = eyetrace;
returned_values.amp = amp;
returned_values.pos = pos;
returned_values.eyeTraceX = eyeX;
returned_values.eyeTraceY = eyeY;
