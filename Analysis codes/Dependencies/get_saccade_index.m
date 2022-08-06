function returned_values = get_saccade_index(eyetrace, params)
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

%% Parameters

stationary = 0.7; % eye movement allowed during stationary period
speed_thr = params.speed_thr;
speed_thr_max = 2000; % saccades very unlikely to exceed 2000 deg/sec
fr=params.fr;
mag = params.mag; % minimum magnitude of saccades to consider
if isfield(params,'edges')
    sBefore = abs(round(params.edges(1)/(1e3/fr)));
    sAfter = abs(round(params.edges(end)/(1e3/fr)));
else
    sBefore = 3;
    sAfter = 2;
end


%% Figure out the indexes in the 1-D array, where there is a change in position

eyetrace_diff = diff(eyetrace); % first order difference
eyetrace_diff_shift = circshift(eyetrace_diff,[0 -1]); % Shift the difference
eyetrace_diff_shift2 = circshift(eyetrace_diff_shift,[0 -1]); % Shift the difference again

idx = find(abs(eyetrace_diff)>speed_thr/fr,1); % this index corresponds to the frame just before the change in the original eyetrace

%% Determine the saccade onset and offset indexes

indexes = []; % Initialize the return value
directions = [];
amp=[];
pos=[];
eyeT = [];

while ~isempty(idx) % keep searching for the saccade as long as idx is found
    if idx>sBefore && idx<length(eyetrace)-sAfter %idx>3 && idx<length(eyetrace)-2 % Make sure the index is bigger than 2 but smaller than length(eyetrace)-2
        if abs(eyetrace(idx+2)-eyetrace(idx)) > stationary % If two points AROUND the change are NOT close to each other
            % Look for the first index after idx where
            % diff<stationary
            idx2 = idx-1+find(abs(eyetrace_diff(idx:end))<stationary & abs(eyetrace_diff_shift(idx:end))<stationary & abs(eyetrace_diff_shift2(idx:end))<stationary,1);
            if ~isempty(idx2) && idx2<length(eyetrace) 
                % Finally, make sure that the movement exceeds the
                % minimum magnitude and that idx and idx2 are no more than
                % 100ms away to avoid artefacts
                direction = sign(eyetrace(idx2)-eyetrace(idx))==sign(eyetrace(idx+1)-eyetrace(idx)) && sign(eyetrace(idx2)-eyetrace(idx))==sign(eyetrace(idx+1)-eyetrace(idx)) && sign(eyetrace(idx2)-eyetrace(idx))==sign(eyetrace(idx+2)-eyetrace(idx))...
                    && abs(eyetrace(idx+2)-eyetrace(idx)) > abs(eyetrace(idx+1)-eyetrace(idx)) && abs(eyetrace(idx+3)-eyetrace(idx)) > abs(eyetrace(idx+2)-eyetrace(idx));
                miv = max([1 idx-0.2*fr]);
                mav = min([idx+0.2*fr length(eyetrace)]);
                if abs(mean([eyetrace(idx-1) eyetrace(idx)]) - mean([eyetrace(idx2) eyetrace(idx2+1)])) > mag && idx2 - idx < 0.1*fr && eyetrace(idx+1) - eyetrace(idx) < speed_thr_max/fr && direction && ~any(isnan(eyetrace(miv:mav)))
                    indexes = [indexes; [idx idx2]];
                    [~,tempIdx] = max(abs(eyetrace(idx:idx2)-eyetrace(idx)));
                    amp = [amp; eyetrace(idx+tempIdx-1)-eyetrace(idx)]; %[amp; eyetrace(idx2)-eyetrace(idx)];
                    pos = [pos; eyetrace(idx+tempIdx-1)]; %[pos; eyetrace(idx2)];
                    eyeT = [eyeT; eyetrace(idx-sBefore:idx+sAfter)];
                end
            else % If it does not exist, end searching
                break;
            end
            
            idx = find(abs(eyetrace_diff(idx2:end))>speed_thr/fr,1)+idx2-1; % update idx
            
        else
            eyetrace(idx+1) = mean([eyetrace(idx) eyetrace(idx+2)]); % if the change is detected to be an artefact, replace the value with the average of the surrounding points
            idx = find(abs(eyetrace_diff(idx+2:end))>speed_thr/fr,1)+idx+1; % update idx
        end
        
    elseif idx<=sBefore % if idx is smaller than or equal to three, look for the next one
        idx = find(abs(eyetrace_diff(idx+1:end))>speed_thr/fr,1)+idx; % If idx is at the beginning or the end of eyetrace, look for the next one
        % idx will return an empty value if no more idx is found
    else % if idx >= length(eyetrace)-2, stop searching
        break;
    end
end
    
returned_values.indexes = indexes;
returned_values.eyetrace = eyetrace;
returned_values.amp = amp;
returned_values.directions = sign(amp);
returned_values.pos = pos;
returned_values.eyeTrace = eyeT;
