function timing = Get_timing_digital(dg_array,specifier)
% This function will obtain the onsets and offsets of digital pulses and
% return their indices. specifier specifies whether to return onset,
% offset, or both. In case of both, first column, onset, second column,
% offset.
% LOG
% 3/10/16 change vec2mat to reshape

dgDiff = diff(dg_array);

ind = find(dgDiff); % return indices of 1s and -1s.

%temp = vec2mat(ind,2);
temp = reshape(ind,2,[])';
temp(:,1)=temp(:,1)+1;

if strcmp(specifier,'onset')
    timing = temp(:,1);
elseif strcmp(specifier,'offset')
    timing = temp(:,2);
elseif strcmp(specifier,'both')
    timing = temp;
else
    error('invalid specifier');
end

