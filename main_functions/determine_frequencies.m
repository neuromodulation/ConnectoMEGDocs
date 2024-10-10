function [freq]=determine_frequencies(frequencyrange, frequencies)
    % frequencyrange = 'all' (3:87 or 3:36), [f1 : f2], or 'theta', 'alpha', 'beta_low','beta_high','beta_total','gamma_narrow'
    % freq = struct()
    % freq.frequency_indices
    % freq.frequencybandnames
    % freq.frequencyband
if nargin < 2
    load('standard_frequencies.mat','frequencies')
end
if isnumeric(frequencyrange)
    assert(frequencyrange(1)>= min(frequencies))
    assert(frequencyrange(2)<= max(frequencies))
    freq.frequency_indices{1} = frequencyrange - 3;
    freq.frequencybandnames={'custom'};
    freq.frequencyband{1} = frequencyrange;
elseif ischar(frequencyrange)
    freq.frequencybandnames={'theta', 'alpha', 'beta low','beta high','beta total'};
    freq.frequencyband{1} = 4:7;
    freq.frequencyband{2} = 7:12;
    freq.frequencyband{3} = 13:21;
    freq.frequencyband{4} = 21:35;
    freq.frequencyband{5} = 13:35;
    if (60 <= max(frequencies) &&  max(frequencies) <= 90)
        freq.frequencybandnames{6} = 'gamma broad';
        freq.frequencyband{6} = 60:90;        
    end
    for i = 1:length(freq.frequencyband)
        freq.frequency_indices{i} = freq.frequencyband{i}-3;
    end
    if strcmp(frequencyrange,'all')
        freq.index=[1:32];
        return;
    else
        ind=ismember(freq.frequencybandnames,frequencyrange);
        freq.frequencybandnames = {frequencyrange};
        temp   = freq.frequencyband{ind};
        freq= rmfield(freq, 'frequencyband');
        freq.frequencyband{1}   = temp;
        temp   = freq.frequency_indices{ind};
        freq= rmfield(freq, 'frequency_indices');
        freq.frequency_indices{1}   = temp;
        freq.index=temp;
    end
          
else
    error("frequencyrange undefined")
end

% freq.index=[];
% for i=freq.frequency_indices
%     freq.index=[freq.index, i{:}]
% end
end