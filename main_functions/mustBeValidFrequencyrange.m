function mustBeValidFrequencyrange(frequencyrange)
    if isnumeric(frequencyrange)
        mustBeInRange(frequencyrange,3,90)
    elseif ischar(frequencyrange)
        mustBeMember(frequencyrange,['theta', 'alpha', 'beta low','beta high','beta total','gamma broad'])
    else
        eidType = 'mustBeValidFrequencyrange:notValidFrequencyrange';
        msgType = 'Input must be a valid frequency range';
        error(eidType,msgType)
    end
end