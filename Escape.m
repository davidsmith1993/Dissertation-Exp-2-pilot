function escapeSequence(key,location)
%function out = escapeSequence(key,location)
%
%5/10/05    swe     modified version of g. baldwin's escapeSequence.m
%


if isequal(upper(key), 'Q')
    clean_up
    while(kbCheck);end %Wait for all keys to be released.
    error_text = ['user quit the experiment in ' location];
    error(error_text); %fprintf('\n%s',error_text);
end