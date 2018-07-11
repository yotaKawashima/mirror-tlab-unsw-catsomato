function [frow, fcol] = a2c_depfreq(data)

% helper function to ANOVA2_CALL and ANP_FIGURES
%   For more information see usage in above call.

% finds the frequencies used in the paradigm
% assumes frequencies and magnitudes between 000 and 999

% Change log:
% 10 Jul 17: automate finding the indicies rather than hardcoding it.

% find indices and error check
f_ind = strfind(data.conditions{1}, 'F');
a_ind = strfind(data.conditions{1}, 'A');

if numel(f_ind) ~= 2 || numel(a_ind) ~= 2
    error('Unconventional naming, cannot process.')
end

% read in the two frequencies.
f1 = data.conditions{1, 1}(f_ind(1)+1:f_ind(1)+3);  % 20/23Hz
f2 = data.conditions{1, 1}(f_ind(2)+1:f_ind(2)+3); % 200Hz

% read the value of the f1 in the top row
% if the value of the frequency in f1 is the same across a row of
% conditions (we know that it is in an increasing grid), then it must be
% the row variable
isf1row = strcmp(data.conditions{1, end}(a_ind(1)+1:a_ind(1)+3), ...
    data.conditions{1, 1}(a_ind(1)+1:a_ind(1)+3));

if isf1row
    frow = f1;
    fcol = f2;
else
    frow = f2;
    fcol = f1;
end