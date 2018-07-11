function [val, ind] = find_closest(A, target)

% find_closest: find the closest value to target in the matrix A
% 
% Usage: [val, ind] = find_closest(A, target)
%        val is the closest value. its index is ind
% Will eventually be upgraded to full find capabilities

if iscell(A)
    warning('Cell detected as comparison. Defaulting to first element of cell.')
    A = A{1};
end

[~, ind] = min(abs(A - target));
val = A(ind);

if nargout < 2
    clear ind
end