function limits = find_lims(datain, argin2)

% find_lims: finds the minimum and maximum values of datain.
%   limits = find_lims(datain) will return a 1x2 double of the format
%   [minimum, maximum]. If datain is a matrix or a cell array, it will find
%   the limits of all values.
%   if datain in a structure, find_lims finds the limits of elements in
%   datain(:).data. 
%
%   limits = find_lims(datain, argin2) where argin2 is a string and datain
%   is a struct finds the limits of datain(:).(argin) instead.

% Written by Rannee Li, Jan 2014. 
% Last updated Jul 2014 (help revamped, changed eval to dynamic fieldname)

% if datain is a struct, argin2 should be the fieldname. Dafault: 'data' 
    
if isstruct(datain) || iscell(datain) % have to iterate
    % find number of iterations
    nIt = numel(datain);
    
    % preallocate
    mins = zeros(nIt, 1);
    maxs = zeros(nIt, 1);
    
    % loop through
    for k = 1:nIt
        if isstruct(datain)
            if nargin == 1
                warning('Looking for field ''data''.')
                argin2 = 'data';
            end
            dat = datain(k).(argin2);
        else % if it is a cell
            dat = datain{k};
            if nargin > 1
                warning('Extra arguments ignored.')
            end
        end
        [tmpl, tmpu] = lim_mat(dat);
        
        mins(k) = tmpl;
        maxs(k) = tmpu;
    end

else 
    % directly find min and max
    [mins, maxs] = lim_mat(datain);
    
    if nargin > 1
        warning('Extra arguments ignored.')
    end    

end


limits = [min(mins), max(maxs)];

% subfunction.
function [bl, bu] = lim_mat(dat)

    bl = min(real(dat(:)));
    bu = max(real(dat(:)));
    
    if isinf(bl) || isinf(bu)
        dat(isinf(dat)) = NaN;
        warning('inf detected and ignored')
        
        bl = min(dat(:));
        bu = max(dat(:));
    end
    
    if ~isreal(bl) || ~isreal(bu)
        warning('complex part detected. considering real parts only')
        
        bl = min(real(dat(:)));
        bu = max(real(dat(:)));
    end
