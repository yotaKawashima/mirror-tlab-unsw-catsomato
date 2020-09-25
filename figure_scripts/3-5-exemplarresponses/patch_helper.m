function [x_data, y_data] = patch_helper(channel_number, relabels)

DEBUG = false;

n_chans = max(relabels(:));
n_perrow = sqrt(n_chans);

c_row = relabels(channel_number, :);
[r1, c1] = ph_find_rc(c_row(1), n_perrow);
[r2, c2] = ph_find_rc(c_row(2), n_perrow);

if DEBUG
    fprintf('row: %i; %i, %i, %i\n', channel_number, ...
        c_row(1), c_row(2), c_row(3))
    fprintf('%i, %i, %i, %i\n', r1, c1, r2, c2)
end

% vertical is 0, horizontal is 1

% assuming no malformed inputs
if c_row(3) == 0
    % vertical. 
    % top and bottom points are defined, c1==c2
    x_data = [c1; c1+0.5; c2; c1-0.5];
    y_data = [r1; r1-0.5; r2; r1-0.5];
        
else
    x_data = [c1; c1-0.5; c2; c1-0.5];
    y_data = [r1; r1+0.5; r2; r1-0.5];   
    
end

function [r, c] = ph_find_rc(chan, n_perrow)

c = floor(chan/n_perrow) + 1;
r = mod(chan, n_perrow);

if r == 0
    r = n_perrow;
    c = c-1;
end

