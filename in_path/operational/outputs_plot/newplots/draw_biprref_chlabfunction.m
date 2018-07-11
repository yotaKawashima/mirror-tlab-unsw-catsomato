function label_out = draw_biprref_chlabfunction(label_in)
% re-labels channels from the output by the function bipolreref to the
% input required by draw_biprref
% mainly used as a helper function

% Written by Rannee Li, Jul 2015

% preallocate
label_out = zeros(numel(label_in), 3);

for k = 1:numel(label_in)
    % find the channel numbers
    label_out(k, 1) = str2double(label_in{k}(10:12)); 
    label_out(k, 2) = str2double(label_in{k}(16:18));
    
    % find the direction
    if label_in{k}(6)=='v'
        label_out(k, 3) = 0;
    elseif label_in{k}(6)=='h'
        label_out(k, 3) = 1;
    end
end
