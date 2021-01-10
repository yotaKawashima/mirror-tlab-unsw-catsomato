% Model fit for top 10% channels
% 
% This script is for those who do not have access to a supercomputer. I have
% not tried this and am not sure how long it would take to finish. I
% recommend you to try only one channel fitting first by commenting out
% for-loops and estimate the time.

% # of chs used for model fitting for each somatosensory area
num_chs = [156, 101]; 

% Start model fitting.
for area_id = 1:2
    num_ch = num_chs(area_id);
    for channel_id = 1:num_ch
        model_rect_fit_massive_v2(channel_id, area_id);
        model_sq_fit_massive_v2(channel_id, area_id);
    end % channel_id = 1:num_ch
end % area_id = 1:2
