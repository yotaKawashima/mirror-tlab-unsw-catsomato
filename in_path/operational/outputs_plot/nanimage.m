function nanimage(data)
% nanimage: plots an image where there is NaN data
%
% nanimage(data) plots data using pcolor to create a imagesc-like plot.

% from https://au.mathworks.com/matlabcentral/answers/81938-set-nan-as-another-color-than-default-using-imagesc

[nr,nc] = size(data);

pcolor([data nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');