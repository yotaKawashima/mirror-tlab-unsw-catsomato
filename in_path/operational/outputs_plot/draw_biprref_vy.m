function draw_biprref_vy(values, ch_label, spatial_cfg, clim)

% draw_biprref: Calls the patch function to plot bipolar re-referenced data
%   draw_biprref(values, ch_label, spatial_cfg) draws the image.
%
%   draw_biprref(values, ch_label, spatial_cfg, clim) uses clim as the
%   limits of the colormap.
%
%   p = draw_biprref(...) returns the handle of the patch.
%
%   [p, FV] = draw_biprref(...) returns the handle and the FV data of the
%   patch.
%
%   INPUT FORMAT 
%   values      : a 1D matrix holding the values of data for ONLY the 
%                 bipolar-rereferenced channels.
%   ch_label    : numel(values)*3 matrix. columns 1 and 2 contain 
%                 the numbers of the channels the re-referenced channel is 
%                 made of. column 3 is the direction 
%                 (0 for vertical, 1 for horizontal)
%   spatial_cfg : [rows, columns] of the spatial arrangement (bipolar ch)
%   clim        : color lim
%

% Written by Rannee Li, Jan 2015 and Modiffied by Yota, Apr 2020
% Change log:
%   Apr 2020 (Yota): Inf patch will be plotted first with gray color. 

% Preallocate patch
% call subfunctions
% xy coordinate of vertices (including electrodes and supporting points)
vertices = make_vertices(spatial_cfg(2), spatial_cfg(1));
% a set of 4 vertex ids for a diamond representing a bipolar channel.
faces = all_faces(ch_label, vertices);

% Copy the original in order to refer to its size when plotting.
faces_ = faces;
%values_ = values;

% Then remove patches whose value is nan.
% First convert inf into nan
if ~isempty(find(isinf(values), 1))
    warning('Infs detected in input values. Assigning +/-Inf = NaN.')
    values(isinf(values)) = NaN;
end

% If NaN is included in values, remove the corresponding patches.
if ~isempty(find(isnan(values), 1))
    warning(['NaNs detected in input values.', ...
            ' Do not show the corresponding patches.'])

    % remove faces which have nan value.
    faces = faces(~isnan(values),:);
    
    % Remove NaN from a matrix 'values'.
    values = values(~isnan(values));
end

% Preallocate patches with the same color first. 
% place into structure as requested by patch

hold on 
FV1.vertices = vertices;
FV1.faces = faces_;
p1 = patch(FV1);
% Set color for NaN. 
set(p1,'FaceColor',[200 200 200]/256, 'LineWidth', 0.25, 'Edgecolor', 'k')
clear cdata1


% Draw patches over the preallocated ones if not NaN.
if ~isempty(faces)
    FV2.vertices = vertices;
    FV2.faces = faces; % NaNs are removed from faces. 
    p2 = patch(FV2);
    
    % Set color accordingly.
    cdata2 = values; % NaNs are removed from values. 
    set(p2,'FaceColor', 'flat', 'LineWidth', 0.25, ...
        'FaceVertexCData', cdata2, 'CDataMapping', 'scaled')
    clear p2
end

set(gca,'CLim', clim)
axis tight
axis equal
axis ij
hold off

clear p1
end

%% END MAIN FUNCTION CALL. SUBFUNCTIONS FOLLOW

%% MAKE_VERTICES
% creates all vertices needed for the patch.
function vertices = make_vertices(spatial_width, spatial_height)

X_gv = 0.5:0.5:(spatial_width + 0.5);
Y_gv = 0.5:0.5:(spatial_height + 0.5);

Y_use = Y_gv(1:2:end);
nummidpts = numel(Y_use);

% make the points that represent the channels
tmp = {1:spatial_height, 1:spatial_width};
[o1, o2] = ndgrid(tmp{:});
vert_1 = [o2(:) o1(:)];

% make vertices in column 1
vert_2 = [X_gv(1)* ones(nummidpts-2, 1), Y_gv(3:2:end-2)'];

% make vertices in column 2 - end-1
tmp = {Y_use, X_gv(3:2:end-2)};
[o1, o2] = ndgrid(tmp{:});
vert_3 = [o2(:) o1(:)];

% make vertices in end column
vert_4 = [X_gv(end)* ones(nummidpts-2, 1), Y_gv(3:2:end-2)'];

vertices = [vert_1; vert_2; vert_3; vert_4];
end

%% ALL_FACES
% maps all of the bipolar channels onto the vertices
function faces = all_faces(ch_label, vertices)
% find the midpoints of all the channels
midps = (vertices(ch_label(:, 1), :) + vertices(ch_label(:, 2), :)) /2;

% find all the horizontally re-referenced channels
horzinds = find(ch_label(:, 3));
vertinds = find(~ch_label(:, 3));

% calculate the location of the other two points
% split into horx and vert as this affects the result
point2 = zeros(size(ch_label, 2), 2); % preallocate
point4 = point2;

point2(horzinds, :) = [midps(horzinds, 1) midps(horzinds, 2)+0.5];
point4(horzinds, :) = [midps(horzinds, 1) midps(horzinds, 2)-0.5];

point2(vertinds, :) = [midps(vertinds, 1)+0.5 midps(vertinds, 2)];
point4(vertinds, :) = [midps(vertinds, 1)-0.5 midps(vertinds, 2)];

% locate the actual vertex numbers
[~, f_2] = ismember(point2, vertices, 'rows');
[~, f_4] = ismember(point4, vertices, 'rows');

faces = [ch_label(:, 1) f_2 ch_label(:, 2) f_4];
end
