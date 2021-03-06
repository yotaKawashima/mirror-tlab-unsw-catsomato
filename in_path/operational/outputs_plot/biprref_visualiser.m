function [p, FV] = biprref_visualiser(data, labels, spatialconfig)


nChan = length(data);

vertices = make_vertices(spatialconfig(2), spatialconfig(1));

% preallocate
faces = zeros(nChan, 4);

for l = 1:nChan
    face = add_face(labels(l, 1:2), vertices, labels(l, 3));
    
    faces(l, :) = face;
end

FV.vertices = vertices;
FV.faces = faces;

clear faces vertices

% draw
% figure
p = patch(FV);


% color
clear cdata 
%set(gca,'CLim', plotcfg.limits)
cdata = data;
set(p,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled')
axis tight
axis equal
axis ij


if nargout < 2
    clear FV
end
if nargout < 1
    clear p
end

% set(gca,'YDir','reverse');

function vertices = make_vertices(spatial_width, spatial_height)
X_gv = 0.5:0.5:(spatial_width + 0.5);
Y_gv = 0.5:0.5:(spatial_height + 0.5);

k = 1;
for x = 1:spatial_width
    for y = 1:spatial_height
        vertices(k, :) = [x, y];  %#ok<AGROW>
        k = k+1;
    end
end

for a = X_gv(1)
    for b = Y_gv(3:2:end-2)
        vertices(k, :) = [a, b]; %#ok<AGROW>
        k = k+1;
    end
end

for c = X_gv(3:2:end-2)
    for d = Y_gv(1:2:end)
        vertices(k, :) = [c, d]; %#ok<AGROW>
        k = k+1;
    end
end

for a = X_gv(end)
    for b = Y_gv(3:2:end-2)
        vertices(k, :) = [a, b]; %#ok<AGROW>
        k = k+1;
    end
end

function face = add_face(anchors, vertices, anchordir)
% find the two vertices given by the anchors
V1 = anchors(1);
V3 = anchors(2);

% find the midpoint between the two anchor points
midps = (vertices(V1, :) + vertices(V3, :)) /2;

% find co-ordinates of the stabilising points
if anchordir == 'h'
    point2 = [midps(1) midps(2)+0.5];
    point4 = [midps(1) midps(2)-0.5];
else
    point2 = [midps(1)+0.5 midps(2)];
    point4 = [midps(1)-0.5 midps(2)];
end

% find values of V2 and V4
for k = anchors(2):length(vertices)
    if vertices(k, :) == point2
        V2 = k;
    end
    if vertices(k, :) == point4
        V4 = k;
    end
end

face = [V1 V2 V3 V4];


