figure(1); clf
[H, S] = venn(10*ones(1, 3), [3*ones(1, 3), 9], 'FaceColor', {'none', 'none', 'none'});
axis equal


zonelab = {'cF1', 'cF2', 'cF1F2', 'CF1\capcF2', 'cF1\capcF1F2', 'cF2\capcF1F2', 'cF1\capcF2\capcF1F2'};
for i = 1:7
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), zonelab{i})
end

