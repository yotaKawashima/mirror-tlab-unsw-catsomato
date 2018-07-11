function [val, ind] = parfind_closest(A, target)

if size(A, 1) > 1
    if size(A, 2) > 1
        error('A has more than one row.')
    else
        A = A';
    end
end

if size(target, 2) > 1
    if size(target, 1) > 1
        error('A has more than one row.')
    else
        target = target';
    end
end


% repeat A
nTarg = numel(target);
A1 = repmat(A, [nTarg, 1]);
target1 = repmat(target, [1, size(A, 2)]);

[~, ind] = min(abs(A1 - target1), [], 2);
val = A1(ind);

