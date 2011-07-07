clear all;
close all;

alpha = 12.5;

G = zeros(9,3);
G(1,:) = [0 0 0];
G(2,:) = [1 1 1];
G(3,:) = [1 1 -1];
G(4,:) = [1 -1 1];
G(5,:) = [-1 1 1];
G(6,:) = [1 -1 -1];
G(7,:) = [-1 -1 1];
G(8,:) = [-1 1 -1];
G(9,:) = [-1 -1 -1];

magicvecs = zeros(8,3);
magicvecs(1,:) = [1 1 1];
magicvecs(2,:) = [1 1 -1];
magicvecs(3,:) = [1 -1 1];
magicvecs(4,:) = [-1 1 1];
magicvecs(5,:) = [1 -1 -1];
magicvecs(6,:) = [-1 -1 1];
magicvecs(7,:) = [-1 1 -1];
magicvecs(8,:) = [-1 -1 -1];

V = zeros(8,8);

for gi = 1:8
    for gj = 1:gi
        tmpv = G(gi,:) - G(gj,:);
        for mi = 1:8
            if (tmpv == magicvecs(mi,:))
                fprintf(1, '(%d,%d)\t[%d %d %d]\n', gi, gj, tmpv);
                V(gi,gj) = -0.125 * alpha^3;
                V(gj,gi) = V(gi,gj);
            end
        end
    end
end
