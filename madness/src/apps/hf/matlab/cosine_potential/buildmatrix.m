function A = buildlaplacian(n)
n = 15;

A = zeros(n^3, n^3);

for i = 1:n
    for j = 1:n
        for k = 1:n
        didx = sub2ind([n n n], i, j, k);
        A(didx,didx) = 6;
            %% do x
            if (i == n)
                iplus = 1;
            else
                iplus = i + 1;
            end
            if (i == 1)
                iminus = n;
            else
                iminus = i - 1;
            end
            ipidx = sub2ind([n n n], iplus, j, k);
            imidx = sub2ind([n n n], iminus, j, k);
            %% do y
            if (j == n)
                jplus = 1;
            else
                jplus = j + 1;
            end
            if (j == 1)
                jminus = n;
            else
                jminus = j - 1;
            end
            jpidx = sub2ind([n n n], i, jplus, k);
            jmidx = sub2ind([n n n], i, jminus, k);
            %% do z
            if (k == n)
                kplus = 1;
            else
                kplus = k + 1;
            end
            if (k == 1)
                kminus = n;
            else
                kminus = k - 1;
            end
            kpidx = sub2ind([n n n], i, j, kplus);
            kmidx = sub2ind([n n n], i, j, kminus);
            idx = sub2ind([n n n],i,j, k);
            A(idx, ipidx) = -1;
            A(idx, imidx) = -1;
            A(idx, jpidx) = -1;
            A(idx, jmidx) = -1;
            A(idx, kpidx) = -1;
            A(idx, kmidx) = -1;
        end
%        fprintf(1, '(%d,%d) (%d) %d %d %d %d\n', i, j, idx, ipidx, imidx, jpidx, jmidx);
    end
end
