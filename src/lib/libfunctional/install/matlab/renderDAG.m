function bg = renderDAG(DAG)

ID = {DAG(:).name};
conn = zeros(length(DAG));

for k = 1:length(DAG)
    depend = DAG(k).depend;
    for l = 1:length(depend)
        index = find(strcmp(ID,depend{l}));
        conn(index, k) = 1;
    end
end

bg = biograph(conn,ID);
set(get(bg,'Edges'), 'LineColor', [0 0 0]);
set(get(bg,'Nodes'), 'Color', [1 1 1]);
set(bg,'LayoutType','radial');
view(bg)
