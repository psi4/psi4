function code = UKSCodeRev(trees, trees_a0, trees_b0) 

prims = trees.prims;
inter0 = trees.inter0;
inter1 = trees.inter1;
inter2 = trees.inter2;
fin0 = trees.fin0;
fin1 = trees.fin1;
fin2 = trees.fin2;
fun = [inter0 fin0];

newline = sprintf('\n');

% HEADER
code = ['/** Reverse accumulation AD code for functional value. **/' newline];
code = [code '/**' newline];
code = [code printTree(fun)];
code = [code '**/' newline];
code = [code newline];

% rho_a > TOL && rho_b > TOL
code = [code 'if (rho_a > cutoff_ && rho_b > cutoff_) { //Open Cutoff' newline];

for k = 1:length(inter0)
    ccode(inter0(k).expr,'file',[inter0(k).name '_0']);
    code = [code '// Zero-th Order Intermediate ' inter0(k).name ':' newline]; 
    code = [code cleanTree([inter0(k).name '_0'],inter0(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin0)
    ccode(fin0(k).expr,'file',[fin0(k).name '_0']);
    code = [code '// Zero-th Order Final Value ' fin0(k).name ':' newline]; 
    code = [code cleanTree([fin0(k).name '_0'],fin0(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 0) { // Open Gradient' newline];

for k = 1:length(inter1)
    ccode(inter1(k).expr,'file',[inter1(k).name '_1']);
    code = [code '// First Order Intermediate ' inter1(k).name ':' newline]; 
    code = [code cleanTree([inter1(k).name '_1'],inter1(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin1)
    ccode(fin1(k).expr,'file',[fin1(k).name '_1']);
    code = [code '// First Order Final Value ' fin1(k).name ':' newline]; 
    code = [code cleanTree([fin1(k).name '_1'],fin1(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 1) { // Open Hessian' newline];

for k = 1:length(inter2)
    ccode(inter2(k).expr,'file',[inter2(k).name '_2']);
    code = [code '// Second Order Intermediate ' inter2(k).name ':' newline]; 
    code = [code cleanTree([inter2(k).name '_2'],inter2(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin2)
    ccode(fin2(k).expr,'file',[fin2(k).name '_2']);
    code = [code '// Second Order Final Value ' fin2(k).name ':' newline]; 
    code = [code cleanTree([fin2(k).name '_2'],fin2(k).name, '=', '')];
    code = [code newline];
end

code = [code '}} // Close Gradient, Close Hessian' newline];

% rho_a < TOL && rho_b > TOL (USE a0)
code = [code '} else if (rho_b > cutoff_) {' newline];

prims = trees_a0.prims;
inter0 = trees_a0.inter0;
inter1 = trees_a0.inter1;
inter2 = trees_a0.inter2;
fin0 = trees_a0.fin0;
fin1 = trees_a0.fin1;
fin2 = trees_a0.fin2;

for k = 1:length(inter0)
    ccode(inter0(k).expr,'file',[inter0(k).name '_0']);
    code = [code '// Zero-th Order Intermediate ' inter0(k).name ':' newline]; 
    code = [code cleanTree([inter0(k).name '_0'],inter0(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin0)
    ccode(fin0(k).expr,'file',[fin0(k).name '_0']);
    code = [code '// Zero-th Order Final Value ' fin0(k).name ':' newline]; 
    code = [code cleanTree([fin0(k).name '_0'],fin0(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 0) { // Open Gradient' newline];

for k = 1:length(inter1)
    ccode(inter1(k).expr,'file',[inter1(k).name '_1']);
    code = [code '// First Order Intermediate ' inter1(k).name ':' newline]; 
    code = [code cleanTree([inter1(k).name '_1'],inter1(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin1)
    ccode(fin1(k).expr,'file',[fin1(k).name '_1']);
    code = [code '// First Order Final Value ' fin1(k).name ':' newline]; 
    code = [code cleanTree([fin1(k).name '_1'],fin1(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 1) { // Open Hessian' newline];

for k = 1:length(inter2)
    ccode(inter2(k).expr,'file',[inter2(k).name '_2']);
    code = [code '// Second Order Intermediate ' inter2(k).name ':' newline]; 
    code = [code cleanTree([inter2(k).name '_2'],inter2(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin2)
    ccode(fin2(k).expr,'file',[fin2(k).name '_2']);
    code = [code '// Second Order Final Value ' fin2(k).name ':' newline]; 
    code = [code cleanTree([fin2(k).name '_2'],fin2(k).name, '=', '')];
    code = [code newline];
end

code = [code '}} // Close Gradient, Close Hessian' newline];

% rho_a > TOL && rho_b < TOL (USE b0)
code = [code '} else if (rho_a > cutoff_) {' newline];

prims = trees_b0.prims;
inter0 = trees_b0.inter0;
inter1 = trees_b0.inter1;
inter2 = trees_b0.inter2;
fin0 = trees_b0.fin0;
fin1 = trees_b0.fin1;
fin2 = trees_b0.fin2;

for k = 1:length(inter0)
    ccode(inter0(k).expr,'file',[inter0(k).name '_0']);
    code = [code '// Zero-th Order Intermediate ' inter0(k).name ':' newline]; 
    code = [code cleanTree([inter0(k).name '_0'],inter0(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin0)
    ccode(fin0(k).expr,'file',[fin0(k).name '_0']);
    code = [code '// Zero-th Order Final Value ' fin0(k).name ':' newline]; 
    code = [code cleanTree([fin0(k).name '_0'],fin0(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 0) { // Open Gradient' newline];

for k = 1:length(inter1)
    ccode(inter1(k).expr,'file',[inter1(k).name '_1']);
    code = [code '// First Order Intermediate ' inter1(k).name ':' newline]; 
    code = [code cleanTree([inter1(k).name '_1'],inter1(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin1)
    ccode(fin1(k).expr,'file',[fin1(k).name '_1']);
    code = [code '// First Order Final Value ' fin1(k).name ':' newline]; 
    code = [code cleanTree([fin1(k).name '_1'],fin1(k).name, '=', '')];
    code = [code newline];
end

code = [code 'if (deriv_ > 1) { // Open Hessian' newline];

for k = 1:length(inter2)
    ccode(inter2(k).expr,'file',[inter2(k).name '_2']);
    code = [code '// Second Order Intermediate ' inter2(k).name ':' newline]; 
    code = [code cleanTree([inter2(k).name '_2'],inter2(k).name, '=', 'double')];
    code = [code newline];
end

for k = 1:length(fin2)
    ccode(fin2(k).expr,'file',[fin2(k).name '_2']);
    code = [code '// Second Order Final Value ' fin2(k).name ':' newline]; 
    code = [code cleanTree([fin2(k).name '_2'],fin2(k).name, '=', '')];
    code = [code newline];
end

code = [code '}} // Close Gradient, Close Hessian' newline];

% rho_a < TOL && rho_b < TOL
code = [code '} // Close Cutoff' newline];

