function tree = getTest()

syms x y z;
syms x2 z3 z6;
tree(1).name = 'x2';
tree(1).expr = x^2;
tree(1).depend = {};
tree(2).name  = 'z3';
tree(2).expr  = z^3;
tree(2).depend  = {};
tree(3).name  = 'z6';
tree(3).expr  = z3^2;
tree(3).depend  = {'z3'};
tree(4).name  = 'f';
tree(4).expr = x2*z3*z6;
tree(4).depend = {'x2', 'z3', 'z6'};
