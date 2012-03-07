function buildDashD()

current_dir = dir('./');
if (any(strcmpi(['dispersion.cc'],{current_dir.name})))
    rm_str = sprintf('rm dispersion.cc');
    system(rm_str);
end


system(sprintf('cp dispersion.cc.template dispersion.cc'));

%-D1 and -D2
syms C6 d_ RvdW real
syms xi xj yi yj zi zj real

r = sqrt((xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2);
D1 = C6/r^6/(1+exp(-d_*(r/RvdW - 1)));
D1_x = diff(D1,xi);
D1_y = diff(D1,yi);
D1_z = diff(D1,zi);
D1_xA_xA = diff(diff(D1,xi),xi);
D1_xA_yA = diff(diff(D1,xi),yi);
D1_xA_zA = diff(diff(D1,xi),zi);
D1_yA_yA = diff(diff(D1,yi),yi);
D1_yA_zA = diff(diff(D1,yi),zi);
D1_zA_zA = diff(diff(D1,zi),zi);
D1_xA_xB = diff(diff(D1,xi),xj);
D1_xA_yB = diff(diff(D1,xi),yj);
D1_xA_zB = diff(diff(D1,xi),zj);
D1_yA_xB = diff(diff(D1,yi),xj);
D1_yA_yB = diff(diff(D1,yi),yj);
D1_yA_zB = diff(diff(D1,yi),zj);
D1_zA_xB = diff(diff(D1,zi),xj);
D1_zA_yB = diff(diff(D1,zi),yj);
D1_zA_zB = diff(diff(D1,zi),zj);

ccode(D1 + 0,'file','D1');
ccode(D1_x + 0,'file','D1_x');
ccode(D1_y + 0,'file','D1_y');
ccode(D1_z + 0,'file','D1_z');
ccode(D1_xA_xA + 0,'file','D1_xA_xA');
ccode(D1_xA_yA + 0,'file','D1_xA_yA');
ccode(D1_xA_zA + 0,'file','D1_xA_zA');
ccode(D1_yA_yA + 0,'file','D1_yA_yA');
ccode(D1_yA_zA + 0,'file','D1_yA_zA');
ccode(D1_zA_zA + 0,'file','D1_zA_zA');
ccode(D1_xA_xB + 0,'file','D1_xA_xB');
ccode(D1_xA_yB + 0,'file','D1_xA_yB');
ccode(D1_xA_zB + 0,'file','D1_xA_zB');
ccode(D1_yA_xB + 0,'file','D1_yA_xB');
ccode(D1_yA_yB + 0,'file','D1_yA_yB');
ccode(D1_yA_zB + 0,'file','D1_yA_zB');
ccode(D1_zA_xB + 0,'file','D1_zA_xB');
ccode(D1_zA_yB + 0,'file','D1_zA_yB');
ccode(D1_zA_zB + 0,'file','D1_zA_zB');

ccD1 = cleanDashD('D1','energy');
ccD1_x = cleanDashD('D1_x','grad[Ax]');
ccD1_y = cleanDashD('D1_y','grad[Ay]');
ccD1_z = cleanDashD('D1_z','grad[Az]');
ccD1_xA_xA = cleanDashD('D1_xA_xA','hess[Ax][Ax]');
ccD1_xA_yA = cleanDashD('D1_xA_yA','hess[Ax][Ay]');
ccD1_xA_zA = cleanDashD('D1_xA_zA','hess[Ax][Az]');
ccD1_yA_yA = cleanDashD('D1_yA_yA','hess[Ay][Ay]');
ccD1_yA_zA = cleanDashD('D1_yA_zA','hess[Ay][Az]');
ccD1_zA_zA = cleanDashD('D1_zA_zA','hess[Az][Az]');
ccD1_xA_xB = cleanDashD('D1_xA_xB','hess[Ax][Bx]');
ccD1_xA_yB = cleanDashD('D1_xA_yB','hess[Ax][By]');
ccD1_xA_zB = cleanDashD('D1_xA_zB','hess[Ax][Bz]');
ccD1_yA_xB = cleanDashD('D1_yA_xB','hess[Ay][Bx]');
ccD1_yA_yB = cleanDashD('D1_yA_yB','hess[Ay][By]');
ccD1_yA_zB = cleanDashD('D1_yA_zB','hess[Ay][Bz]');
ccD1_zA_xB = cleanDashD('D1_zA_xB','hess[Az][Bx]');
ccD1_zA_yB = cleanDashD('D1_zA_yB','hess[Az][By]');
ccD1_zA_zB = cleanDashD('D1_zA_zB','hess[Az][Bz]');

replaceInFile('ED1_0', ccD1, 'dispersion.cc');
replaceInFile('ED1_xC', ccD1_x, 'dispersion.cc');
replaceInFile('ED1_yC', ccD1_y, 'dispersion.cc');
replaceInFile('ED1_zC', ccD1_z, 'dispersion.cc');
replaceInFile('ED1_xA_xA', ccD1_xA_xA, 'dispersion.cc');
replaceInFile('ED1_xA_yA', ccD1_xA_yA, 'dispersion.cc');
replaceInFile('ED1_xA_zA', ccD1_xA_zA, 'dispersion.cc');
replaceInFile('ED1_yA_yA', ccD1_yA_yA, 'dispersion.cc');
replaceInFile('ED1_yA_zA', ccD1_yA_zA, 'dispersion.cc');
replaceInFile('ED1_zA_zA', ccD1_zA_zA, 'dispersion.cc');
replaceInFile('ED1_xA_xB', ccD1_xA_xB, 'dispersion.cc');
replaceInFile('ED1_xA_yB', ccD1_xA_yB, 'dispersion.cc');
replaceInFile('ED1_xA_zB', ccD1_xA_zB, 'dispersion.cc');
replaceInFile('ED1_yA_xB', ccD1_yA_xB, 'dispersion.cc');
replaceInFile('ED1_yA_yB', ccD1_yA_yB, 'dispersion.cc');
replaceInFile('ED1_yA_zB', ccD1_yA_zB, 'dispersion.cc');
replaceInFile('ED1_zA_xB', ccD1_zA_xB, 'dispersion.cc');
replaceInFile('ED1_zA_yB', ccD1_zA_yB, 'dispersion.cc');
replaceInFile('ED1_zA_zB', ccD1_zA_zB, 'dispersion.cc');

replaceInFile('ED2_0', ccD1, 'dispersion.cc');
replaceInFile('ED2_xC', ccD1_x, 'dispersion.cc');
replaceInFile('ED2_yC', ccD1_y, 'dispersion.cc');
replaceInFile('ED2_zC', ccD1_z, 'dispersion.cc');
replaceInFile('ED2_xA_xA', ccD1_xA_xA, 'dispersion.cc');
replaceInFile('ED2_xA_yA', ccD1_xA_yA, 'dispersion.cc');
replaceInFile('ED2_xA_zA', ccD1_xA_zA, 'dispersion.cc');
replaceInFile('ED2_yA_yA', ccD1_yA_yA, 'dispersion.cc');
replaceInFile('ED2_yA_zA', ccD1_yA_zA, 'dispersion.cc');
replaceInFile('ED2_zA_zA', ccD1_zA_zA, 'dispersion.cc');
replaceInFile('ED2_xA_xB', ccD1_xA_xB, 'dispersion.cc');
replaceInFile('ED2_xA_yB', ccD1_xA_yB, 'dispersion.cc');
replaceInFile('ED2_xA_zB', ccD1_xA_zB, 'dispersion.cc');
replaceInFile('ED2_yA_xB', ccD1_yA_xB, 'dispersion.cc');
replaceInFile('ED2_yA_yB', ccD1_yA_yB, 'dispersion.cc');
replaceInFile('ED2_yA_zB', ccD1_yA_zB, 'dispersion.cc');
replaceInFile('ED2_zA_xB', ccD1_zA_xB, 'dispersion.cc');
replaceInFile('ED2_zA_yB', ccD1_zA_yB, 'dispersion.cc');
replaceInFile('ED2_zA_zB', ccD1_zA_zB, 'dispersion.cc');

% -D3
syms s6_ s8_ C6 C8 sr6_ sr8_ a6 a8 RAB 
syms xi yi xj yj zi zj

r = sqrt((xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2);
D3_6 = -s6_*C6/r^6/(1+6*(r/(sr6_*RAB)^(-a6));
D3_8 = -s8_*C8/r^8/(1+6*(r/(sr8_*RAB)^(-a8));


