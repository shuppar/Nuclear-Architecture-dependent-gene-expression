%% To fit linear curves


A = load('Count_compare.dat');
x = A(:,3); y = A(:,4);
figure, scatter(A(:,3),A(:,4),40,'MarkerEdgeColor',[0 .4 .6], 'MarkerFaceColor',[0 .6 .8], 'LineWidth',1.5);
hold on;
plot(f);
[f, g] = fit(x, y, 'poly1');

%% For Bar Graphs comparing Cell Cycle Distribution in Control and Damage Cells.

 U = load('G1.dat');
 V = load('S.dat');
 W = load('G2_M.dat');
 cd ../../
 cd Control/CycleData/
 X = load('G1.dat');
 Y = load('S.dat');
 Z = load('G2_M.dat');
 a = mean(U); b= mean(V); c= mean(W);
 d = std(U)/sqrt(length(U)); e = std(V)/sqrt(length(V)); f = std(W)/sqrt(length(W));
 g = mean(X); h = mean(Y); i = mean(Z);
 j = std(X)/sqrt(length(X)); k = std(Y)/sqrt(length(Y)); l = std(Z)/sqrt(length(Z));
 p = [a(:,4) g(:,4); b(:,4) h(:,4); c(:,4) i(:,4)];
 q = [d(:,4) j(:,4); e(:,4) k(:,4); f(:,4) l(:,4)];
 barwitherr(q, p);
