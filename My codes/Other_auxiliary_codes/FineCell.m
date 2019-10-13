%% 

function out = FineCell(FileName, ColumnNumber, Resolution, G1_Peak_file)


%%
F = Resolution;
G1 = load(G1_Peak_file);
G2 = G1(1,2)/G1(1,1);
G1 = G1(1,1);
A = load(FileName);
A = sortrows(A, 1);
DNA = A(:,1);
DNA = DNA/G1;
G1 = G1(1,1)/G1(1,1);
B = A(:, ColumnNumber);
l = length(DNA);
n = l/F;
n = floor(n);


D = zeros(n,1);
P = D;

for i = 1:n
    
    r1 = (i-1)*F + 1;
    r2 = i*F;
    
    D(i) = mean(DNA(r1:r2));
    P(i) = mean(B(r1:r2));


end


f = fit(D, P, 'smoothingspline', 'SmoothingParam', 0.9);
plot(f, D, P);
hold on;
plot(D, P, 'o', 'LineWidth', 2,'MarkerSize', 5, 'MarkerEdgeColor', [0.3, 0.3, 0.6], 'MarkerFaceColor', [0.2, 0.2, 0.5]);
plot(G1, 0, 'o', 'MarkerSize', 13, 'MarkerEdgeColor', [0.5, 0.6, 0.4], 'MarkerFaceColor', [0.4, 0.3, 0.2]);
plot(G2, 0, 'o', 'MarkerSize', 13, 'MarkerEdgeColor', [0.5, 0.6, 0.4], 'MarkerFaceColor', [0.4, 0.3, 0.2]);
hold off;


