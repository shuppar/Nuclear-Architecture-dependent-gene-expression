%% This function finds two points on a given closed object (S) such that those points are 
%% closest to a given point A(x, y) outside or inside the closed object and that those two 
%% points share either coordinates of the given point A(p, q).
%% That is set of all (x, y) such that (x, y) falls in to the S and
%% either x = p or y = q.

%% It is called crosshair because it looks like that if you imagine S as a cricle
%% and some point A(p, q) inside the circle. Then by constructing a crosshair passing 
%% through the point A(p, q) will give you the two points of interest.


function res = crosshair(tuple, Set)


%%
res = zeros(round(length(Set)/3), 2);
A = tuple(1); B = tuple(2);
X_set = Set(:,1); Y_set = Set(:,2);

for i = 1:length(Set)
   
    if X_set(i) == A || Y_set(i) == B
        res(i,:) = [X_set(i), Y_set(i)]; % res for result.
    else
    end
    
    % Keeping just the non-zero elements (which are indices op points).
    res(res(:,1)==0, :) = []; % because in Matlab index starts with 1.
    
end

%%
d = zeros(length(res), 1); sl = d;

for i = 1:length(res)
   
    tx = res(i, 1); ty = res(i, 2);
    d(i) = sqrt((A-tx)^2 + (B-ty)^2);
    clear tx; clear ty;
    
end

bada = sort(d, 'ascend');
chota = bada(1); bada = bada(2);

for i = 1:length(res)
   
    if d(i) == chota || d(i) == bada
        sl(i) = i;
    else
    end
    sl(sl==0) = [];
end

%%

for i = 1:length(res)
    if ismember(i, sl)
    else
        res(i,:) = [0, 0];
    end
end
    
    res(res(:,1)==0, :) = [];

end