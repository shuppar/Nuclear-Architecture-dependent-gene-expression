function Ind = FokInd(Image_Stack)


%%
A = Image_Stack;
[~, ~, zs] = size(A);
ind = zeros(zs,1);

        % LAPV from fmeasure looks like the following:
        % LAP = fspecial('laplacian');
        % ILAP = imfilter(Image, LAP, 'replicate', 'conv');
        % FM = std2(ILAP)^2;

for i = 1:zs
   
    
    ind(i) = fmeasure(A(:,:,i), 'LAPV');
%     ind(i) = fmeasure(A(:,:,i), 'LAPE');
%     ind(i) = fmeasure(A(:,:,i), 'LAPD');
%     ind(i) = fmeasure(A(:,:,i), 'TENG');
%     ind(i) = fmeasure(A(:,:,i), 'TENV');
    
end

[fm, Ind] = max(ind);
%%

end