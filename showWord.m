function showWord(B, terms)

[~,I] = sort(B,2,'descend');
K = size(B,1);
for i = 1:K
C = terms(I(i,1:10));
display(strjoin(C', ', '));

end



end
