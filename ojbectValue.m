function fv = ojbectValue(X, y, gamma, theta, phi, B, alpha, eta, lambda, K)

D = size(X,1);
V = size(X,2);
fa = 0;
for d = 1:D
    Ck = B*exp(phi*y(d));
    Ck = Ck*ones(1,V);
    B_yd = B.*exp(ones(K,1)*y(d)*phi')./Ck;
    fa = fa+(alpha-1)*sum(log(theta(d,:)))+ sum(sum((ones(K,1)*X(d,:)).*gamma.*log(B_yd)));
end

fv = (eta - 1)*sum(sum(log(B)))-lambda*sum(abs(phi))+ fa;



end
