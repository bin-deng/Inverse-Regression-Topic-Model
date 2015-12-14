function [y, e] = predictY (X, y_test, phi, B, alpha, K)

rng(10725);
D = size(X,1);
V = size(X,2);
T = 10;
e = 0;
y = rand(D,1)-0.5;
theta = drchrnd(alpha*ones(1,K),D);

for t = 1:T
   
    fa = 0;
for d = 1:D
    Ck = B*exp(phi*y(d));
    Ck = Ck*ones(1,V);
    B_yd = B.*exp(ones(K,1)*y(d)*phi')./Ck;
    gamma = (theta(d,:)'*ones(1,V)).*B_yd;
    gamma_s = sum(gamma,1);
    gs = ones(K,1)*gamma_s;
    gamma = gamma./gs;
    theta(d,:) = X(d,:)*(gamma')+alpha-1;
    th = theta(d,:);
    I = (th<=0);
    th(I) = 10^(-6);
    theta(d,:) = th/sum(th);
    Bs = B*(exp(phi*y(d)).*phi)*ones(1,V);
    g = sum(sum((ones(K,1)*X(d,:)).*gamma.*(ones(K,1)*phi'-Bs./Ck)));
    h = -sum(sum((ones(K,1)*X(d,:)).*gamma.*(Ck.*(B*(exp(phi*y(d).*(phi.^2)))*ones(1,V)-Bs.^2)./(Ck.^2))));
   
    delta = 0.1*(abs(y(d)));
    y_star = y(d)-1/h*g;
    if (y_star < y(d)-delta)
        y(d) = y(d)-delta;
    elseif (y_star > y(d)+delta)
        y(d) = y(d) + delta;
    else
        y(d) = y_star;
    end
    Ck = B*exp(phi*y(d));
    Ck = Ck*ones(1,V);
    B_yd = B.*exp(ones(K,1)*y(d)*phi')./Ck;
    fa = fa+(alpha-1)*sum(log(theta(d,:)))+ sum(sum((ones(K,1)*X(d,:)).*gamma.*log(B_yd)));
   
end
    
    
     
    
    fv = fa;
    e = sum((y-y_test).^2);
    fprintf('%s\t%f\n', fv, e);
   
end







end
