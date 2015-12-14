function [gamma, theta, phi, B, phi0, B0]=IRTM(X, y, terms)
rng(10725);
D = size(X,1);
V = size(X,2);
K = 10;
T = 100;
delta = zeros(V,1);
%deltab = 0;
alpha = 0.1;
eta = 0.1;
lambda = 1;

B = drchrnd(eta*ones(1,V), K);
B0 = B;

gamma = zeros(K,V);
%theta = rand(D,K);
%theta_s = sum(theta,2);
%ds = theta_s*ones(1,K);
%theta = theta./ds;
theta = drchrnd(alpha*ones(1,K),D);
phi = rand(V,1);
phi0 = phi;


S = zeros(K,V);
for t = 1:T
    g = zeros(V,1);
    h = zeros(V,1);
    gb = zeros(K,V);
    hb = zeros(K,V);
   % fa = 0;
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
    g = g + y(d)*sum(gamma.*(ones(K,1)*X(d,:))-B_yd.*(gamma*X(d,:)'*ones(1,V)),1)';
    
    s_w = Ck - B.*exp(ones(K,1)*phi'*y(d));
    fd = s_w;
    delta = 0.1*(abs(phi));
    E1 = ones(K,1)*exp(phi'*y(d)+delta'*abs(y(d)));
    E2 = ones(K,1)*exp(phi'*y(d)-delta'*abs(y(d)));
    I = (E1<s_w);
    fd(I) = E1(I);
    I = (E2>s_w);
    fd(I) = E2(I);
    Fd = 2 + fd./s_w + s_w./fd;
    h = h+ (y(d)^2)*sum(((gamma*X(d,:)')*ones(1,V))./Fd, 1)';
    gb = gb + (ones(K,1)*X(d,:)).*(gamma.*(1-B_yd));
    hb = hb + (ones(K,1)*X(d,:)).*gamma;
    
   % fa = fa+(alpha-1)*sum(log(theta(d,:)))+ sum(sum((ones(K,1)*X(d,:)).*gamma.*log(B_yd)));
   
end
    
    H_inv = diag(-1./h);
    s = H_inv*g;
    sl = H_inv*ones(V,1);
    p1 = s-lambda*sl;
    p2 = s+lambda*sl;
    phi_star = (phi-p1).*(phi>p1)+(phi-p2).*(phi<p2);
    phi = projectPhi(phi_star, phi, delta);
    
    gb = (eta-1+gb)./B;
    deltab = 0.1*abs(B);
    hb = -(eta-1+hb)./((B-deltab).^2);
   
    
    for i = 1:K
        H_inv = diag((1./hb(i,:))');
        s = H_inv*gb(i,:)';
        S(i,:) = s';
    end
    
    B_star = B-S;
    B1 = B-deltab;
    B2 = B+deltab;
    B = B_star;
     I = (B_star<B1);
    B(I) = B1(I);
    I = (B_star>B2);
     B(I) = B2(I);
    I = (B<=0);
    B(I) = 10^(-9);
    BN = sum(B,2);
    BN = BN * ones(1,V);
    B = B./BN;
    
    
    fv = ojbectValue(X, y, gamma, theta, phi, B, alpha, eta, lambda, K);
    fprintf('%f\n', fv);
    %showWord(B, terms);
end




end
