function [L,R,Q] = CGMM_EM(X, init)
% Perform EM algorithm for fitting the Complex Gaussian mixture model.

% Written by Bo Ren (justnow.ren@gmail.com).
%% init
fprintf('EM for Complex Gaussian mixture: running ... \n');
tol = 1e-6;
maxiter = 50;

Q = -inf(1,maxiter);
[L,R]= initialization(X,init);
for iter = 2:maxiter
    %[~,label(1,:)] = max(L,[],2);
    %L = L(:,unique(label));   % remove empty clusters
    [phi, R] = maximization(X,L,R);
    [L,Q(iter)] = expectation(X, phi, R);
    if abs(Q(iter)-Q(iter-1)) < tol*abs(Q(iter)); break; end;
end
Q = Q(2:iter);

function [L, R] = initialization(X, init)
[n,m,f] = size(X);

k = init;
label = ceil(k*rand(1,n));
L = full(sparse(1:n,label,1,n,k,n));
L = repmat(L,1,1,f);
R = zeros(m,m,2,f);
for i = 1:f
    R(:,:,1,i) = X(:,:,i)'*X(:,:,i)/n;
    R(:,:,2,i) = eye(m)+eps;
end

function [L, Q] = expectation(X, phi, R)

[n,v,f] = size(phi);
Lp = zeros(n,v,f);
parfor k = 1:f
    for i = 1:v
        for j = 1:n
            Lp(j,i,k) = logCNpdf(X(j,:,k),phi(j,i,k)*R(:,:,i,k));
        end
    end
end 
T = logsumexp(Lp,2);
L = exp(bsxfun(@minus,Lp,T));
Q = sum(sum(sum(Lp.*L)));

function [phi, R] = maximization(X, L, R)
[n,m,f] = size(X);
v = size(L,2);
nl = permute(sum(L,1),[2,3,1]);
phi = zeros(n,v,f)+eps;

parfor k = 1:f
    for i = 1:v
        tR = zeros(m,m);
        for j = 1:n
            phi(j,i,k)=trace(X(j,:,k)'*X(j,:,k)/R(:,:,i,k))/m;
            tR = tR + phi(j,i,k)\(X(j,:,k)'*X(j,:,k))*L(j,i,k);
        end
        R(:,:,i,k) = tR/nl(i,k);
    end
end

%similar to guassian distribution
%p(x) = 1/(pi*Sigma)*exp[-(x-mu)^2/Sigma]
function y = logCNpdf(X,Sigma)
    c = log(pi)+log(det(Sigma));
    d = X/Sigma*X';
    y = real(- (c+d));


function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;

