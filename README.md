# iRFDA
completely intrinsic framework for Riemannian functional data analysis

# Example

clear;

% define the 2D hyperbolic manifold

mfd = hyperbolic();

d = mfd.d;

% synthetic mean function and Riemannian slope function on the manifold

m = 100; % #observations per curve

t = linspace(0.005,0.995,m);

theta_t = 2*t.^3+t.^2+t-2;

varphi_t = t.^3-t.^2+pi*t+1;

mu = [sinh(theta_t).*cos(varphi_t);sinh(theta_t).*sin(varphi_t);cosh(theta_t)];

K = 20;

phi = mdfourier(m,K,diag(ones(1,d)),t);

beta_coef = 1.5*(1./(1:K).^2);

betaZ = zeros(d,m); % coefficients of beta in the frame

for k = 1:K

    betaZ = betaZ + beta_coef(k) * phi(:,:,k);
    
end

betaV = mfd.coef_to_log(mu,betaZ); % tangent-vector-valued beta

beta = mfd.Exp(mu,betaV); % beta on the manifold

%% generate Riemannian functional data

n = 100;

eigfunc = mdfourier(m,K,diag(ones(1,d)),t); % coefficients of eigenfunctions in the frame

lam = 1 ./ ((1:K).^2); % eigenvalues

xi = randn(n,K).* repmat(sqrt(lam(:)'),n,1); % RFPC scores

Z = zeros(d,m,n); % realizations of the coefficient process

for k = 1:K

    xik = xi(:,k);
    
    phik = eigfunc(:,:,k); % 
    
    for i = 1:n
    
        Z(:,:,i) = Z(:,:,i) + (xik(i)*phik);
        
    end
    
end

V = mfd.coef_to_log(mu,Z); % tangent-vector-valued functions

X = mfd.Exp(mu,V); % manifold-valued functions

%%

iRFPCA(X,mfd) % performance iRFPCA on X

