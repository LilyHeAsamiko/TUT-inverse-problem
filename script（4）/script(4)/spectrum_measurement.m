clear

sigma = 0.05;
n_iter = 4;
theta0 = 0.000125;
beta = 1.5;
eta = beta - 1.5;
kappa = beta + 1.5;
hypermodel = 'InverseGamma';
%Concentration function
d_t = 1e-3;
t = [-0.5 : d_t: 1.5];
g = exp(-((t-0.7)/0.2).^2) + exp(-((t-0.3)/0.1).^2);

figure(1); plot(t, g)
hold on

%A = randn(size(g));
A=2;
y = A.*g+sigma*randn(size(g));
figure(2); clf; plot(t,y);


theta = theta0*ones(size(A,2),1);

 for k = 1 : n_iter
  k
  
  d = sqrt(theta);  
  A_d = A*diag(d);

  x= d.*(A_d'*(( A_d * A_d' + (sigma^2)*eye(size(A_d,1),size(A_d,1))).\y));
  
  if strcmp(hypermodel,'Gamma')
      
  theta = 0.5*theta0*(eta + sqrt(eta^2 + 2*x.^2/theta0));
  theta = theta(1);
  
  elseif strcmp(hypermodel,'InverseGamma')
  
  theta = (theta0+0.5*x.^2)./kappa
  theta = theta(1);
  
  end
  
  figure(1); plot(t,x(1,:))
  
  end
  
