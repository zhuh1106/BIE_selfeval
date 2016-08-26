% test the DLP Laplace 2D evaluator evalDLP, convergence at a point. 
% Hai 05/17/16

clear
close

% set up boundary
a = 2; b = 1;
G = ellipses(a,b);                

% set up target point to observe convergence
t = 0.2*a-.3i*b;    

% set up number of source points
Ns = 10:10:200; 
e = 0*Ns;

% loop over different numbers of source point 
for i=1:numel(Ns)
  N = Ns(i);
  G = curvquad(G, 'ptr', N);   
  sigma = ones(N,1);  
  u = smoothquad(t,G,sigma);
  e(i) = u+1; 
end

% plot
figure
semilogy(Ns,abs(e),'+-'); 
title('smooth kernel quadrature error convergence');

