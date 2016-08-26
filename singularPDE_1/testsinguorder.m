% test the 2D Helmholtz evaluator, convergence at a point. 
% Hai 05/19/16

clear
close

% set up boundary
a = 2; b = 1;
k = 10; eta = k;
G = ellipses(a,b);                

% set up target point to observe convergence
t = 2.5*a-0.1*i*b;    

% set up number of source points
Ns = 10:10:150; 
e = 0*Ns;

% loop over different numbers of source point 
for k=1:numel(Ns)
  
    N = Ns(k);
    G = curvquad(G, 'ptr', N);   
    f = @(z) besselh(0,k*abs(z-(0.2+1i*0.3)));  % known solution
    rhs = 2*f(G.x); % boundary data
    A = nan(N,N);   % initilize Nystrom matrix
    for i = 1:N
        for j = 1:N
            A(i,j) = 2*CFIEnystKR(G,i,j,k,eta); 
        end
    end
    sigma = (eye(N)+A)\rhs;
  
    u = evalCFIEhelm(t,G,sigma,k,eta);
    e(k) = u-f(t); 
end

% plot
figure
semilogy(Ns,abs(e),'+-'); 
%axis tight; 
title('singular kernel quadrature error convergence');
hold on
lsq = ([Ns;ones(1,numel(Ns))]*[Ns;ones(1,numel(Ns))]')\([Ns;ones(1,numel(Ns))]*log10(abs(e))');
plot(Ns,10.^(lsq'*[Ns;ones(1,numel(Ns))]),'+-');

