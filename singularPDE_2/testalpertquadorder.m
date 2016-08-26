% Test code for order of 2D singular quadrature with Alpert correction. 
%
% The signular kernel I use is 2D fundamental solution of Helmholtz
% One known solution is f = @(z) besselh(0,k*abs(z-(0.2+0.3i)))
% Hai 05/25/16

%%
clear
close all
a = 2; b = 1;
k = 10; eta = k; ord = pi;

% quadrature rule
G = ellipses(a,b);                 % curve shape

% set up target point to observe convergence
t = 2.5*a-1.5i*b;    

% set up number of source points
Ns = 10:10:150; 
e = 0*Ns;

%% loop over Ns, dimension of underlying nodes
for k=1:numel(Ns)
  
    N = Ns(k);
    G = curvquad(G, 'ptr', N);   % curve quadrature rule

    % get density function for prescribed solution
    f = @(z) besselh(0,k*abs(z-(0.2+1i*0.3)));   % known solution
    rhs = 2*f(G.x); % boundary data
    A = nan(N,N);   % initilize Nystrom matrix

    for i = 1:N
        for j = 1:N
            A(i,j) = 2*CFIEnystAL2(G,i,j,k,eta,ord); 
        end
    end
    sigma = (eye(N)+A)\rhs;
    
    % evaluation of kernel and error with density function sigma
    u = evalCFIEhelm(t,G,sigma,k,eta); 
    e(k) = u-f(t); 
end

%% plot of error
figure()
semilogy(Ns,abs(e),'+-'); 
title('smooth kernel quadrature error convergence');
hold on
% least square approximation for semilogy scale
lsq = ([Ns;ones(1,numel(Ns))]*[Ns;ones(1,numel(Ns))]')\([Ns;ones(1,numel(Ns))]*log10(abs(e))');
plot(Ns,10.^(lsq'*[Ns;ones(1,numel(Ns))]),'+-');
hold off

figure()
loglog(Ns,abs(e),'+-')
title('smooth kernel quadrature error convergence, loglog scale');
hold on
% least square approxiamtion for loglog scale
lsq = ([log10(Ns);ones(1,numel(Ns))]*[log10(Ns);ones(1,numel(Ns))]')\([log10(Ns);ones(1,numel(Ns))]*log10(abs(e))');
plot(Ns,10.^(lsq'*[log10(Ns);ones(1,numel(Ns))]),'+-');
hold off
