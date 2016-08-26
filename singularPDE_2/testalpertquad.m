% Test code for 2D singular quadrature with Alpert correction. Hai 05/24/16
% The signular kernel I use is 2D fundamental solution of Helmholtz
% One known solution is f = @(z) besselh(0,k*abs(z-(0.2+0.3i)))

%%
clear
close all
a = 2; b = 1;
k = 10; eta = k; ord = 16;

% quadrature rule
G = ellipses(a,b);                 % curve shape
N = 100; 
G = curvquad(G, 'ptr', N);   % curve quadrature rule

%% get density function for prescribed solution
f = @(z) besselh(0,k*abs(z-(0.2+1i*0.3)));   % known solution
rhs = 2*f(G.x); % boundary data
A = nan(N,N);   % initilize Nystrom matrix

t1 = 0;
t2 = 0;
for i = 1:N
    for j = 1:N
        tic;
        A(i,j) = 2*CFIEnystAL(G,i,j,k,eta,ord); 
        time1 = toc;
        t1 = t1+time1;
        
        tic
        A(i,j) = 2*CFIEnystAL2(G,i,j,k,eta,ord);
        time2 = toc;
        t2 = t2+time2;
    end
end
t1,t2

sigma = (eye(N)+A)\rhs;

%% set up target points
haxis = -2*a:.05:2*a; 
vaxis = -2*b:.05:2*b;
[ xx, yy] = meshgrid(haxis, vaxis); 
t = xx + 1i*yy; % targets

%% evaluation of kernel with density function sigma
tic
u = evalCFIEhelm2(t,G,sigma,k,eta); 
toc

tic
u = evalCFIEhelm(t,G,sigma,k,eta);
toc
%max(max(abs(u-u2)))

% plot of error
figure
imagesc(haxis,vaxis,log10(abs(u-f(t)))); 
caxis([-16 0]); 
colorbar;
hold on 
plot([G.x;G.x(1)], '-'); 
plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; 
title('ext Dirichlet 2D Helm BVP log_{10} error')
%set(gcf,'paperposition',[0 0 4 3]); print -depsc2 ../helmbvperr.eps

figure
imagesc(haxis,vaxis,real(u)); 
caxis(.4*[-1 1]); colorbar;
hold on; plot([G.x;G.x(1)], '-'); 
plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; 
title('ext Dirichlet 2D Helm BVP Re u')
%set(gcf,'paperposition',[0 0 4 3]); print -depsc2 ../helmbvp.eps
