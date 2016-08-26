% test the smooth kernel quadrature. Hai 05/17/16

clear
close
a = 2; b = 1;
G = ellipses(a,b);                 % curve shape
N = 150; 
G = curvquad(G, 'ptr', N);   % curve quadrature rule

% set up target points
haxis = -a:.01:a; 
vaxis = -b:.01:b;
[ xx, yy] = meshgrid(haxis, vaxis); 
t = xx + 1i*yy; % targets

% density function
sigma = ones(N,1);  % with density function 1, the solution is supposed to be -1

% evaluation of smooth kernel with density function sigma
u = smoothquad(t,G,sigma); 

% plot of error
figure
imagesc(haxis,vaxis,log10(abs(u+1))); 
caxis([-16 0]); 
colorbar;
hold on; 

% plot of geometry of boundary
plot([G.x;G.x(1)], '-'); 
plot(G.x, 'k.','markersize',5); % show geom
axis xy equal tight; 
title('interior Laplace DLP log_{10} error');
set(gcf,'paperposition',[0 0 4 3]);
print -depsc2 ../testevalDLP.eps
