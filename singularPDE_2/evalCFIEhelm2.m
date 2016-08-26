function u = evalCFIEhelm2(t,G,sigma,k,eta)
% evaluate 2D Helmholtz CFIE (D-i.eta.S) potential at targets
% same function as eval CFIEhelm, instead helm2 uses matrix computation
% u = evalCFIEhelm(t,G,sigma,k,eta) where t is a list of target points
%  (points in the complex plane), G is a curve struct (see curvquad.m) with
%  N quadrature nodes, sigma is a N-element column vector, k is the wavenumber,
%  and eta (~k) controls the amount of SLP in the CFIE, returns the potential
%  at the target points evaluated using the quadrature rule in G.
%  A direct sum is used, vectorized over target points only.
%
% Hai 05/26/16

N = numel(G.x);
[ m, n] = size(t);

% reshape t, G.x, G.nx to perform matrix operation
tt = repmat( reshape(t,m*n,1), N, 1);   % column vector
GxGx = reshape( repmat( (conj(G.x))',m*n,1), m*n*N, 1);  % column vector, with first m*n entries equal G.x(1)
GnxGnx = reshape( repmat( (conj(G.nx))',m*n,1), m*n*N, 1);
GwGw = reshape( repmat( G.w',m*n,1), m*n*N, 1);
GspGsp = reshape( repmat( G.sp',m*n,1), m*n*N, 1);
sigmas = reshape( repmat( (conj(sigma))',m*n,1), m*n*N, 1);

% vector computation
d = tt - GxGx;
kr = k*abs(d);
costheta = real(conj(GnxGnx).*d)./abs(d);
uu = (k*costheta.*besselh(1,kr) - 1i*eta*besselh(0,kr)) .* ...
      ((1i/4) * GwGw .* GspGsp .* sigmas);
  
u = sum( reshape( uu, m*n, N), 2);
u = reshape( u, m, n);
