function a = CFIEnystAL2(G,i,j,k,eta,ord)
% compute element of 2D Helmholtz CFIE Nystrom matrix with Aplert correction
% same function as CFIEnystAL, instead AL2 uses vector computation for each
% A_ij
% Aij = CFIEnystAL(G,i,j,k,eta) returns A_{ij} given i,j, and a curve struct G
% whose boundary quadrature must be PTR, wavenumber k, and mixing
% parameter.
% Hai 05/26/16.

% speed, weight, problem dimension, and distance
sw = G.sp(j)*G.w(j); 
N = numel(G.x); 
l = mod(i-j,N); 
if l>N/2
    l=N-l; 
end   

% get auxiliary nodes for different order
[tex,wex,nskip] = QuadLogExtraPtNodes(ord);

% corresponding offset nodes that will be used to do interpolation
boff = -(nskip-1):(nskip-1);    

% check if  A_ij is within correction band of A_ii
if l>=0 && l<nskip
    
    a = 0*tex;  % initilize A_ij
    
    % distance between i and j source nodes mode N
    Dist = [ j-i, j-i-N, j-i+N];
    [ ~, idx] = min(abs(Dist));
    dist = Dist(idx);
    
    % add the " contribution of jth interpolation node to each auxi node"
    % to get coeff of A_ij
    %
    % auxiliary nodes, rhs of ith node
    Lpjp = barycoefs( boff,  tex, dist);    % column vector stores Lagrange interpolation coeff
    s = 2*pi*mod(i+tex,N)/N;    % auxi node coordinates in [0,2*pi)
    Gx = G.Z(s);    % auxi node coordinates in complex plain
    d = G.x(i)-Gx;  
    kr = k*abs(d);  % for CFIE kernel         
    Gsp = abs(G.Zp(s)); % speeds at auxi nodes
    Gnx = -1i*G.Zp(s)./Gsp; % normal vectors at auxi nodes
    costheta = real(conj(Gnx).*d)./abs(d);  % cos between normal vector and distance
    a = a + Lpjp .* (wex * G.w(j)) .* Gsp * (1i/4) .* ...
            (k*costheta.*besselh(1,kr) - 1i*eta*besselh(0,kr)); 
    
    % auxiliary nodes, lhs of ith node
    Lpjm = barycoefs( boff, -tex, dist);
    s = 2*pi*mod(i-tex,N)/N; 
    Gx = G.Z(s);  
    d = G.x(i)-Gx;
    kr = k*abs(d);           
    Gsp = abs(G.Zp(s));   
    Gnx = -1i*G.Zp(s)./Gsp;   
    costheta = real(conj(Gnx).*d)./abs(d); 
    a = a + Lpjm .* (wex * G.w(j)) .* Gsp * (1i/4) .* ...
            (k*costheta.*besselh(1,kr) - 1i*eta*besselh(0,kr)); 
    
    % sum over all 2*numel(tex) auxi nodes
    a = sum(a,1);
    
else 
        % j is far from i
        d = G.x(i)-G.x(j); kr = k*abs(d);           
        costheta = real(conj(G.nx(j)).*d)./abs(d);  
        a = (1i/4) * (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)) * sw;
end  

