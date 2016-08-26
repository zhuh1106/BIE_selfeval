function a = CFIEnystAL(G,i,j,k,eta,ord)
% compute element of 2D Helmholtz CFIE Nystrom matrix with Aplert correction
%
% Aij = CFIEnystAL(G,i,j,k,eta) returns A_{ij} given i,j, and a curve struct G
% whose boundary quadrature must be PTR, wavenumber k, and mixing
% parameter.
% Hai 05/20/16.

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
    
    a = 0;  % initilize A_ij
    
    % distance between i and j source nodes mode N
    Dist = [ j-i, j-i-N, j-i+N];
    [ ~, idx] = min(abs(Dist));
    dist = Dist(idx);
    
    % add the " contribution of jth interpolation node to each auxi node"
    % to get coeff of A_ij
    for p=1:numel(tex)  
        
        % auxiliary node and corresponding weight
        kappa  = tex(p);    
        wkappa = wex(p);    
        
        % Lagrange interpolation coefficients, plus and minus kappa
        Lpjp = barycoefs( boff, kappa, dist);
        Lpjm = barycoefs( boff, -kappa, dist);
        
        % parametrized coordinate of auxiliary node in [0,2pi), plus and
        % minus kappa
        skappap = 2*pi*mod(i+kappa,N)/N;  
        skappam = 2*pi*mod(i-kappa,N)/N; 
        
        % evaluation of A_ij, including both plus and minus kappa
        Gxkappa = G.Z(skappap);  
        d = G.x(i)-Gxkappa;
        kr = k*abs(d);           
        Gspkappa = abs(G.Zp(skappap));   
        Gnxkappa = -1i*G.Zp(skappap)/Gspkappa;   
        costheta = real(conj(Gnxkappa).*d)/abs(d); 
        a = a + Lpjp * wkappa * G.w(j) * Gspkappa * (1i/4) * ...
            (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)); 
        
        Gxkappa = G.Z(skappam);  
        d = G.x(i)-Gxkappa;
        kr = k*abs(d);
        Gspkappa = abs(G.Zp(skappam));   
        Gnxkappa = -1i*G.Zp(skappam)/Gspkappa;   
        costheta = real(conj(Gnxkappa).*d)/abs(d);  
        a = a + Lpjm * wkappa * G.w(j) * Gspkappa * (1i/4) * ...
            (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)); 
    end
    
else 
        % j is far from i
        d = G.x(i)-G.x(j); kr = k*abs(d);           
        costheta = real(conj(G.nx(j)).*d)./abs(d);  
        a = (1i/4) * (k*costheta*besselh(1,kr) - 1i*eta*besselh(0,kr)) * sw;
end  

