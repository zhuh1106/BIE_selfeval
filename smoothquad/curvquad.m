function G = curvquad( G, rule, N)
% CURVQUAD - set up undelying quadrature for a closed curve struct

if strcmp(rule,'ptr')
    G.s = 2*pi*(1:N)'/N;    % column vector of points for parameterized curve
    G.w = 2*pi/N*ones(N,1);
end

s = G.s;
G.x = G.Z(s);
G.sp = abs(G.Zp(s));
G.nx = -1i*G.Zp(s)./G.sp;   % normal vector, requires dot operation
G.cur = -real(conj(G.Zpp(s)).*G.nx)./G.sp.^2;  % curvature, requires dot operation

end

