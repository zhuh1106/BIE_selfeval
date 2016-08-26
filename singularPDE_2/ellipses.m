function G = ellipses( a, b)
% ELLIPSE - set up [0,2pi) parametrization of smooth elliptic closed curve

G.Z = @(t) (a+b)/2*exp(1i*t) + (a-b)/2*exp(-1i*t);
G.Zp = @(t) 1i*(a+b)/2*exp(1i*t) - 1i*(a-b)/2*exp(-1i*t);
G.Zpp = @(t) -(a+b)/2*exp(1i*t) - (a-b)/2*exp(-1i*t);

end

