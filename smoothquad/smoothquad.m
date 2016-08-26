function u = smoothquad(t,G,sigma)
% EVALDLP - evaluate 2D integration of smooth kernel 
% for smooth kernel, just use the kernal function value at those points
% the difference for singluar kernel is that we need to cook up a matrix A

N = numel(G.x);
u = 0*t;        % initialize output
for j=1:N
  d = t-G.x(j);  
  u = u + (G.w(j) * G.sp(j) * sigma(j)) * real(conj(G.nx(j)).*d)./abs(d).^2;
end
u = u/(2*pi);
