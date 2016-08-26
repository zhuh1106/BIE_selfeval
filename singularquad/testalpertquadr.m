% test Alpert log-quadratures on periodic log singularity's Fourier coeffs
% Hai 05/20/16

clear
close
ords = [2,3,4,5,6,8,10,12,14,16,pi];

% set up singular integration kernel
P = 2*pi;    
m = 10;
f = @(t) log(4*sin(t/2).^2) .* cos(2*pi*m*t/P)/P;

% exact solution of integration
Iex = 0; 
if m~=0
    Iex = -1/abs(m); 
end

% number of node in evaluation of quadrature
Ns = 10:2:70;
IN = nan(numel(Ns), numel(ords));

% test for different order of Alpert
for j=1:numel(ords)
    ord = ords(j); 
    fprintf('\n ord = %d:\n', ord);
  
  % test for different number of nodes 
  for i=1:numel(Ns)
      N = Ns(i);
      h = P/N;
      [tex,wex,nskip] = QuadLogExtraPtNodes(ord); % extra nodes
      
      if N>2*nskip
        t = [tex; (nskip:(N-nskip))'; N-tex(end:-1:1)] * h; % nodes correction (easy because we only need to correct end points )
        w = [wex; ones(N-2*nskip+1,1); wex(end:-1:1)] * h; % weights
        IN(i,j) = sum(f(t).*w);
      end

    fprintf('N=%d: err = %g\n', N, abs(IN(i,j)-Iex))
  end
end
figure
semilogy(Ns, abs(IN-Iex),'+-')
xlabel('N'); ylabel('error'); title('Alpert log endpoint correction errors');
axis([min(Ns) max(Ns) 1e-16 1]);
