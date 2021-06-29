m = mshDisk(1e4,1);
Gamma = dom(m,3);
Vh = fem(m,'P1');
green= '[H0(kr)]';
tol = 1e-3;
k = 10;
lambda = 1; % Paramètre de la méthode, tu peux jouer avec pour voir si ça améliore la perf. Cela n'affecte pas la précision.

S = integralEbd(Gamma,Gamma,Vh,green,k,Vh,tol,lambda);
S = S + -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[1/r]',Vh); % Changer cette ligne avec [log(r)] pour utiliser le code d'ignacio.
u = randn(size(S,2),1);
S*u;