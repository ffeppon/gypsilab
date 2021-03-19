close all;
clear all;

%% Geometry
% Create a non radially symetric mesh mesh whose boundary has two connected components:
N = 50; rad1 = 1; rad2 = 2;
% m_Omega = mshAssymetricRing(N,rad1,rad2);
% m_Gamma = m_Omega.bnd;
% c = m_Gamma.ctr;
% indGamma1 = c(:,1).^2 + c(:,2).^2 > (rad1 + rad2)^2/4;
% indGamma2 = ~indGamma1;
% m_Gamma1 = m_Gamma.sub(indGamma1);
% m_Gamma2 = m_Gamma.sub(indGamma2);

N = 800;
m_Omega1 = mshDisk(N,0.1); 
m_Omega2 = translate(mshDisk(N,0.3),[-0.5,0,0]);
m_Omega = union(m_Omega1,m_Omega2);
m_Gamma = m_Omega.bnd;
[~,ind_Gamma1] = intersect(m_Gamma,bnd(m_Omega1));
[~,ind_Gamma2] = intersect(m_Gamma,bnd(m_Omega2));
m_Gamma1 = m_Gamma.sub(ind_Gamma1);
m_Gamma2 = m_Gamma.sub(ind_Gamma2);

m_Space = mshDisk(3*N,2);

% Domains of integration
gss = 3; % Number of Gauss points
Gamma = dom(m_Gamma,gss);
Gamma2 = dom(m_Gamma2,gss);

Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',0); % 0 wave number
S1_Gamma = fem(m_Gamma,'P1');
S1_Gamma2 = fem(m_Gamma2,'P1');

% Single layer potential
V = -1/(2*pi)*integral(Gamma,Gamma,S1_Gamma,Gxy,S1_Gamma);
V = V + -1/(2*pi)*regularize(Gamma,Gamma,S1_Gamma,'[log(r)]',S1_Gamma);

B = integral(Gamma,S1_Gamma); % Array containing the mean value of each basis function 
sys = [V B;B' 0]; % Solve with zero mean value constraing

% right hand side:
P = restriction(S1_Gamma,m_Gamma2); % Operator restrictinng dofs from S^1,0(Gamma) to S^{1,0}(Gamma2). P' extends by 0 from Gamma2 to Gamma. 
g = P'*ones(size(P,1),1); % Function equal to 1 on Gamma2 and 0 on Gamma1
F = integral(Gamma,S1_Gamma,S1_Gamma)*g;
rhs = [F;0];
sol = sys\rhs;
lambda = sol(1:end-1); % Normal trace
mu = sol(end); % Lagrange multiplier. 
fprintf('Lagrange Multiplier: %s\n',num2str(mu))
dnu = P*lambda; % Restriction to Gamma2
figure;
plot(m_Omega,'w');
hold on
surf(S1_Gamma,lambda);

axis equal;
title('Charge distribution on the conductor Omega2')

SL = (-1/(2*pi))*integral(m_Space.vtx,Gamma,Gxy,S1_Gamma);
% Array of int_{Gamma}G(x - y)phi_h(y) for x = x1 ... xN
SL = SL + (-1/(2*pi))*regularize(m_Space.vtx,Gamma,'[log(r)]',S1_Gamma);


u = SL*lambda+mu; 
S1_Omega = fem(m_Space,'P1');
figure
graph(S1_Omega,u);
title('Electrostatic potential')


I = integral(Gamma2,S1_Gamma2,ntimes(S1_Gamma2));
F = 1/2*dnu'*I{1}*dnu;
 
fprintf('Force (in the x direction) %s\n',num2str(F))