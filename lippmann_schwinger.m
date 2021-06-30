%% Loading gypsilab
clear
close all
prefix = '/home/florian/polybox/latex/draft_nystrom/numerics/gypsilab/martin/gypsilab/';
addpath([[prefix,'openMsh/']])
addpath([[prefix,'openDom/']])
addpath([[prefix,'openFem/']])
addpath([[prefix,'openHmx/']])
addpath([[prefix,'openFe/']])
addpath(genpath([[prefix,'openEbd/']]))
addpath([[prefix,'miscellaneous/']])

%% Output

%% Mesh of the disk

N = 1000;
mesh = mshDisk(N,1);
Omega = dom(mesh,3);

Vh = fem(mesh,'P1');

%% Fine mesh for visualisation
N=5000;
mesh2 = translate(mshDisk(N,4),[1,0,0]);
Omega2 = dom(mesh2,3);
Vh2 = fem(mesh2,'P1');

%% Integral equation
k=5;

theta=0;
% Incident wave
F = @(X) exp(1i*k*(cos(theta)*X(:,1)+sin(theta)*X(:,2)));

% Green Kernel
Gxy = @(X,Y) -1i/(4)*femGreenKernel(X,Y,'[H0(kr)]',k);

tol=1e-3;
lambda =1;
Id = integral(Omega, Vh,Vh);
A = integral(Omega,Omega,Vh,Gxy,Vh);
A = A+1/(2*pi)*regularize(Omega,Omega,Vh,'[log(r)]',Vh);

% A = -1i/4*integralEbd(Omega,Omega,Vh,'[H0(kr)]',k,Vh,tol,lambda);
% A = A+1/(2*pi)*regularize(Omega,Omega,Vh,'[log(r)]',Vh);

%%
coeff = 10;
K = Id - (1-coeff)*k^2*A;
b = integral(Omega,Vh,F);

x = K\b;
% x = gmres(@(x)(K*x),b,[],1e-6,100);
system('mkdir output')
%clc;
%axis equal;
%plot(mesh);

%figure;
%graph(Vh,real(x));
%colormap('jet')
%axis equal;

%% Visualisation

pts2 = mesh2.vtx;
Interp = integral(pts2,Omega,Gxy,Vh); % J'ai enlevé -1i/4 en facteur ici (il y est déjà dans Gxy)
Interp = Interp+1/(2*pi)*regularize(pts2,Omega,'[log(r)]',Vh);


xInterp = F(mesh2.vtx)+(1-coeff)*k^2*Interp*x;

figure
axis equal
axis tight
colormap('jet')
graph(Vh2,real(xInterp));

%% Saving

pts = mesh.vtx;
elems = mesh.elt;
elems2 = mesh2.elt;
save('output/lippmann.mat','x','xInterp','pts','pts2','elems','elems2')
