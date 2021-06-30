%% Loading gypsilab
clear all
close all
addpathGypsilab()

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
tol = 1e-3;
theta=0;
% Incident wave
F = @(X) exp(1i*k*(cos(theta)*X(:,1)+sin(theta)*X(:,2)));

% Green Kernel
Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);


lambda =1;
Id = integral(Omega, Vh,Vh);
AR = regularize(Omega,Omega,Vh,'[log(r)]',Vh);

tic;
A1 = integral(Omega,Omega,Vh,Gxy,Vh,tol);
toc;

tic;
A2 = integralEbd(Omega,Omega,Vh,'[H0(kr)]',k,Vh,tol,lambda);
toc;

A1 = -1i/4*A1 + 1/(2*pi)*AR;
A2 = -1i/4*A2 + 1/(2*pi)*AR;

%%
coeff = 10;
K1 = Id - (1-coeff)*k^2*A1;
K2 = Id - (1-coeff)*k^2*A2;
b = integral(Omega,Vh,F);

tic;
x1 = K1\b;
toc;
tic;
x2 = gmres(@(x)(K2*x),b,[],1e-8,size(K2,1));
toc;


%% Visualisation

pts2 = mesh2.vtx;
tic;
Interp1 = integral(pts2,Omega,Gxy,Vh);
toc;

tic; 
Interp2 = integralEbd(pts2,Omega,'[H0(kr)]',k,Vh,tol,lambda);
toc; 

InterpReg = regularize(pts2,Omega,'[log(r)]',Vh);
Interp1 = -1i/4*Interp1+1/(2*pi)*InterpReg;
Interp2 = -1i/4*Interp2+1/(2*pi)*InterpReg;

xInterp1 = F(mesh2.vtx)+(1-coeff)*k^2*Interp1*x1;
xInterp2 = F(mesh2.vtx)+(1-coeff)*k^2*Interp2*x2;

figure
axis equal
axis tight
colormap('jet')
graph(Vh2,real(xInterp1));


figure
axis equal
axis tight
colormap('jet')
graph(Vh2,real(xInterp2));
