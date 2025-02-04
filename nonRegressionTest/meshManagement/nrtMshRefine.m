%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtMshRefine.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Mesh refinement algorithm for triangular      |
%|  `---'  |                and segment meshes                            |                  
%+========================================================================+

% Cleaning
clear all
close all
clc

% Create mesh
Nvtx = 50;
mesh = mshSquare(Nvtx,[2 2]);

% Colours
mesh.col(1:20)                 = -1;
mesh.col(mesh.nelt+(-20:0)) = 1;

% Graphical representation
figure
plot(mesh)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refine all element with midpoint algorithm
tic
meshr = midpoint(mesh);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Midpoint algorithm for selected indices
I = find(mesh.col~=0);
tic
meshr = midpoint(mesh,I);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and fixed order
ord = (1 + mesh.vtx(mesh.elt(:,1),1))/2;
ord = floor(4*ord);
tic
meshr = refine(mesh,ord);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and function order
fct = @(X) floor(5*(1+X(:,1))/2);
tic
meshr = refine(mesh,fct);
toc

% Graphical representation
figure
plot(meshr)
hold on
plotNrm(meshr,'r')
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and fixed edge length
lambda = 0.08;
tic
meshr = refine(mesh,lambda);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')


%% Same tests with segment mesh


% Cleaning
clear all
close all
clc

% Create mesh
N = 20;
mesh = meshCurve(normalParam(spirale),N);

% Colours
mesh.col(1:3)                 = -1;
mesh.col(mesh.nelt+(-3:0)) = 1;

% Graphical representation
figure
plot(mesh)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refine all element with midpoint algorithm
tic
meshr = midpoint(mesh);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Midpoint algorithm for selected indices
I = find(mesh.col~=0);
tic
meshr = midpoint(mesh,I);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and fixed order
ord = (1 + mesh.vtx(mesh.elt(:,1),1))/2;
ord = floor(4*ord);
tic
meshr = refine(mesh,ord);
toc

% Graphical representation
figure
plot(meshr)
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and function order
fct = @(X) floor(5*(1+X(:,1))/2);
tic
meshr = refine(mesh,fct);
toc

% Graphical representation
figure
plot(meshr)
hold on
plotNrm(meshr,'r')
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z')

% Refinement with recursive midpoint algorithm and fixed edge length
lambda = 0.08;
tic
meshr = refine(mesh,lambda);
toc
assert(max(meshr.ndv) <= lambda);

% Graphical representation
figure
plot(meshr)
axis equal 
xlabel('X'); ylabel('Y'); zlabel('Z')

disp('~~> Michto gypsilab !')

