
% Cleaning
clear all
close all
clc

%% Simple illustration

vtx = [0 0 0;0 1 0; 10 20 30; 1 0 1]; 
elt = [1 2 4]; % Third vertex unused
tri = msh(vtx,elt);

assert(isequal(tri.elt,[1,2,3]));


%% Advanced example

% Load mesh
load tetmesh

% Build mesh
mesh1 = msh(X,tet);

% Graphical representation
figure
plot(mesh1)
axis equal
view(45,45)

% Add vertices to tetrahedral mesh
tmp = [zeros(10,3) ; 10*rand(10,3)];
tet = tet + size(tmp,1);
X   = [tmp ; X ; 10*rand(20,3)];

% Build and clean mesh
mesh2 = msh(X,tet);

% Test mesh egality
b = isequal(mesh1,mesh2);

% Graphical representation
figure
plot(mesh2)
axis equal
view(45,45)

% Colours
mesh = mshSquare(20,[1 1]);
ctr  = mesh.ctr;
mesh.col(ctr(:,1)<0)  = 1;    
mesh.col(ctr(:,1)>=0) = 2;

% Graphical representation
figure
plot(mesh)
axis equal
view(45,45)

% Sub-meshing with small translation 
delta     = 1e-2;
mesh1     = mesh.sub(mesh.col==1);
mesh2     = mesh.sub(mesh.col==2);
mesh2.vtx = mesh2.vtx + delta;
meshT     = union(mesh1,mesh2);

% Graphical representation
figure
plot(meshT)
axis equal
view(45,45)

% Clean with specified range
stp   = mesh.stp;
meshT = clean(meshT,0.5*stp(1));

figure;
plot(meshT);
axis equal
view(45,45);


disp('~~> Michto gypsilab !')


