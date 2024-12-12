
%%Created on 2/16/2023
%%This file solves diffusive SI model using finite element scheme

tic
clc
clear
clear model
model = createpde(2);
g = geometryFromEdges(model,@circleg);
%g = geometryFromEdges(model,@squareg);


hmax=0.1; %% Mesh size
generateMesh(model, 'Hmax', hmax);   %%% Generate the mesh 


syms u1(t,x,y) u2(t,x,y)
d2=1; d1=0.00001; q=0.5;p=0.5;

pdeeq = [diff(u1,t) - d1*laplacian(u1,[x,y]) + 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p - 1 + u1 - u2; diff(u2,t) - d2*laplacian(u2,[x,y]) + 2*u2 - 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p];

coeffs=pdeCoefficients(pdeeq,[u1 u2]);
%coeffs = pdeCoefficientsToDouble(scoeffs);


specifyCoefficients(model,'m',coeffs.m,'d',coeffs.d, 'c',coeffs.c,'a',coeffs.a,'f',coeffs.f)

applyBoundaryCondition(model,'neumann', 'Edge', 1:4, 'q', zeros(2,2),'g', zeros(2,1))
setInitialConditions(model,[0.8;0.2]);


T=200
tlist=0:1:T;   %%Solution time points

sol=solvepde(model,tlist);

h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,1,T),"ZData",sol.NodalSolution(:,1,T))
pdeplot(model,"XYData",sol.NodalSolution(:,1,T))
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/model2dsSp5', 'epsc');


h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,2,T),"ZData",sol.NodalSolution(:,2,T))
pdeplot(model,"XYData",sol.NodalSolution(:,2,T))
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/model2dsIp5', 'epsc');

toc
