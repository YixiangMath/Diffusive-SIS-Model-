
%%Created on 2/16/2023
%%This file solves diffusive SI model using finite element scheme

tic
clc
clear
clear model
model = createpde(2);
g = geometryFromEdges(model,@circleg);
%g = geometryFromEdges(model,@squareg);


hmax=0.02; %% Mesh size
generateMesh(model, 'Hmax', hmax);   %%% Generate the mesh 

% h=figure
% %pdeplot(pdem, 'NodeLabels', 'on') %%%Plot the mesh with nodes label
% pdeplot(model) %% Plot mesh
% saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/mesh', 'epsc');

syms u1(t,x,y) u2(t,x,y)
d1=1; d2=10^-7; q=0.5;p=1;

%pdeeq = [diff(u1,t) - d1*laplacian(u1,[x,y]) + 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p - (1.345+y+x^2) + u1 - u2; diff(u2,t) - d2*laplacian(u2,[x,y]) + 2*u2 - 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p];

pdeeq = [diff(u1,t) - d1*laplacian(u1,[x,y]) + 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p - 1 + u1 - u2; diff(u2,t) - d2*laplacian(u2,[x,y]) + 2*u2 - 2*(1.5+sin(pi*x)*sin(pi*y))*u1^q*u2^p];

coeffs=pdeCoefficients(pdeeq,[u1 u2]);
%coeffs = pdeCoefficientsToDouble(scoeffs);


specifyCoefficients(model,'m',coeffs.m,'d',coeffs.d, 'c',coeffs.c,'a',coeffs.a,'f', coeffs.f)

applyBoundaryCondition(model,'neumann', 'Edge', 1:4, 'q', zeros(2,2),'g', zeros(2,1))
setInitialConditions(model,[0.8;0.2]);


T=600
tlist=0:1:T;   %%Solution time points

sol=solvepde(model,tlist);


nodes=sol.Mesh.Nodes;
[m,n]=size(nodes);
ind=zeros(n,1);
H=zeros(n, 1);
X=sol.NodalSolution(:,1,T);

for I=1:n
    x=nodes(1,I);
    y=nodes(2,I);
    H(I)=(1/(1.5+sin(pi*x)*sin(pi*y)))^(1/q);
    if H(I)-X(I)<0.00001
        ind(I)=1;
    end
end




h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,1,T),"ZData",sol.NodalSolution(:,1,T))
%pdeplot(model,"XYData",H-sol.NodalSolution(:,1,T))
%pdeplot(model,"XYData",H-sol.NodalSolution(:,1,T),"ZData",H-sol.NodalSolution(:,1,T))
pdeplot(model,"XYData",ind)
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/model2diS', 'epsc');


h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,2,T),"ZData",sol.NodalSolution(:,2,T))
pdeplot(model,"XYData",sol.NodalSolution(:,2,T))
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/model2diI', 'epsc');



%%%%The following code draws R, high and low risk domains. Use a finer mesh
%%%%to draw the graphs
% % nodes=sol.Mesh.Nodes;
% % [m,n]=size(nodes);
% % ind=zeros(n,1);
% % R=zeros(n, 1);
% % for I=1:n
% %     x=nodes(1,I);
% %     y=nodes(2,I);
% %     ind(I)=sign(1.5+sin(pi*x)*sin(pi*y)-1);
% %     R(I)=1.5+sin(pi*x)*sin(pi*y);
% % end
% % 
% % h=figure
% % %pdeplot(model,"XYData",sol.NodalSolution(:,1,T),"ZData",sol.NodalSolution(:,1,T))
% % pdeplot(model,"XYData",ind)
% % saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/ind', 'epsc');
% % 
% % 
% % h=figure
% % %pdeplot(model,"XYData",sol.NodalSolution(:,1,T),"ZData",sol.NodalSolution(:,1,T))
% % pdeplot(model,"XYData",R)
% % saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/R', 'epsc');


toc
