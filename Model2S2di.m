
%%Created on 2/16/2023
%%This file solves diffusive SI model using finite element scheme

tic
clc
clear
clear model
model = createpde(2);
g = geometryFromEdges(model,@circleg);



hmax=0.02; %% Mesh size
generateMesh(model, 'Hmax', hmax);   %%% Generate the mesh 

% h=figure
%pdegplot(model, 'FaceLabels', 'on') %%%Plot the mesh with nodes label
% pdeplot(model) %% Plot mesh
% saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/mesh', 'epsc');

syms u1(t,x,y) u2(t,x,y)
d1=1; d2=10^-5; 

c=[d1; 0; d1; d2; 0; d2];



specifyCoefficients(model,'m',0,'d', 1, 'c',c,'a',0,'f',@fcoeffunction)

applyBoundaryCondition(model,'neumann', 'Edge', 1:4, 'q', zeros(2,2),'g', zeros(2,1))

setInitialConditions(model,[0.8;0.2]);


T=1000
tlist=0:1:T;   %%Solution time points

sol=solvepde(model,tlist);

q=0.5;
p=1;

nodes=sol.Mesh.Nodes;
[m,n]=size(nodes);
ind=zeros(n,1);
H=zeros(n, 1);

X=sol.NodalSolution(:,1,T);

for I=1:n
    x=nodes(1,I);
    y=nodes(2,I);
    %H(I)=(1/(1.5+sin(pi*x)*sin(pi*y)))^(1/q);
   if  x<0 
       A=0.5+0.4*x^2;
   elseif x< 0.25 && x>=0
       A=0.5;
   elseif x<0.5 && x>=0.25
       A=0.5+0.4*(x-0.25)^2;
   elseif x>=0.5
       A=0.5+1.6*(x-0.625)^2;
   end

   if  y<0
       B=0.5+0.4*y^2;
   elseif y< 0.25 && y>=0
       B=0.5;
   elseif y<0.5 && y>=0.25
       B=0.5+0.4*(y-0.25)^2;
   elseif y>=0.5
       B=0.5+1.6*(y-0.625)^2;
   end
   H(I)=((A*B+0.1)/0.5)^(1/q);
    if H(I)-X(I)<0.01
        ind(I)=1;
    end
end




h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,1,T),"ZData",sol.NodalSolution(:,1,T))
%pdeplot(model,"XYData",H-sol.NodalSolution(:,1,T))
%pdeplot(model,"XYData",H-sol.NodalSolution(:,1,T),"ZData",H-sol.NodalSolution(:,1,T))

pdeplot(model,"XYData",ind)
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/M2S2ind1', 'epsc');
%saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/M2S2ind2', 'epsc');


h=figure
pdeplot(model,"XYData",sol.NodalSolution(:,1,T))
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/M2S2diS', 'epsc');


%pdeplot(model,"XYData",H-sol.NodalSolution(:,1,T))


h=figure
%pdeplot(model,"XYData",sol.NodalSolution(:,2,T),"ZData",sol.NodalSolution(:,2,T))
pdeplot(model,"XYData",sol.NodalSolution(:,2,T))
saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-2-SI with nonlinear incidence-PengRui-Salako/M2S2diI', 'epsc');

save('Amodel2di')


pdeplot(model,"XYData",sol.NodalSolution(:,2,T), "ZData",sol.NodalSolution(:,2,T))



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


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Right hand side of SI model

toc



function f = fcoeffunction(region,state)

N = 2; % Number of equations
q=0.5;
p=1;
nr = length(region.x); % Number of columns
f = zeros(N,nr); % Allocate f

X=region.x;
Y=region.y;
gammamu=ones(1, nr);
A=0;B=0;
 for I=1:nr
   % beta(I)=3+2*sin(pi*X(I)).*sin(pi*Y(I));
  
   if  X(I)<0 
       A=0.5+0.4*X(I)^2;
   elseif X(I)< 0.25 && X(I)>=0
       A=0.5;
   elseif X(I)<0.5 && X(I)>=0.25
       A=0.5+0.4*(X(I)-0.25)^2;
   elseif X(I)>=0.5
       A=0.5+1.6*(X(I)-0.625)^2;
   end

   if  Y(I)<0
       B=0.5+0.4*Y(I)^2;
   elseif Y(I)< 0.25 && Y(I)>=0
       B=0.5;
   elseif Y(I)<0.5 && Y(I)>=0.25
       B=0.5+0.4*(Y(I)-0.25)^2;
   elseif Y(I)>=0.5
       B=0.5+1.6*(Y(I)-0.625)^2;
   end

   gammamu(I)=A*B;

%     if X(I)>0 &&  Y(I)>0
%         gammamu(I)=3+32*(X(I)-X(I)*X(I))*(Y(I)-Y(I)*Y(I));
%     end
    
end

% Now the particular functional form of f
f(1,:) = 1- state.u(1,:)  - 0.5*state.u(1,:).^q.*state.u(2,:).^p + gammamu.*state.u(2,:);   %-0.015*state.u(1,:).*state.u(2,:)./(1+0.025*state.u(2,:));
f(2,:) = 0.5*state.u(1,:).^q.*state.u(2,:).^p - (gammamu+0.1).*state.u(2,:);
end


