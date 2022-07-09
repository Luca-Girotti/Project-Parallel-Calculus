clear all
close all
clc
%Solving the Steady State 2D Heat Conduction Equation
%Length of Domain in x and y directions (unit square)
Lx=input("enter value of a");
Ly=input("enter value of b");  
%No. of grid points 
nx=1+input("enter no.of grids along x direction");
ny=1+input("enter no.of grids along y direction");
%Creating the mesh
x=linspace(0,Lx,nx);
y=linspace(0,Ly,ny);
dx=x(2)-x(1);
dy=abs(y(2)-y(1));
%Initial Conditions
T=298*ones(nx,ny);
K1 = dy^2/(2*((dx)^2 + (dy^2)));
K2 = dx^2/(2*((dx^2) + (dy^2)));
%Boundary Conditions
T(1,:) = 0; %First row or Bottom row of nodes in the matrix 
T(end,:) = 0; %Last row or Top row of nodes in the matrix 
T(:,1) = 100; %First column or Left side column of nodes in the matrix 
T(:,end) = 0; %Last column or Right side column of nodes in the matrix
T_old=T; %Initializing the values for T_old for the first loop
tol=1e-4; %tolerance limit to attain convergence
error=9e9; %initial error to execute while loop
%Code for Solvers
%Jacobi Method
count=1;
while(error>tol)
for i=2:(nx-1)
for j=2:(ny-1)
h =K1*(T_old(i-1,j)+T_old(i+1,j));
v = K2*(T_old(i,j-1)+T_old(i,j+1));
T(i,j) = h+v;
end
end
error=max(abs(T(i,j)-T_old(i,j)));
T_old=T;
count=count+1;
end
% Plotting the Results
figure(1)
[a,b] = contourf(x,y,T); 
clabel(a,b); colorbar;
colormap(jet); 
xlabel('X Axis');
ylabel('Y Axis');
text = sprintf('2D steady state Heat Conduction n No. of iterations: %d',count);
title(text);
%analytical sol1
a=Lx;
b=Ly;
TA=298*ones(nx,ny);
%Boundary Conditions
TA(1,:) = 0; %First row or Bottom row of nodes in the matrix 
TA(end,:) = 0; %Last row or Top row of nodes in the matrix 
TA(:,1) = 100; %First column or Left side column of nodes in the matrix 
TA(:,end) = 0; %Last column or Right side column of nodes in the matrix
sum=0;
      for f=2:(nx-1)       
            for g=2:(ny-1) 
                for n=0:10
                    th=((b-g)*(2*n+1)*pi)/a;
                    tg=((2*n+1)*pi*f)/a;
                    tj=((2*n+1)*pi*b)/a;
                    tr=1/(2*n+1);
sum=sum+(sin(tg)*sinh(th)*tr)/sinh(tj)
                end        
        TA(f,g)= (((4*100)*sum)/pi);
 sum=0;
            end  
      end 
figure(2)
[a,b] = contourf(x,y,TA); 
clabel(a,b); colorbar;
colormap(jet); 
xlabel('X Axis');
ylabel('Y Axis');
text = sprintf('2D Transient Heat Conduction using analytical Method n No. of iterations: %d',count);
title(text);
TAA=298*ones(nx,ny);
%analytical sol2
%Boundary Conditions
TAA(1,:) = 0; %First row or Bottom row of nodes in the matrix 
TAA(end,:) = 0; %Last row or Top row of nodes in the matrix 
TAA(:,1) = 100; %First column or Left side column of nodes in the matrix 
TAA(:,end) = 0; %Last column or Right side column of nodes in the matrix
sumA=0;
      for fA=2:(nx-1)       
            for gA=2:(ny-1) 
                for nA=1:10
                    thA=(gA*nA*pi)./a;
                    tgA=(nA*pi*fA)./a;
                    tjA=(nA* (pi)*b)/a;
                    trA=1/nA;
sumA=sumA+(sin(tgA)*sinh(thA)*trA*(1-(1^nA)))/sinh(tjA);
                end        
        TAA(fA,gA)= (((2*100)*sumA)/pi);
 sumA=0;
            end  
      end 
figure(3)
[a,b] = contourf(x,y,TAA); 
clabel(a,b); colorbar;
colormap(jet); 
xlabel('X Axis');
ylabel('Y Axis');
text = sprintf('2D Transient Heat Conduction using analytical Method n No. of iterations: %d',count);
title(text);