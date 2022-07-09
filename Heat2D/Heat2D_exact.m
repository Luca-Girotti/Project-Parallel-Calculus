
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Exact solution of the 1D heat conduction equation.              %
%             Written by W. Boscheri, 03/02/2020                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Heat conduction coefficient
kappa = 5e-3;

% Set boundary conditions
TL = 1; % left  temperature
TR = 0.1;  % right temperature

% Set final time at which we want to compute the solution
t = 1.0;

% Computational domain
xL   = -0.5;   % left  boundary
xR   = +0.5;   % right boundary
xD   = 0.0;    % position of the initial discontinuity
IMAX = 1000;   % total number of cells
x    = magic(IMAX);  % coordinates of cell barycenters

% EXACT SOLUTION
Te = zeros(IMAX,1);  % solution vector
for i = 1:IMAX
    for j = 1:IMAX
        Te(i,j) =  0.5 * (TR+TL) + 0.5 * erf((x(i,j)-xD)/(2*sqrt(kappa*t))) * (TR-TL);
    end
end



% Write to file for export to Tecplot data structure
fileID = fopen('Heat1D-exact-5e-3.dat','w');
fprintf(fileID, ' VARIABLES = "x" "Te" \n');
fprintf(fileID, ' ZONE I = %d, DATAPACKING=BLOCK \n',IMAX);
fprintf(fileID,'%f \n',x(1:IMAX));
fprintf(fileID,'%f \n',Te(1:IMAX));
fclose(fileID);