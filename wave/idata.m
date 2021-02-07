warning('off','MATLAB:dispatcher:InexactMatch');

nnodel = 4;                % Number of nodes per Element
nnod = 22;                 % Total number of nodes in the system
nel = (nnod-2)/2;          % Total number of elements in the system
nvar = 2;                  % Number of DOFs per node (x-, and y-)(2D)
nplot = 44;                % Node at which displacement is to be known

H = 10.0;                  % Height of Global Element
h = H/nel;                 % Height of Local Element

E(1)    = 1000;            % Young's Modulus of the Material
nu(1) = 0;                 % Poisson's ratio
gama(1) = 10.0;            % Unit weight of Element 9.81 kg ~ 10 kg
rho(1)  = gama(1)/10.0;    % Density of Material
dr = 0.0;                  % Damping coefficient


nstep   = 500;            % Number of steps (for implicit scheme)
timef   = 5.0;            % Final time
time    = 0.0;            % Initialize time


alpha = 1.25;

dtime = timef/nstep;
% fprintf('Time step for implicit scheme:   %g\n',dtime) 

%---------------------------------------------
%  input data for nodal coordinate values (yg)
%---------------------------------------------

ix(:,1) = zeros(nnod*2,1)';

%Vertical Bar Problem
temp = 1;
while(temp<nnod*2)
    ix(temp,1) = 1;
    temp = temp+2;
end
%Constraining the bottom nodes fully
ix(2,1) = 1;
ix(4,1) = 1;

%Nodal Coordinates

ygtable = readtable('InputCoord.dat');
ygconv = table2array(ygtable);
yg = ygconv';
clear ygtable; %Free up memory
clear ygconv;  %Free up memory
%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  ico(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------
icotable = readtable('Icon.txt');
icoconv = table2array(icotable);
ico = icoconv';
clear icotable; %Free up memory
clear icoconv;  %Free up memory

%Finding area of each element
for iel = 1:nel
    inod = ico(:,iel);
    area(iel,1) = abs(yg(1,inod(3))-yg(1,inod(1)))*abs(yg(2,inod(3))-yg(2,inod(1)));
end

is       = ones(nel,1);
 
nnet=nnod*nvar;            % total system dofs  
nvel=nnodel*nvar;          % degrees of freedom per element
 
dtmax = 10000000.00;
for iel=1:nel           % loop over elements to determine maximum time step
%     nd=ico(:,iel);   
%     y=yg(:,iel); 
   iss = is(iel);
   dt = h/sqrt(E(iss)/rho(iss));  
   if(dt < dtmax)
      dtmax = dt;
   end
end

% fprintf('The maximum permissible time step:   %g\n',dtmax) 

