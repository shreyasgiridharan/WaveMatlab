
                             %2D FEM Analysis
                        %Institut für Geotechnik

                                %wave.m

                            %Sub Programs BOM
% 1. Idata.m    For initialising input data
%   1.1 ICon.txt    Text file where element connectivity data is stored
%   1.2 InputCoord.dat  Data file where the nodal coordinates are stored
% 2. map.m      Compute system DOFs
% 3. mass.m     Initialize mass matrix
%   3.1 massmatrix.txt  Text file where the entries of matrix are stored
% 4. placez.m   For applying displacement boundary condition
% 5. relax.m    For computing the relaxation term S
% 6. stiff.m    For computing K,B,Jacobian using 2X2 Gauß quadrature
% 7. wave.m     Main run file for calling all the above functions and plot


clear
clc
run_time=cputime;

% Initialise Input data
% Change data in idata.m for different results
idata
% End of Initialisation

M=sparse(nnet,nnet);     % allocate global mass matrix
K=sparse(nnet,nnet);     % allocate global stiffness matrix
Kbar=sparse(nnet,nnet);  % allocate tangential stiffness matrix
S=zeros(nnet,1);
f = zeros(1,nnet)';      % allocate global load vector 
u = zeros(1,nnet)';      % allocate displacement vector 
p0= zeros(1,nnet)';      % allocate momentum vector 
stress = zeros(3,nel );  % allocate stresses 
index=zeros(nvel,1);   
nd=zeros(nnodel,1);

%define the center of the element
%    nel = 1;
for iel=1:nel
    nd = ico(:,iel);
    index=map(nd,nnodel,nvar);
    iss =is(iel);
    [Me] = mass(rho(iss));
    [Se] = relax(h,rho(iss),E(iss),dtime,alpha);
    [Ke,Kebar,Bm,C,Jacobian] = stiff(h,E(iss),nu(iss),rho(iss),dtime,dr,iel,ico,yg,Me);
    M(index(:),index(:)) = M(index(:),index(:)) + Me(:,:); 
    S(index(:)) = S(index(:)) + Se(:);    
    K(index(:),index(:)) = K(index(:),index(:)) + Ke(:,:); 
    Kbar(index(:),index(:)) = Kbar(index(:),index(:)) + Kebar(:,:); 
end

[Kbar] = placez(Kbar); % Applying displacement boudary conditions

t = [time];  
ut = [u(nplot)];    % checking for displacement at the last node 
sig = [stress(2,1)];   % stress at 1st element y-stress 

for i  = 1:nstep
    f(42)=-1;
    f(44)=-1;
    q = K*u; 
    R0 = dtime/2*(f-q) + p0; 
    R0(1:4)=0;
    du = R0./S; %Initial estimate for change in displacement psi
    du(1:4) = 0.0; %Base is fixed
    
    iter = 0;     
    
    while (iter < 100)
        iter = iter + 1;
        ddu = (R0-Kbar*du)./S;               % Here the diagonal term of S is divided into  the corresponding term of the residual load,  (scalar operation nnod times)
        du = du +ddu;                        % Correction to change in displacement
        du(1:4) = 0.0;                       % Base being fixed
%         if(sqrt(ddu'*ddu/(du'*du+1.e-04)) < 1.e-04) 
%            iter = 100;
%         end     
    end
    u = u + du;
    p0 = 2*M*du/dtime - p0;  
    p0(1:4) = 0.0;
    
    for iel=1:nel           % Updating stresses in elements
        elecod = ico(:,iel);
         temp = 1;
         for pointer = 1:4
            Id = (elecod(pointer)-1)*2;
            unew(temp) = du(Id+1);
            temp = temp+1;
            unew(temp) = du(Id+2);
            temp = temp+1;
         end
         strain = Bm*unew'*det(Jacobian); 
         dstress = C*strain;
         stress(:,iel) = stress(:,iel) + dstress(:);
    end
    
    time = time + dtime;        % Incremental time
    t = [ t time];              % Time vector appended
    ut = [ut u(nplot)];         % Displacement vector appended
    sig = [sig stress(2,1)];    % Stress (ele 1, sigy) appended
    
end
figure
plot(t,ut);
xlabel('time t [s]')
ylabel('displacement u [m]')
title ('Displacement over time')

figure
plot(t,sig);
xlabel('time t [s]')
ylabel('stress \sigma [Pa]')
title ('Stress distribution over time')

%END OF PROGRAM