function [Ke,Kebar,Bm,C,Jacobian] = stiff(h,E,nu,rho,dtime,dr,iel,ico,yg,Me)

%Compute Elasticity matrix
%Assuming Plane Strain Case
 c11 = 1-nu;
 c33= 1/(2-nu);
 
 C = zeros(3);
 Ce = zeros(3);
 Ce(1,1) = c11;
 Ce(2,2) = c11;
 Ce(1,2) = nu;
 Ce(2,1) = nu;
 Ce(3,3) = c33;
 C = (E/((1+nu)*(1-2*nu)))*Ce;
 clear Ce; 
 clear c11;
 clear c33;
 
 %Gauss Quadrature points and weights
 %We use 2 GP for rank sufficient integration
 
 point(1,1) = -sqrt(1/3);
 point(1,2) = -sqrt(1/3);
 point(2,1) = sqrt(1/3);
 point(2,2) = -sqrt(1/3);
 point(3,1) = -sqrt(1/3);
 point(3,2) = sqrt(1/3);
 point(4,1) = sqrt(1/3);
 point(4,2) = sqrt(1/3);
 
 weight(1) = 1.0;
 weight(2) = 1.0;
 weight(3) = 1.0;
 weight(4) = 1.0;
 
 %Shape Functions Bilinear shape functions are used
 for i = 1:4
     psi = point(i,1);
     eta = point(i,2);
    
     %Defining the shape function matrix
     
     N(1,i) = 0.25*(1-psi)*(1-eta);
     N(2,i) = 0.25*(1+psi)*(1-eta);
     N(3,i) = 0.25*(1+psi)*(1+eta);
     N(4,i) = 0.25*(1-psi)*(1+eta);
     
     %Derivative of N with respect to psi
     
     dNxe(1,1,i) = -0.25*(1-eta);
     dNxe(1,2,i) =  0.25*(1-eta);
     dNxe(1,3,i) =  0.25*(1+eta);
     dNxe(1,4,i) = -0.25*(1+eta); 
    
     %derivative of N with respect to eta
     dNxe(2,1,i) = -0.25*(1-psi);
     dNxe(2,2,i) = -0.25*(1+psi);
     dNxe(2,3,i) =  0.25*(1+psi);
     dNxe(2,4,i) =  0.25*(1-psi);
     
 end
 
 elecod = ico(:,iel);
 Jacobian = zeros(2);
 
 ke = zeros(8);
 for i = 1:4
     x = 0;
     y = 0;
     
     for j = 1:4
        elecord(:,j) = yg(:,elecod(j));
        x = x + N(j,i)*elecord(1,j);
        y = y + N(j,i)*elecord(2,j);
     end
      Jacobian(:,:) = dNxe(:,:,i)*elecord(:,:)';
      InvJacobian(:,:) = inv(Jacobian(:,:));
      Be = zeros(2,4);
      for k = 1:4
          for j = 1:4
              Be(1,j,k) = dNxe(1,j,k)*InvJacobian(1,1) + dNxe(2,j,k)*InvJacobian(1,2);
              Be(2,j,k) = dNxe(1,j,k)*InvJacobian(2,1) + dNxe(2,j,k)*InvJacobian(2,2);
          end
      end
 end
 
 for i = 1:4
     t = 1;
     for j = 1:4
         B(1,t,i) = Be(1,j,i);
         B(3,t,i) = Be(2,j,i);
         t = t+1;
         B(2,t,i) = Be(2,j,i);
         B(3,t,i) = Be(1,j,i);
         t = t+1;
     end
     b(:,:) = B(:,:,i);
     btc = b'*C;
     btcb = btc*b;
     ke = ke + btcb;
    
 end
 Ke = ke*det(Jacobian);
 clear ke;
 clear b;
 
 Kebar = (1.0/dtime+dr/2)*Me+dtime*Ke/4.0;
 bnew = zeros(3,8);
 for i=1:4
     b(:,:) = B(:,:,i);
     bnew(:,:) = bnew(:,:) + b(:,:);
 end
 Bm = bnew; %Strain Displacement Matrix
end
