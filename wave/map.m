function [index]=map(nd,nnodel,nvar)
%----------------------------------------------------------
%     Compute system dofs associated with each element 
%-----------------------------------------------------------
 
   k=0;
   for i=1:nnodel
     start = (nd(i)-1)*nvar;
       for j=1:nvar
         k=k+1;
         index(k)=start+j;
       end
   end

 
