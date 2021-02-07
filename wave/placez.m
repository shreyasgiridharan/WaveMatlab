
%  placez.m

function [gstif]=placez(gstif)

 sdof=size(gstif);  
 sdof = sdof(1);
 for i = 1:4
     gstif(i,:)=zeros(1,sdof);
     gstif(:,i)=zeros(1,sdof);
     gstif(i,i)=1;
 end 
end
 
