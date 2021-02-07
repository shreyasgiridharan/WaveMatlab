function [Se] = relax(h,rho,E,dtime,alpha) 

s = alpha*(rho/dtime+dtime*E/(2*h*h));
con = h*s/8;
for i = 1:8
    Se(i) = con;
end
end