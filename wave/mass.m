
function [Me] = mass(rho);   
    Te = readtable('massmatrix.txt','Delimiter',',');
    Tn = table2array(Te);
    Me = 1/72.*Tn;
    clear Te;
    clear Tn;
end
   
