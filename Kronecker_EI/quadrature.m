function [z,w] = quadrature(n)
%
% QUADRATURE - Function that returns the points and weights of  
%              quadrature rules
%
% INPUT: 
%   n - number of quadrature points         
%
% OUTPUT: 
%   z - quadrature point  
%   w - quadrature weight 
%

if n==1
    z=0;
    w=2;
elseif n==2
    %z=[-sqrt(1/3) sqrt(1/3)]; %Gauss
    z=[-1 1]; %Lobatto
    w=[1 1];
elseif n==3
    z=[-sqrt(3/5) 0 sqrt(3/5)];
    w=[5/9 8/9 5/9];
elseif n==4
    z=[-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7) sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];
    w=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36];
elseif n==5
    z=[-sqrt(5+2*sqrt(10/7))/3 -sqrt(5-2*sqrt(10/7))/3 0 sqrt(5-2*sqrt(10/7))/3 sqrt(5+2*sqrt(10/7))/3];
    w=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
end

end
