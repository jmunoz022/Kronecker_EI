function C = action(A,B,p)
% ACTION - Function that returns the action (A^p)*B
%
% INPUT: 
%   A - a square matrix
%   B - a vector or a square matrix
%   p - an integer
%
% OUTPUT: 
%   C - the action (A^p)*B         

C=A*B;
for i=2:p
    C=A*C;
end

end

