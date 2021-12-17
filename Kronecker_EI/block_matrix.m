function At = block_matrix(A,u,dimx,order)
%
% BLOCK_MATRIX - Function that returns the block-wise sparse matrix 
%                necessary  to compute the Varphi_p(A)*u
%
% INPUT: 
%   A     - full 2D matrix
%   u     - vector 
%   dimx  - size of vector u
%   order - order of the Varphi function 
%
% OUTPUT: 
%   At - amplified matrix 
%

%Allocate memory
At=spalloc(dimx+order,dimx+order,nnz(A)+nnz(u));

%Assign the blocks
At(1:dimx,1:dimx)=A;
At(1:dimx,dimx+1)=u;
if order>1
  At(dimx+1:dimx+order-1,dimx+2:dimx+order)=eye(order-1);
end
  
end
