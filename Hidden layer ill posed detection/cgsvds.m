function [U,sm,X,V] = cgsvds(A,L) 
%CGSVD Compact generalized SVD of a matrix pair in regularization problems. 
% 
% sm = cgsvd(A,L) 
% [U,sm,X,V] = cgsvd(A,L) ,  sm = [sigma,mu] 
% 
% Computes the generalized SVD of the matrix pair (A,L): 
%    [ A ] = [ U  0 ]*[ diag(sigma)      0    ]*inv(X) 
%    [ L ]   [ 0  V ] [      0       eye(n-p) ] 
%                     [  diag(mu)        0    ] 
% where 
%    U  is  m-by-n ,    sigma  is  p-by-1 
%    V  is  p-by-p ,    mu     is  p-by-1 
%    X  is  n-by-n . 
% 
% It is assumed that m >= n >= p, which is true in regularization problems. 
 
% Reference: C. F. Van Loan, "Computing the CS and the generalized 
% singular value decomposition", Numer. Math. 46 (1985), 479-491. 
 
% Per Christian Hansen, IMM, 12/19/97. 
 
% Initialization. 
[m,n] = size(A); [p,n1] = size(L); 
if (n1 ~= n | m < n | n < p) 
  error('Incorrect dimensions of A and L') 
end 
 
% Call Matlab's GSVD routine. 
[U,V,W,C,S] = gsvd(full(A),full(L),0); 
sm = [diag(C(1:p,1:p)),diag(S(1:p,1:p))]; 
 
% Finalize. 
if (nargout < 2) 
   U = sm; 
else 
   % Full decomposition. 
   X = inv(W'); 
end