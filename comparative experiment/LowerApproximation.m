function lowAppSet = LowerApproximation(U,R,X)
% 返回X的R下近似，即R_(X)
% Input:  Name   description         type
%          U        论域             array
%          R       等价关系（即属性） array
%          X        属性             array
% Output: Name   description         tpye
%       lowAppSet  下近似            array

% =========================================================================
lowAppSet = [];
if(~IsSub(X,U))
    warning('RS:warning','set X is not sub set of U.');
return;
end
U_R_cell = EquivalentClassSet(U,R);
for i=1:length(U_R_cell)
    if(IsSub(U_R_cell{i},X))
      lowAppSet  = union(lowAppSet,U_R_cell{i}); 
    end
end
        
end
