function lowAppSet = LowerApproximation(U,R,X)
% ����X��R�½��ƣ���R_(X)
% Input:  Name   description         type
%          U        ����             array
%          R       �ȼ۹�ϵ�������ԣ� array
%          X        ����             array
% Output: Name   description         tpye
%       lowAppSet  �½���            array

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
