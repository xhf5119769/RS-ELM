function posSet = PositiveRegion(U,C,D)
% PositiveRegion 返回D的C正域，即pos_C(D)
% Input:  Name   description         type
%          U        论域             array
%          C       条件属性          array
%          D       决策属性          array
% Output: Name   description         tpye
%        posSet  正域pos_C(D)        array
% Example:
%     >> U = [1 2 3 4 5]'
%     >> C = [0 1 2;2 1 0; 3 2 4;2 2 1; 3 2 4]
%     >> D = [0 2 2 1 1 ]'
%     >> pos = PositiveRegion(U,C,D)      


%==========================================================================

posSet = [];
U_D_cell = EquivalentClassSet(U,D);
for i = 1:length(U_D_cell)
    X = U_D_cell{i};
    posSet = union(posSet,LowerApproximation(U,C,X));
end
end