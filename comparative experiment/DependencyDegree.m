function gama = DependencyDegree(U,C,D)
% DependencyDegree ����D��C��������gama
% Input:  Name     description         type
%          U           ����            array
%          C         ��������          array
%          D         ��������          array
% Output: Name     description         tpye
%         gama        ������           double

% Author: Neptune_zx 
% Email:  553680533@qq.com
% Time:   2011/5/19
% =========================================================================

if(isempty(C))
    gama = 0;
    return;
end
posCD = PositiveRegion(U,C,D);
gama = length(posCD)/length(U);
end