function gama = DependencyDegree(U,C,D)
% DependencyDegree 返回D对C的依赖度gama
% Input:  Name     description         type
%          U           论域            array
%          C         条件属性          array
%          D         决策属性          array
% Output: Name     description         tpye
%         gama        依赖度           double

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