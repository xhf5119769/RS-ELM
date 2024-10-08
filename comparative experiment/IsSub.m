function bool = IsSub(A,B)
% IsSub 判断集合A是否是集合B的子集
% Input:  Name     description         type
%          A        集合A              array
%          B        集合B              array
% Output: Name     description         type
%         bool       真假             ligical

% Author: Neptune_zx 
% Email:  553680533@qq.com
% Time:   2011/5/19
% =========================================================================
if(all(ismember(A,B)))
    bool = true;
else
    bool = false;
end