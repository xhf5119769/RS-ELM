function bool = IsSub(A,B)
% IsSub �жϼ���A�Ƿ��Ǽ���B���Ӽ�
% Input:  Name     description         type
%          A        ����A              array
%          B        ����B              array
% Output: Name     description         type
%         bool       ���             ligical

% Author: Neptune_zx 
% Email:  553680533@qq.com
% Time:   2011/5/19
% =========================================================================
if(all(ismember(A,B)))
    bool = true;
else
    bool = false;
end