function equSet = EleEquivalentSet(U,x,R)
% EleEquivalentSet ���ذ���Ԫ��x���һ���R�����µĵȼۼ�[x]_R
% Input:  Name   description         type
%          U        ����             array
%          x        Ԫ��             array
%          R        �ȼ۹�ϵ�����ԣ�  array
% Output: Name   description         tpye
%        equSet   �ȼۼ�[x]_R        array

% =========================================================================

if(~ismember(x,U))
    equSet = [];
    warning('RS:warning',['the element: ' num2str(x) ' is not in set U.']);
    return;
end
equSetCell = EquivalentClassSet(U,R);
for i = 1:length(equSetCell)
    if(ismember(x,equSetCell{i}))
        equSet = equSetCell{i};
        return;
    end
end
end

