function equSet = EleEquivalentSet(U,x,R)
% EleEquivalentSet 返回包含元素x的且基于R分类下的等价集[x]_R
% Input:  Name   description         type
%          U        论域             array
%          x        元素             array
%          R        等价关系（属性）  array
% Output: Name   description         tpye
%        equSet   等价集[x]_R        array

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

