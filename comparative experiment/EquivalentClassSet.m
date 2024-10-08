function equSetCell = EquivalentClassSet(U,R)
% EquivalentClassSet 返回基于知识R对U进行分类的所有等价类构成的集合，即U/R
% Input:  Name     description         type
%          U          论域             array
%          R         等价关系（属性）   array
% Output: Name      description        tpye
%       equSetCell  等价类集合[x]_R    cell array

% =========================================================================

% R中若有多个属性则通过编码将其转化成一个整数值以便比较
% 编码值相同的将被归为一类
if(size(R,1)~=size(U,1))
    error('RS:error','wrong input argument size => size(R,1)~=size(U,1)');
end
N = max(R(:))+1;
Rcode = R*N.^[size(R,2)-1:-1:0]';
uniqRcode = unique(Rcode);
numClass =  size(uniqRcode,1);
% 将基于R的各个分类组成一个集合（用cell array数据类型存储）
equSetCell = cell(1,numClass);
for i = 1:numClass
    ind = (Rcode==uniqRcode(i));
    equSetCell{i} = U(ind);
end   
end
