function equSetCell = EquivalentClassSet(U,R)
% EquivalentClassSet ���ػ���֪ʶR��U���з�������еȼ��๹�ɵļ��ϣ���U/R
% Input:  Name     description         type
%          U          ����             array
%          R         �ȼ۹�ϵ�����ԣ�   array
% Output: Name      description        tpye
%       equSetCell  �ȼ��༯��[x]_R    cell array

% =========================================================================

% R�����ж��������ͨ�����뽫��ת����һ������ֵ�Ա�Ƚ�
% ����ֵ��ͬ�Ľ�����Ϊһ��
if(size(R,1)~=size(U,1))
    error('RS:error','wrong input argument size => size(R,1)~=size(U,1)');
end
N = max(R(:))+1;
Rcode = R*N.^[size(R,2)-1:-1:0]';
uniqRcode = unique(Rcode);
numClass =  size(uniqRcode,1);
% ������R�ĸ����������һ�����ϣ���cell array�������ʹ洢��
equSetCell = cell(1,numClass);
for i = 1:numClass
    ind = (Rcode==uniqRcode(i));
    equSetCell{i} = U(ind);
end   
end
