function targetValue = TargetOptFcn(r,U,C,D)
% Input:  Name     description         type
%          r         ���Ա���          array
%          U           ����            array
%          C         ��������          array
%          D         ��������          array
% Output: Name      description        tpye
%      targetValue    Ŀ��ֵ           double

%==========================================================================

    if(length(r)~= size(C,2))
        error('RS:error','wrong binary input r.');
    end
    r_bin =round(r);
    C_reduct = C(:,r_bin==1);
    % ��һ��Ŀ�������Ŀ��D����C_reduct��������
    target1 = DependencyDegree(U,C_reduct,D);
    % �ڶ���Ŀ���Ǿ�����С����C��1�ĸ�����Ҳ���Ǿ�����Լ��C
    target2 = (length(r_bin) - sum(r_bin))/length(r_bin);
    % ������Ŀ������Ⱦɫ��Lr�е����־����ܵĽӽ�1����0
    target3 = sum(abs(r.*(r-1)))/length(r);
    % ע������Ϊ�˸��õ��Ż������3��Ŀ������˼�Ȩ����ʵ����target1 +target2
    % �����������Ŀ�꺯��
    % ga����Ĭ������С��Ŀ�꺯������������ǰ2��Ŀ��Ҫ�Ӹ�����תΪ��С��
    targetValue = -5*target1 - 2*target2 + target3;   
end
