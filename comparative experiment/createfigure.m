function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  x ���ݵ�����
%  YMATRIX1:  y ���ݵľ���

%  �� MATLAB �� 15-May-2023 21:03:50 �Զ�����

% ���� figure
figure1 = figure;

% ���� axes
axes1 = axes('Parent',figure1,'Position',[0.092 0.177 0.875 0.8]);
hold(axes1,'on');

% ʹ�� plot �ľ������봴������
plot(X1,YMatrix1,'Parent',axes1);

% ���� ylabel
ylabel('Reflectance (%)');

% ���� xlabel
xlabel(['Wavelength (nm)';'(d)            ']);

% ȡ�������е�ע���Ա����������� X ��Χ
 xlim(axes1,[300 2600]);
% ȡ�������е�ע���Ա����������� Y ��Χ
% ylim(axes1,[-Inf Inf]);
box(axes1,'on');
% ������������������
set(axes1,'FontSize',24,'XTick',[300 600 900 1200 1500 1800 2100 2400 2600]);
