%events function
function [value,isterminal,direction] = events_control(~,y)
value=y(12);%���嵥Ԫ��Ŀ
isterminal=1;%�����굥��ֹͣ��������ֹͣ��Ϊ0
direction=-1;%����������Ϊ0��Ϊ1�����ɸ�������0��������
end