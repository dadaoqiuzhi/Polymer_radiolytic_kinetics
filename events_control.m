%events function
function [value,isterminal,direction] = events_control(~,y)
value=y(12);%单体单元数目
isterminal=1;%消耗完单体停止迭代，不停止则为0
direction=-1;%由正到负减为0，为1则是由负到正，0则是任意
end