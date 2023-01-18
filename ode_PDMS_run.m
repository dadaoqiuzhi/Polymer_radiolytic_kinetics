%run sode of analysis degradation dynamics of PDMS
%version 1, 2020.12.17
%��ǿ���������뻯ѧ�о���210�ң�forliubinqiang@163.com,15828636974

%%
%options = odeset('RelTol',1e-4,'AbsTol',1e-8);%default1e-3,1e-6 respectively,AbsTol�ɶ�ÿ��΢�ַ�����ָ��[1e-4 1e-4 1e-5]
%p0act�Խ��Ӱ��ܴ�
fprintf('\n����Ӧ�޸�PDMS_ode.m�ļ���ode_PDMS_run.m�ļ���\n');
fprintf('\n��������Ϊ1e-2��1e-3�����¶�Ϊ298Kʱ�����ͼ��������ݲ���ͼ\n');
eventans=input('\n�Ƿ���е�Ч���嵥Ԫ��Ϊ0��ֹͣģ����ƣ������򲻻���з����Ӧģ�⣺y/n\n','s');
Dosemax=input('\n����������Ч���տ������Gy��������Ҫ����ķ����ЧӦ��Ч��������\n');

global I T tmax on_off1 on_off2 on_off3 on_off4 on_off5 on_off6 on_off7 on_off8 on_off8 on_off10 on_off11 on_off12 on_off13
on_off1=1;on_off2=1;on_off3=1;on_off4=1;on_off5=1;on_off6=1;on_off7=1;on_off8=1;on_off9=1;on_off10=1;on_off11=1;on_off12=1;on_off13=1;
I=input('\n��������������Gy/s��\n');% Gy/s
T=input('\n��������䳡ʵ���¶�K��\n');% K
if T>=600
    fprintf('\n����ĭ����600 K�������Ƚ��⣬�������ʱ������ܲ�����\n')
    warndlg('�¶ȳ���600 K�Ƚ��¶ȣ�������������')
end
Dosemax2=input('\n�����������ռ���Gy��������ʱ�佫�رշ��տ����ЧӦ��\n');
tmax=Dosemax2/I;

fprintf('\n���������У���ȴ�......\n')
load('ini_condition.mat')%�����ֵ����,save('ini_condition', '-tabs')�޸����ݺ�洢
ini_conditiondata=cell2mat(ini_condition(2,:));%��΢�ַ���������أ�dy=PDMS_ode(t,y)

if strcmpi(eventans,'y')
    %ָ��������ֲ����и���������������ض��������д˹���,ode23sû�У�������refine�����ݵ㣬ͨ����Ϊ0�¼��������Ƶ�������
    options = odeset('NonNegative',[1:1:length(ini_conditiondata)],'Refine',5,'Events',@events_control);
else
    options = odeset('NonNegative',[1:1:length(ini_conditiondata)],'Refine',5);
end
%��ѡ�����ode15s,ode23s,ode23t,ode23tb,ode45,ode23,ode113;��������ʱ��s����ֵ�����ĸ�����΢�ַ��������һ��,����ѡ��
tmaxgamma=Dosemax/I;
[Time,Y]=ode15s(@PDMS_ode,[0 tmaxgamma],ini_conditiondata,options);
Dose=Time*I/1000;%kGy����

%%
%����Y�е�0
for i=1:size(Y,2)
    for j=2:size(Y,1)
        if Y(j,i)==0 && j>1
            if j==2
                Y(j,i)=Y(j+1,i);
            elseif j==size(Y,1)
                Y(j,i)=Y(j-1,i);
            else
                Y(j,i)=Y(j-1,i)+(Y(j+1,i)-Y(j-1,i))*(Dose(j)-Dose(j-1))/(Dose(j+1)-Dose(j-1));
            end
        end
    end
end
%%
figure(1);
plot(Dose,Y(:,1),'-+r','Linewidth',2,'markersize',8);
hold on;
plot(Dose,Y(:,2),'-ok','Linewidth',2,'markersize',8);
plot(Dose,Y(:,3),':*g','Linewidth',2,'markersize',8);
plot(Dose,Y(:,4),'-.xm','Linewidth',2,'markersize',8);
plot(Dose,Y(:,5),'-sc','Linewidth',2,'markersize',8);
plot(Dose,Y(:,6),'--db','Linewidth',2,'markersize',8);
plot(Dose,Y(:,13),'-vy','Linewidth',2,'markersize',8);
hold off;
legend('H','CH3','P1','P2','P3','P4','P0act');
xlabel('Dose (kGy)');
ylabel('Yield (mol/L)')

figure(2);
plot(Dose,Y(:,27),':pg','Linewidth',2,'markersize',8,'MarkerEdgeColor','g');
hold on
plot(Dose,Y(:,28),'-.hr','Linewidth',2,'markersize',8);
plot(Dose,Y(:,29),'-+m','Linewidth',2,'markersize',8);
plot(Dose,Y(:,25),':vy','Linewidth',2,'markersize',8);
plot(Dose,Y(:,26),'-sc','Linewidth',2,'markersize',8);
plot(Dose,Y(:,10),'-ok','Linewidth',2,'markersize',8);
plot(Dose,Y(:,24),'-.<r','Linewidth',2,'markersize',8);
plot(Dose,Y(:,11),':*b','Linewidth',2,'markersize',8);
hold off
legend('CH4','H2','C2H6','Si-H','Si-OH','crosslinking-H','crosslinking-Y','scission');
xlabel('Dose (kGy)');
ylabel('Yield (mol/L)')

figure (3)
plot(Dose,Y(:,12),':pg','Linewidth',2,'markersize',8,'MarkerEdgeColor','g');
hold on
plot(Dose,Y(:,23),'-.hr','Linewidth',2,'markersize',8);
plot(Dose,Y(:,22),'-^b','Linewidth',2,'markersize',8);
hold off
legend('P0-equivalent','P0-neat','P3-P4');
xlabel('Dose (kGy)');
ylabel('Yield (mol/L)')

figure(4);
plot(Dose,Y(:,14),'-+r','Linewidth',2,'markersize',8);
hold on;
plot(Dose,Y(:,15),'-ok','Linewidth',2,'markersize',8);
plot(Dose,Y(:,16),':*g','Linewidth',2,'markersize',8);
plot(Dose,Y(:,17),'-.xm','Linewidth',2,'markersize',8);
plot(Dose,Y(:,18),'-sc','Linewidth',2,'markersize',8);
plot(Dose,Y(:,19),'--db','Linewidth',2,'markersize',8);
plot(Dose,Y(:,20),':pg','Linewidth',2,'markersize',8);
plot(Dose,Y(:,21),'-vb','Linewidth',2,'markersize',8);
hold off;
legend('P1-P1','P1-P2','P1-P3','P1-P4','P2-P2','P2-P3','P2-P4','P3-P3');
xlabel('Dose (kGy)');
ylabel('Yield (mol/L)')

figure(5);
plot(Dose,Y(:,7),':pg','Linewidth',2,'markersize',8,'MarkerEdgeColor','g');
hold on
plot(Dose,Y(:,8),'-.hr','Linewidth',2,'markersize',8);
plot(Dose,Y(:,9),'-+m','Linewidth',2,'markersize',8);

hold off
legend('CH4-in-bulk','H2-in-bulk','C2H6-in-bulk');
xlabel('Dose (kGy)');
ylabel('Yield (mol/L)')

%%
%������ʵ�����ݶԱ�
load('crosslinking_exp.mat')%�ױ��Ͱ�ˮ��ý����ܶȣ��ų����
load('ini_condition')

%��ʵ�齻���ܶȱ仯�����㽻�������Ũ�ȱ仯��
chain_exp=[];
chain_exp2=[];
Mw=72;%��ѧ�ظ���Ԫ������g/mol
for i=1:size(ini_condition,2)
    if strcmp(ini_condition{1,i},'P0')
        P0=ini_condition{2,i};%����Ũ��
        break
    end
end
M_distdata=[];
chain_exp(:,1)=crosslinking_exp(:,1);
chain_exp2(:,1)=crosslinking_exp(:,1);
for i=1:size(crosslinking_exp,1)
    for j=2:size(crosslinking_exp,2)
    Mc=Mw/crosslinking_exp(i,j);
        var=Mc/10;
        point=1000;
        dd=(Mc*5-200)/(point-1);
        Mc_dist=200:dd:Mc*5;
        Mc_dist=Mc_dist';
        M_dist=normpdf(Mc_dist,Mc,var);%�������γ��ָ�˹�ֲ�/��̬�ֲ�
        for m=1:length(M_dist)
            if M_dist(m)<1e-20
                M_dist(m)=0;
            end
        end
        M_distdata(:,size(M_distdata,2)+1:size(M_distdata,2)+2)=[Mc_dist,M_dist];
        %figure(9)
        %plot(Mc_dist,M_dist,'Linewidth',3)
        chain_exp2(i,j)=0;
        for k=2:length(Mc_dist)
            chain_exp2(i,j)=chain_exp2(i,j)+((Mc_dist(k)-Mc_dist(k-1))*M_dist(k))*(Mw/Mc_dist(k))*P0;%���Ǹ�˹�ֲ����������Ũ��
        end
        chain_exp2(i,j)=chain_exp2(i,j)-chain_exp2(1,j);
        chain_exp(i,j)=(crosslinking_exp(i,j)-crosslinking_exp(1,j))*P0;%��ֱֵ�Ӽ��������Ũ��mol/L
    end
end

%��ģ�⽻���ܶȱ仯�����㽻�������Ũ�ȱ仯��

Mc0=14400;%ƽ����ʼ�������������mol/g
C_H0=P0*Mw/Mc0;
C_Y0=0;
C_chain0=3*C_H0;
C_chaint=3*(C_H0+Y(:,10))+2*(C_Y0+Y(:,24));
Mct=Mc0*C_chain0./C_chaint;
Ne=25;%���᳤�ȣ���ѧ�ظ���Ԫ��Ŀ���߾ۺ϶�
point_entangle=floor(Mct/(Ne*Mw));%����������
for j=1:length(point_entangle)%ͳ�Ƽ�����������ֵ
    point_stat=0;
    if point_entangle(j)>=1
        for i=1:point_entangle(j)
            point_stat=point_stat+factorial(point_entangle(j))/(factorial(i-1)*factorial(point_entangle(j)-i));
        end
    end
    point_stat=(1/2)^point_entangle(j)*point_stat;
    point_entangle(j,2)=point_stat;%
end
%�������¹�ʽ���㣺delt=4X_Y+6X_H-6S
chain_simu=[];
chain_simu=4*Y(:,24)+6*Y(:,10)-6*Y(:,11);%δУ��������Ӱ��
for i=1:length(point_entangle)
    chain_simu2=4*Y(:,24)*(1+point_entangle(i,2))+6*Y(:,10)*(1+point_entangle(i,2))-6*Y(:,11);
end
figure(6);
plot(Dose,chain_simu2(:,1),'-b','Linewidth',3,'markersize',8,'MarkerEdgeColor','g');
hold on
plot(chain_exp2(:,1),chain_exp2(:,2),'hr','Linewidth',2,'markersize',10);
plot(chain_exp2(:,1),chain_exp2(:,3),'ok','Linewidth',2,'markersize',10);
plot(chain_exp2(:,1),chain_exp2(:,4),'*g','Linewidth',2,'markersize',10);
hold off
%xlim([-20 max(Dose)+20])
chainname=strcat('Chain-simu-',num2str(T),'K');
legend(chainname,'Chain-exp293K','Chain-exp313K','Chain-exp343K');
xlabel('Dose (kGy)');
ylabel('Changes of chain segments (mol/L)')

%%
%ģ�����������ʵ�����Ա�
load('gas_exp.mat')%�����������20,40��70�����������ʵ������
figure(7);
plot(Dose,Y(:,27),'-b','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��CH4
hold on
plot(gas_exp(:,1),gas_exp(:,2),'hr','Linewidth',2,'markersize',10);
plot(gas_exp(:,1),gas_exp(:,3),'ok','Linewidth',2,'markersize',10);
plot(gas_exp(:,1),gas_exp(:,4),'*g','Linewidth',2,'markersize',10);
hold off
% xlim([-20 max(Dose)+20])
ch4name=strcat('CH4-simu-',num2str(T),'K');
legend(ch4name,'CH4-exp293K','CH4-exp313K','CH4-exp343K');
xlabel('Dose (kGy)');
ylabel('CH4 (mol/L)')

figure(8);
plot(Dose,Y(:,28),'-b','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��H2
hold on
plot(gas_exp(:,1),gas_exp(:,5),'hr','Linewidth',2,'markersize',10);
plot(gas_exp(:,1),gas_exp(:,6),'ok','Linewidth',2,'markersize',10);
plot(gas_exp(:,1),gas_exp(:,7),'*g','Linewidth',2,'markersize',10);
hold off
% xlim([-20 max(Dose)+20])
h2name=strcat('H2-simu-',num2str(T),'K');
legend(h2name,'H2-exp293K','H2-exp313K','H2-exp343K');
xlabel('Dose (kGy)');
ylabel('H2 (mol/L)')

%%
%��ͬ��������CH4�����ʵ���ģ��Ա�
if T==298
    if I==1e-2 || I==1e-3
        load('CH4_exp.mat')%�����������20,40��70�����������ʵ������
        if I==1e-2
            if exist('CH4_simu2','var') && length(CH4_simu2)~=length(Dose)
                CH4_simu2=[];
            end
            CH4_simu2(:,1)=Dose(:,1);
            CH4_simu2(:,2)=Y(:,27);
        elseif I==1e-3
            if exist('CH4_simu3','var') && length(CH4_simu3)~=length(Dose)
                CH4_simu3=[];
            end
            CH4_simu3(:,1)=Dose(:,1);
            CH4_simu3(:,2)=Y(:,27);
        end
    end
end

if exist('CH4_simu2','var') || exist('CH4_simu3','var')
    figure(9);
    if exist('CH4_simu2','var') && ~exist('CH4_simu3','var')
        plot(CH4_simu2(:,1),CH4_simu2(:,2),'-b','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��I=1e-2��CH4
        hold on
        plot(CH4_exp(:,1),CH4_exp(:,2),'hg','Linewidth',2,'markersize',10);%ʵ��ģ��I=1e-2��CH4
        hold off
%         xlim([-0.2 1.2])
        if T==298
            ch4name2=strcat('CH4-simu-1E-2-',num2str(T),'K');
        end
        legend(ch4name2,'CH4-exp-1E-2-298K');
        xlabel('Dose (kGy)');
        ylabel('CH4@1E-2 (mol/L)')
    elseif exist('CH4_simu3','var') && ~exist('CH4_simu2','var')
        plot(CH4_simu3(:,1),CH4_simu3(:,2),'-r','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��I=1e-3��CH4
        hold on
        plot(CH4_exp(:,1),CH4_exp(:,3),'ok','Linewidth',2,'markersize',10);%ʵ��ģ��I=1e-3��CH4
        hold off
%         xlim([-0.2 1.2])
        if T==298
            ch4name3=strcat('CH4-simu-1E-3-',num2str(T),'K');
        end
        legend(ch4name3,'CH4-exp-1E-3-298K');
        xlabel('Dose (kGy)');
        ylabel('CH4@1E-3 (mol/L)')
    else
        plot(CH4_simu2(:,1),CH4_simu2(:,2),'-b','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��I=1e-2��CH4
        hold on
        plot(CH4_simu3(:,1),CH4_simu3(:,2),'-r','Linewidth',4,'markersize',8,'MarkerEdgeColor','g');%ģ��I=1e-3��CH4
        plot(CH4_exp(:,1),CH4_exp(:,2),'hg','Linewidth',2,'markersize',10);%ʵ��ģ��I=1e-2��CH4
        plot(CH4_exp(:,1),CH4_exp(:,3),'ok','Linewidth',2,'markersize',10);%ʵ��ģ��I=1e-3��CH4
        hold off
%         xlim([-0.2 1.2])
        if T==298
            ch4name2=strcat('CH4-simu-1E-2-',num2str(T),'K');
            ch4name3=strcat('CH4-simu-1E-3-',num2str(T),'K');
        end
        legend(ch4name2,ch4name3,'CH4-exp-1E-2-298K','CH4-exp-1E-3-298K');
        xlabel('Dose (kGy)');
        ylabel('CH4 (mol/L)')
    end
end
%%
clear i eventans j options Dosemax Dosemax2 h2name I T tmax tmaxgamma dd k point C_H0 C_Y0 m Ne 

