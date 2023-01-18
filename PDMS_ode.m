%Dynamics differential functions for the radiolysis of PDMS.
function dy=PDMS_ode(t,y)
dy=zeros(29,1);%΢�ַ��̸���

%%
%��������Ԥ��,�����ʡ��¶ȡ�Ӧ����ʪ�ȡ�����
% I=1e-3;% Gy/s
% T=293;% k
% tmax=1e6;%s,����ʱ�䣬������պ�Ӧ
global I T tmax on_off1 on_off2 on_off3 on_off4 on_off5 on_off6 on_off7 on_off8 on_off8 on_off10 on_off11 on_off12 on_off13
p0=9.72;%��ʼ���嵥ԪŨ��mol/L�����Ǿۺ�������ϱ���
IP=9.6;%�����ܣ�eV
density=0.99;%����ĭ�����ܶȣ�g/cm^3����ĭ0.54g/cm3
porosity_stress=0.54;%ʵ��������0.47���ܶ��Ƶ���0.58
fai=1;%����Ч�ʣ���qy=1-[1-exp(-BI)]^n
R_yield=7.2*fai;%���ɻ���ѧ����,��/100eV���߷��ӷ��仯ѧ�飬δ��������Ч�ʵ�Ӱ�����أ�Ĭ��ֵ7.2

%%
%298k�¶��¾ۺ��������̬���������ܶȡ��ݿ��ʡ����ϡ��ᾧ�ȡ�ȡ����ɢ��͸ϵ����������������ѧ���������,ʱ���¶�������
%ˮ�ܽ��ϵ����ˮ�����ܺͷ���

%������ת��Ϊ��ѧ��
R_gama_sink=I/(1.602e-19)/density;%������������eV/L.s
P0act_gama_ini=R_gama_sink/IP/(6.02e23);%����̬��������mol/Ls�����뼤��ʱ��1e-17~1e-16 s,�߼���̬��ԥʱ��1e-15~1e-13 s
R_gen=R_gama_sink/100*R_yield/(6.02e23);%���ɻ��������ʣ�mol/Ls�����������ӷֽ�ʱ��1e-14 s
quantum_fai=R_gen/P0act_gama_ini/2;%����Ч�ʣ�һ���������ӷֽ���������ɻ������߾���ʽ��fai=1-(1-exp(-1e-9*I))^0.003;%��������Ч��qy=[1-exp(-BI)]^n��������ЧӦ���ֲ��Ȼ���Ӱ��
if on_off12==1
    fprintf('\n298Kʱ���ۼ���̬���ӷֽ������Ч��Ϊ%0.4f',quantum_fai);
    on_off12=on_off12+1;
end

%����ЧӦ
ratio_P1H_cage=0.6;
ratio_P2CH3_cage=0.5;
ratio_P3P4_cage=0.4;
%��Ӧ��֧��
branch_P1H=0.20;
branch_P2CH3=0.30;
branch_P3P4=0.50;
%���ֽⷴӦ���µĳ�ʼ���ɻ���֧��
ratio_P1H=branch_P1H*ratio_P1H_cage;
ratio_P2CH3=branch_P2CH3*ratio_P2CH3_cage;
ratio_P3P4=branch_P3P4*ratio_P3P4_cage;
quantum_fai2=quantum_fai*(ratio_P1H+ratio_P2CH3+ratio_P3P4);
if on_off11==1
    fprintf('\n298Kʱ��������ЧӦ����̬���ӷֽ������Ч��Ϊ%0.4f',quantum_fai2);
    on_off11=on_off11+1;
end

%%
%298k�¶��¶���ѧ����,�����¶ȡ�Ӧ������ɢ���ƹ���,������/s��˫����L/mol.s
k14=2e6;%����̬�˼�����,k11-k13�ں���
k21=1.45E13;k22=5.82E9;k23=5.58E8;k24=1.94E-10;k25=2.59E10;k26=7.18E-2;k27=2.93E-24;k28=5.65E2;k29=3.37E8;k210=2.74E-26;k211=1.24E-6;k212=1.04E9;%k25�Ż�
k31=7.53E10;k32=1.67E10;k33=1.69e10;k34=2.65e9;k35=3.12e11;k36=6.68e11;k37=2.06e8;k38=9.75e10;k39=7.10e10;k310=5.56E-2;k311=3.73E-5;k312=5.18E-6;k313=2.22E5;
%k31=1E7;k32=1.67E7;k33=1.69e5;k34=2.65e6;k35=3.12e8;k36=6.68e8;k37=2.06e5;k38=9.75e7;k39=7.10e7;
% k41=1.70e12;k42=5e7;k43=3.24e11;k44=1e4;k45=1.14e12;k46=3.15e10;k47=1e4;k48=6.44e10;
k41=0;k42=0;k43=0;k44=0;k45=0;k46=0;k47=0;k48=0;
k51=0;k52=4.13E7;k53=5.40E1;
%k51=2.68E10;k52=4.13E7;k53=5.40E1;

%��������˹��ʽ�¶Ȼ��J/mol���������չ������ʽ
E11=5000;E12=5000;E13=10000;
E21=600;E22=-1740;E23=5110;E24=88450;E25=750;E26=32920;E27=155630;E28=5320;E29=2860;E210=170710;E211=65170;E212=37170;
E31=-2910;E32=-2170;E33=-760;E34=-2970;E35=-850;E36=-1310;E37=-2000;E38=-1360;E39=-2260;E310=-970;E311=27000;E312=34050;E313=890;
E41=-810;E42=15000;E43=-2180;E44=35000;E45=-410;E46=-2730;E47=5000;E48=-2080;
E51=18750;E52=3940;E53=4860;

%�����¶�У����Ӧ����,ע����Խת���¶�
T_corr11=exp(-E11/8.314/T)/exp(-E11/8.314/298);
T_corr12=exp(-E12/8.314/T)/exp(-E12/8.314/298);
T_corr13=exp(-E13/8.314/T)/exp(-E13/8.314/298);
T_corr21=(T/298)^2.12*exp(-E21/8.314/T)/exp(-E21/8.314/298);
T_corr22=(T/298)^0.69*exp(-E22/8.314/T)/exp(-E22/8.314/298);
T_corr23=(T/298)^1.59*exp(-E23/8.314/T)/exp(-E23/8.314/298);
T_corr24=exp(-E24/8.314/T)/exp(-E24/8.314/298);
T_corr25=(T/298)^-0.46*exp(-E25/8.314/T)/exp(-E25/8.314/298);
T_corr26=exp(-E26/8.314/T)/exp(-E26/8.314/298);
T_corr27=exp(-E27/8.314/T)/exp(-E27/8.314/298);
T_corr28=(T/298)^5.8*exp(-E28/8.314/T)/exp(-E28/8.314/298);
T_corr29=exp(-E29/8.314/T)/exp(-E29/8.314/298);
T_corr210=(T/298)^4.82*exp(-E210/8.314/T)/exp(-E210/8.314/298);
T_corr211=(T/298)^0.14*exp(-E211/8.314/T)/exp(-E211/8.314/298);
T_corr212=(T/298)^1.89*exp(-E212/8.314/T)/exp(-E212/8.314/298);
T_corr31=(T/298)^-1.70*exp(-E31/8.314/T)/exp(-E31/8.314/298);
T_corr32=(T/298)^-1.84*exp(-E32/8.314/T)/exp(-E32/8.314/298);
T_corr33=(T/298)^0.08*exp(-E33/8.314/T)/exp(-E33/8.314/298);
T_corr34=(T/298)^-0.68*exp(-E34/8.314/T)/exp(-E34/8.314/298);
T_corr35=(T/298)^-2.87*exp(-E35/8.314/T)/exp(-E35/8.314/298);
T_corr36=(T/298)^-2.53*exp(-E36/8.314/T)/exp(-E36/8.314/298);
T_corr37=(T/298)^-0.91*exp(-E37/8.314/T)/exp(-E37/8.314/298);
T_corr38=(T/298)^-2.51*exp(-E38/8.314/T)/exp(-E38/8.314/298);
T_corr39=(T/298)^-1.92*exp(-E39/8.314/T)/exp(-E39/8.314/298);
T_corr310=(T/298)^15.50*exp(-E310/8.314/T)/exp(-E310/8.314/298);
T_corr311=(T/298)^8.35*exp(-E311/8.314/T)/exp(-E311/8.314/298);
T_corr312=(T/298)^5.48*exp(-E312/8.314/T)/exp(-E312/8.314/298);
T_corr313=(T/298)^2.05*exp(-E313/8.314/T)/exp(-E313/8.314/298);
T_corr41=(T/298)^-1.39*exp(-E41/8.314/T)/exp(-E41/8.314/298);
T_corr42=exp(-E42/8.314/T)/exp(-E42/8.314/298);
T_corr43=(T/298)^-1.95*exp(-E43/8.314/T)/exp(-E43/8.314/298);
T_corr44=exp(-E44/8.314/T)/exp(-E44/8.314/298);
T_corr45=(T/298)^-1.65*exp(-E45/8.314/T)/exp(-E45/8.314/298);
T_corr46=(T/298)^-1.55*exp(-E46/8.314/T)/exp(-E46/8.314/298);
T_corr47=exp(-E47/8.314/T)/exp(-E47/8.314/298);
T_corr48=(T/298)^2.00*exp(-E48/8.314/T)/exp(-E48/8.314/298);
T_corr51=(T/298)^1.04*exp(-E51/8.314/T)/exp(-E51/8.314/298);
T_corr52=(T/298)^2.74*exp(-E52/8.314/T)/exp(-E52/8.314/298);
T_corr53=(T/298)^0.00*exp(-E53/8.314/T)/exp(-E53/8.314/298);
%%
%298k�ۺ����������ɢϵ��ָǰ����m^2/s����ɢ���J/mol
load('ini_dynamics_bulk.mat');
%298k�ۺ�������϶����ɢϵ��ָǰ����m^2/s����ɢ���J/mol
load('ini_dynamics_void.mat');
%����298k�¶��µ���ɢϵ��
if size(ini_dynamics_bulk,1)~=size(ini_dynamics_void,1)
    error('����ini_dynamics_bulk��ini_dynamics_void������Ŀ�Ƿ���ͬ');
end
for i=1:size(ini_dynamics_bulk,1)
    ini_dynamics_bulk{i,5}=ini_dynamics_bulk{i,2}*exp(-ini_dynamics_bulk{i,4}/298/8.314);
    ini_dynamics_void{i,5}=ini_dynamics_void{i,2}*exp(-ini_dynamics_void{i,4}/298/8.314);
end

%298k�¶��������ݿ��ʶԱ���ۺ���ɢ����D=Dp+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(Dh-Dp)
for i=1:size(ini_dynamics_bulk,1)
    ini_dynamics_bulk{i,6}=ini_dynamics_bulk{i,5}+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(ini_dynamics_void{i,5}-ini_dynamics_bulk{i,5});
    ini_dynamics_void{i,6}=ini_dynamics_bulk{i,5}+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(ini_dynamics_void{i,5}-ini_dynamics_bulk{i,5});
end

%�ۺ�������кͿ�϶��ָ���¶��µ���ɢϵ��
for i=1:size(ini_dynamics_bulk,1)
    ini_dynamics_bulk{i,7}=ini_dynamics_bulk{i,2}*exp(-ini_dynamics_bulk{i,4}/T/8.314);
    ini_dynamics_void{i,7}=ini_dynamics_void{i,2}*exp(-ini_dynamics_void{i,4}/T/8.314);
end

%ָ���¶������±���ۺ���ɢ����D=Dp+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(Dh-Dp)
for i=1:size(ini_dynamics_bulk,1)
    ini_dynamics_bulk{i,8}=ini_dynamics_bulk{i,7}+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(ini_dynamics_void{i,7}-ini_dynamics_bulk{i,7});
    ini_dynamics_void{i,8}=ini_dynamics_bulk{i,7}+(1.3*porosity_stress^2-0.3*porosity_stress^3)*(ini_dynamics_void{i,7}-ini_dynamics_bulk{i,7});
end

%���ڻ�Ϸ����������ɢϵ��D=(D0*D1)^0.5����У����Ӧ����

D_cor11=(ini_dynamics_bulk{2,8}*ini_dynamics_bulk{2,8})^0.5/(ini_dynamics_bulk{2,6}*ini_dynamics_bulk{2,6})^0.5;
D_cor12=D_cor11;
D_cor13=D_cor11;
D_cor21=(ini_dynamics_bulk{8,8}*ini_dynamics_bulk{8,8})^0.5/(ini_dynamics_bulk{8,6}*ini_dynamics_bulk{8,6})^0.5;
D298_22=(ini_dynamics_bulk{8,6}*ini_dynamics_bulk{9,6})^0.5;D_22=(ini_dynamics_bulk{8,8}*ini_dynamics_bulk{9,8})^0.5;D_cor22=D_22/D298_22;
D298_23=(ini_dynamics_bulk{1,6}*ini_dynamics_bulk{8,6})^0.5;D_23=(ini_dynamics_bulk{1,8}*ini_dynamics_bulk{8,8})^0.5;D_cor23=D_23/D298_23;
D298_24=(ini_dynamics_bulk{7,6}*ini_dynamics_bulk{8,6})^0.5;D_24=(ini_dynamics_bulk{7,8}*ini_dynamics_bulk{8,8})^0.5;D_cor24=D_24/D298_24;
D_cor25=(ini_dynamics_bulk{9,8}*ini_dynamics_bulk{9,8})^0.5/(ini_dynamics_bulk{9,6}*ini_dynamics_bulk{9,6})^0.5;
D298_26=(ini_dynamics_bulk{7,6}*ini_dynamics_bulk{9,6})^0.5;D_26=(ini_dynamics_bulk{7,8}*ini_dynamics_bulk{9,8})^0.5;D_cor26=D_26/D298_26;
D_cor27=D_cor26;
D298_28=(ini_dynamics_bulk{9,6}*ini_dynamics_bulk{13,6})^0.5;D_28=(ini_dynamics_bulk{9,8}*ini_dynamics_bulk{13,8})^0.5;D_cor28=D_28/D298_28;
D298_29=(ini_dynamics_bulk{8,6}*ini_dynamics_bulk{13,6})^0.5;D_29=(ini_dynamics_bulk{8,8}*ini_dynamics_bulk{13,8})^0.5;D_cor29=D_29/D298_29;
D298_210=(ini_dynamics_bulk{9,6}*ini_dynamics_bulk{14,6})^0.5;D_210=(ini_dynamics_bulk{9,8}*ini_dynamics_bulk{14,8})^0.5;D_cor210=D_210/D298_210;
D_cor211=D_cor210;
D298_212=(ini_dynamics_bulk{8,6}*ini_dynamics_bulk{14,6})^0.5;D_212=(ini_dynamics_bulk{8,8}*ini_dynamics_bulk{14,8})^0.5;D_cor212=D_210/D298_212;
D_cor31=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{3,8})^0.5/(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{3,6})^0.5;
D298_32=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{4,6})^0.5;D_32=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{4,8})^0.5;D_cor32=D_32/D298_32;
D298_33=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{5,6})^0.5;D_33=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{5,8})^0.5;D_cor33=D_33/D298_33;
D298_34=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{6,6})^0.5;D_34=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{6,8})^0.5;D_cor34=D_34/D298_34;
D_cor35=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{4,8})^0.5/(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{4,6})^0.5;
D298_36=(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{5,6})^0.5;D_36=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{5,8})^0.5;D_cor36=D_36/D298_36;
D298_37=(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{6,6})^0.5;D_37=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{6,8})^0.5;D_cor37=D_37/D298_37;
D_cor38=(ini_dynamics_bulk{5,8}*ini_dynamics_bulk{5,8})^0.5/(ini_dynamics_bulk{5,6}*ini_dynamics_bulk{5,6})^0.5;
D298_39=(ini_dynamics_bulk{5,6}*ini_dynamics_bulk{6,6})^0.5;D_39=(ini_dynamics_bulk{5,8}*ini_dynamics_bulk{6,8})^0.5;D_cor39=D_39/D298_39;
D298_310=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{7,6})^0.5;D_310=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{7,8})^0.5;D_cor310=D_310/D298_310;
D298_311=(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{7,6})^0.5;D_311=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{7,8})^0.5;D_cor311=D_311/D298_311;
D298_312=(ini_dynamics_bulk{5,6}*ini_dynamics_bulk{7,6})^0.5;D_312=(ini_dynamics_bulk{5,8}*ini_dynamics_bulk{7,8})^0.5;D_cor312=D_312/D298_312;
D298_313=(ini_dynamics_bulk{6,6}*ini_dynamics_bulk{7,6})^0.5;D_313=(ini_dynamics_bulk{6,8}*ini_dynamics_bulk{7,8})^0.5;D_cor313=D_313/D298_313;
D298_41=(ini_dynamics_bulk{5,6}*ini_dynamics_bulk{8,6})^0.5;D_41=(ini_dynamics_bulk{5,8}*ini_dynamics_bulk{8,8})^0.5;D_cor41=D_41/D298_41;
D298_42=(ini_dynamics_bulk{6,6}*ini_dynamics_bulk{8,6})^0.5;D_42=(ini_dynamics_bulk{6,8}*ini_dynamics_bulk{8,8})^0.5;D_cor42=D_42/D298_42;
D298_43=(ini_dynamics_bulk{5,6}*ini_dynamics_bulk{9,6})^0.5;D_43=(ini_dynamics_bulk{5,8}*ini_dynamics_bulk{9,8})^0.5;D_cor43=D_43/D298_43;
D298_44=(ini_dynamics_bulk{6,6}*ini_dynamics_bulk{9,6})^0.5;D_44=(ini_dynamics_bulk{6,8}*ini_dynamics_bulk{9,8})^0.5;D_cor44=D_44/D298_44;
D298_45=(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{8,6})^0.5;D_45=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{8,8})^0.5;D_cor45=D_45/D298_45;
D298_46=(ini_dynamics_bulk{4,6}*ini_dynamics_bulk{9,6})^0.5;D_46=(ini_dynamics_bulk{4,8}*ini_dynamics_bulk{9,8})^0.5;D_cor46=D_46/D298_46;
D298_47=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{8,6})^0.5;D_47=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{8,8})^0.5;D_cor47=D_47/D298_47;
D298_48=(ini_dynamics_bulk{3,6}*ini_dynamics_bulk{9,6})^0.5;D_48=(ini_dynamics_bulk{3,8}*ini_dynamics_bulk{9,8})^0.5;D_cor48=D_48/D298_48;
D298_51=(ini_dynamics_bulk{9,6}*ini_dynamics_bulk{9,6})^0.5;D_51=(ini_dynamics_bulk{9,8}*ini_dynamics_bulk{9,8})^0.5;D_cor51=D_51/D298_51;
D298_52=(ini_dynamics_bulk{9,6}*ini_dynamics_bulk{12,6})^0.5;D_52=(ini_dynamics_bulk{9,8}*ini_dynamics_bulk{12,8})^0.5;D_cor52=D_52/D298_52;
D298_53=(ini_dynamics_bulk{16,6}*ini_dynamics_bulk{8,6})^0.5;D_53=(ini_dynamics_bulk{16,8}*ini_dynamics_bulk{8,8})^0.5;D_cor53=D_53/D298_53;
%D_corqy=(D_cor23*D_cor26*D_cor39)^(1/3);%�������������Ч�ʵ���ɢ����У��

%%
%٤��������������������̬���ٶȣ��������¶��޹أ��о����ִ��ں�С�Ļ�ܣ�Frank�CRabinowitch��s principle (cagee ffect).����һ����5-20 kJ/mol
k10=P0act_gama_ini/p0;%����̬�������ʳ���,s-1
if t>tmax %�رշ��պ�Ӧ
    P0act_gama_ini=0;
    k10=0;
    if on_off1==1
        fprintf('\nģ��ʱ�䳬�����տ���ʱ�䣬���ڷ����Ӧ����ʱ P0act_gama_ini=%f\n', P0act_gama_ini);
        on_off1=on_off1+1;
    end
end
%�����¶Ȼ�ܡ���ɢ���ƵȵĶ���ѧ����

P0act_gama_ini=k10*y(12);%����̬���Ӳ�������
R_gen_all=P0act_gama_ini*quantum_fai;%298K�������ɻ����������ʣ��ɷֽ�����Ч��ȷ��
K11=R_gen_all*ratio_P1H*D_cor11*T_corr11;
K12=R_gen_all*ratio_P2CH3*D_cor12*T_corr12;
K13=R_gen_all*ratio_P3P4*D_cor13*T_corr13;
k14=k14*D_cor35*T_corr35;
R_gen_all=K11+K12+K13;
if P0act_gama_ini>0
    quantum_fai_now=R_gen_all/P0act_gama_ini;
    if on_off2==1
        fprintf('\n%dK�¶��¿�������ЧӦ����ɢЧӦ���¶�ЧӦ����̬���ӷֽ������Ч��Ϊ%0.4f',T,quantum_fai_now);
        on_off2=on_off2+1;
    end
end
if t>tmax %�رշ��պ�Ӧ
    R_gen_all=0;
    K11=R_gen_all*ratio_P1H*D_cor11*T_corr11;K12=R_gen_all*ratio_P2CH3*D_cor12*T_corr12;K13=R_gen_all*ratio_P3P4*D_cor13*T_corr13;k14=k14*D_cor35*T_corr35;
    if on_off3==1
        fprintf('\nģ��ʱ�䳬�����տ���ʱ�䣬���ڷ����Ӧ����ʱ R_gen_all=%f\n', R_gen_all);
        on_off3=on_off3+1;
    end
end
k21=k21*D_cor21*T_corr21;k22=k22*D_cor22*T_corr22;k23=k23*D_cor23*T_corr23;k24=k24*D_cor24*T_corr24;k25=k25*D_cor25*T_corr25;k26=k26*D_cor26*T_corr26;k27=k27*D_cor27*T_corr27;k28=k28*D_cor28*T_corr28;k29=k29*D_cor29*T_corr29;k210=k210*D_cor210*T_corr210;k211=k211*D_cor211*T_corr211;k212=k212*D_cor212*T_corr212;
k31=k31*D_cor31*T_corr31;k32=k32*D_cor32*T_corr32;k33=k33*D_cor33*T_corr33;k34=k34*D_cor34*T_corr34;k35=k35*D_cor35*T_corr35;k36=k36*D_cor36*T_corr36;k37=k37*D_cor37*T_corr37;k38=k38*D_cor38*T_corr38;k39=k39*D_cor39*T_corr39;k310=k310*D_cor310*T_corr310;k311=k311*D_cor311*T_corr311;k312=k312*D_cor312*T_corr312;k313=k313*D_cor313*T_corr313;
k41=k41*D_cor41*T_corr41;k42=k42*D_cor42*T_corr42;k43=k43*D_cor43*T_corr43;k44=k44*D_cor44*T_corr44;k45=k45*D_cor45*T_corr45;k46=k46*D_cor46*T_corr46;k47=k47*D_cor47*T_corr47;k48=k48*D_cor48*T_corr48;
k51=k51*D_cor51*T_corr51;k52=k52*D_cor52*T_corr52;k53=k53*D_cor53*T_corr53;
%%
for i=1:29
    if y(i)<=0
        y(i)=0;
    end
end
if y(12) <=0
    fprintf('\nע�⣺�ڸ����Ļ������²�����ȫ���⣡����\n');
end
if y(13)<=0 && on_off13==1
    fprintf('\n����̬���Ӻľ���\n');
    on_off13=on_off13+1;
end

%%
%���Է�Χ����ѧ΢�ַ���
dy(1)=K11-2*k21*y(1)^2+k210*y(2)*y(26)-k22*y(1)*y(2)-(k23+k24)*y(1)*y(12)-k41*y(5)*y(1)-k42*y(6)*y(1)-k45*y(4)*y(1)-k47*y(1)*y(3)-k29*y(1)*y(25)+k51*y(7)*y(2)+k52*y(8)*y(2)-k53*y(9)*y(1)-k212*y(1)*y(26);%H
dy(2)=K12-k22*y(1)*y(2)-2*k25*y(2)^2-(k26+k27)*y(2)*y(12)-k43*y(5)*y(2)-k44*y(6)*y(2)-k46*y(2)*(4)-k48*y(2)*y(3)-k28*y(2)*y(25)-k210*y(2)*y(26)-k211*y(2)*y(26)-k51*y(7)*y(2)-k52*y(8)*y(2)+k53*y(9)*y(1);%CH3
dy(3)=K11+k23*y(1)*y(12)+k26*y(2)*y(12)-2*k31*y(3)^2-k32*y(3)*y(4)-k33*y(3)*y(5)-k34*y(3)*y(6)-k47*y(1)*y(3)-k48*y(2)*y(3)+k311*y(4)*y(23)+k312*y(5)*y(23)+k313*y(6)*y(23);%p1
dy(4)=K12+k24*y(1)*y(12)+k27*y(2)*y(12)-k32*y(3)*y(4)-2*k35*y(4)^2-k36*y(4)*y(5)-k37*y(4)*y(6)-k45*y(4)*y(1)-k46*y(2)*(4)-k311*y(4)*y(23);%p2
dy(5)=K13-k36*y(4)*y(5)-2*k38*y(5)^2-k39*y(5)*y(6)-k41*y(5)*y(1)-k43*y(2)*y(5)-k33*y(3)*y(5)+k28*y(2)*y(25)+k29*y(1)*y(25)-k312*y(5)*y(23);%p3
dy(6)=K13-k37*y(4)*y(5)-k39*y(5)*y(6)-k42*y(6)*y(1)-k44*y(6)*y(2)-k34*y(3)*y(6)+k211*y(2)*y(26)+k212*y(1)*y(26)-k313*y(6)*y(23);%p4

if y(7)>=1.95E-2%������CH4����
    y(7)=1.95E-2;
    dy(7)=0;
    if on_off4==1
        fprintf('\n������CH4���ͣ�\n');
        on_off4=on_off4+1;
    end
end
if y(8)>=2.65E-3%������H2����
    y(8)=2.65E-3;
    dy(8)=0;
    if on_off5==1
        fprintf('\n������H2���ͣ�\n');
        on_off5=on_off5+1;
    end
end
if y(9)>=1.12E-1%������C2H6����
    y(9)=1.12E-1;
    dy(9)=0;
    if on_off6==1
        fprintf('\n������C2H6���ͣ�\n');
        on_off6=on_off6+1;
    end
end

dy(7)=k22*y(1)*y(2)+k24*y(1)*y(12)+k26*y(2)*y(12)+k28*y(2)*y(25)+k211*y(2)*y(26)-k51*y(7)*y(2)+k52*y(8)*y(2)+k53*y(9)*y(1);%������CH4����
if y(7)>=1.95E-2%������CH4����
    y(7)=1.95E-2;
    dy(7)=0;
    if on_off7==1
        fprintf('\n������CH4���ͣ�');
        on_off7=on_off7+1;
    end
end
dy(8)=k21*y(1)^2+k23*y(1)*y(12)+k29*y(1)*y(25)-k52*y(8)*y(2)+k212*y(1)*y(26);%������H2����
if y(8)>=2.65E-3%������H2����
    y(8)=2.65E-3;
    dy(8)=0;
    if on_off8==1
        fprintf('\n������H2���ͣ�');
    end
    on_off8=on_off8+1;
end
dy(9)=k25*y(2)^2+k27*y(2)*y(12)+k51*y(7)*y(2)-k53*y(9)*y(1);%������C2H6����
if y(9)>=1.12E-1%������C2H6����
    y(9)=1.12E-1;
    dy(9)=0;
    if on_off9==1
        fprintf('\n������C2H6���ͣ�');
        on_off9=on_off9+1;
    end
end
if on_off10==1
    fprintf('\nH2������%f',y(8))
    on_off10=on_off10+1;
end

dy(10)=k31*y(3)^2+k32*y(3)*y(4)+k35*y(4)^2;%H�ͽ�����
dy(11)=K13-k39*y(5)*y(6);%������
dy(12)=-P0act_gama_ini-k23*y(1)*y(23)-k24*y(1)*y(23)-k26*y(2)*(23)-k27*y(2)*(23)+k39*y(5)*y(6)+k46*y(2)*(4)+k47*y(1)*y(3)+k14*y(13)-k311*y(4)*y(23)-k312*y(5)*y(23)-k313*y(6)*y(23);%+(k31*y(3)^2)*2/3+(k32*y(3)*y(4))*2/3+(k33*y(3)*y(5))*2/3+(k34*y(3)*y(6))*2/3+(k35*y(4)^2)*2/3+(k36*y(4)*y(5))*2/3+(k37*y(4)*y(6))*2/3+(k38*y(5)^2)+(k39*y(5)*y(6))+k41*y(5)*y(1)*1/2+dy(25)+k42*y(6)*y(1)+k43*y(2)*y(5)*3/2+k44*y(6)*y(2)*3/2+k45*y(4)*y(1)*1/2+k46*y(2)*(4)+k47*y(1)*y(3)+k48*y(2)*y(3)*3/2+k14*y(13);%��ЧP0���嵥Ԫ(P0-equivalent)�����ʱ��仯
dy(13)=P0act_gama_ini-R_gen_all-k14*y(13);%p0act����̬���������ʱ��仯
dy(14)=k31*y(3)^2;%p1-p1
dy(15)=k32*y(3)*y(4);%p1-p2
dy(16)=k33*y(3)*y(5);%p1-p3
dy(17)=k34*y(3)*y(6);%p1-p4
dy(18)=k35*y(4)^2;%p2-p2
dy(19)=k36*y(4)*y(5)-k10*y(19)*2/3;%p2-p3
dy(20)=k37*y(4)*y(6)-k10*y(20)*2/3;%p2-p4
dy(21)=k38*y(5)^2;%p3-p3
dy(22)=k39*y(5)*y(6)-K13;%p3-p4
dy(23)=-P0act_gama_ini-k23*y(1)*y(23)-k24*y(1)*y(23)-k26*y(2)*(23)-k27*y(2)*(23)+k39*y(5)*y(6)+k46*y(2)*(4)+k47*y(1)*y(3)+k14*y(13)-k311*y(4)*y(23)-k312*y(5)*y(23)-k313*y(6)*y(23);%��ʼ���嵥ԪP0�����ʱ��仯
dy(24)=k33*y(3)*y(5)+k34*y(3)*y(6)+k36*y(4)*y(5)+k37*y(4)*y(6);%Y�ͽ�����
dy(25)=k41*y(5)*y(1)+k45*y(4)*y(1)-k28*y(2)*y(25)-k29*y(1)*y(25)+k311*y(4)*y(23)+k312*y(5)*y(23);%����
dy(26)=k42*y(6)*y(1)-k210*y(2)*y(26)-k211*y(2)*y(26)-k212*y(1)*y(26)+k313*y(6)*y(23);%�贼
dy(27)=k22*y(1)*y(2)+k24*y(1)*y(12)+k26*y(2)*y(12)+k28*y(2)*y(25)+k211*y(2)*y(26)-k51*y(7)*y(2)+k52*y(8)*y(2)+k53*y(9)*y(1);%CH4��������
dy(28)=k21*y(1)^2+k23*y(1)*y(12)+k29*y(1)*y(25)-k52*y(8)*y(2);%H2��������
dy(29)=k25*y(2)^2+k27*y(2)*y(12)+k51*y(7)*y(2)-k53*y(9)*y(1);%C2H6����in bulk
%%
%����������ĺ���
%��ɢ����
%��ֽṹ����

end


