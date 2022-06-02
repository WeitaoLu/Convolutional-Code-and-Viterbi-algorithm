clear;
total=2000;%����֡��
SNR = -2:1:13; 
for n = 1:length(SNR)%��ͬSNR
    BER=0;%BERX�е�X�������õ�L0���
    BER2=0;
    BER4=0;
    BER8=0;
    BER10=0;
    BER25=0;
    BERadam=0;%����Ӧ�㷨��BER
    Ladam=0;%����Ӧ�㷨��ƽ��L
    %ע�Ͳ��ֵĴ���ΪMatlab�⺯��������д�ĺ���ִ�н����ͬ
    pb=0;
    for k=1:total%total��Ӧ�����֡����ÿһ֡���棩
    D=[0,0];%(D(i)Ϊ��i���Ĵ���)
    input_data = [randi([0,1],[1,498]),[0,0]]; %498λ�������+2λ������
    %trellis = poly2trellis(3,[7 5]);
    %code_data=convenc(input_data,trellis);
    code_data2=mycc2_1_2(input_data,D);%��������
    [QPSK_data,~]=QPSK_Modulation(SNR(n),code_data2);%QPSK���ơ��ŵ����桢QPSK���
    [~,Pb]=QPSK_Modulation(SNR(n),input_data);%�޾����ֱ�ӵ���
    %tmp=issame(code_data,code_data2);
    %decode = vitdec(code_data2,trellis,15,'term','hard');
    %decode4=vitdec(QPSK_data,trellis,15,'term','hard');
    decode3=viterbi_hard(QPSK_data);%ʹ���Լ���Viterbi�㷨����
    decode2=viterbi_short_hard(QPSK_data,2);%ʹ���Լ���Viterbi�ض��㷨����
    decode4=viterbi_short_hard(QPSK_data,4);
    decode8=viterbi_short_hard(QPSK_data,8);
    decode10=viterbi_short_hard(QPSK_data,10);
    decode25=viterbi_short_hard(QPSK_data,25);
    [decode_adam,L]=adam_viterbi(QPSK_data);%ʹ���Լ�������ӦViterbi�ض��㷨����
    %correct=issame(decode,input_data);
    %correct2=issame(decode2,input_data);
    correct3=issame(decode3,input_data);
    correct2=issame(decode2,input_data);
    pb=pb+Pb;
    BER=BER+length(input_data(:))-correct3;%�ֱ�������еĽ�������ӻ�
    BER2=BER2+length(input_data(:))-correct2;
    BER4=BER4+length(input_data(:))-issame(decode4,input_data);
    BER8=BER8+length(input_data(:))-issame(decode8,input_data);
    BER10=BER10+length(input_data(:))-issame(decode10,input_data);
    BER25=BER25+length(input_data(:))-issame(decode25,input_data);
    BERadam=BERadam+length(input_data(:))-issame(decode_adam,input_data);
    Ladam=Ladam+L;
    end
    BER=BER/(total*length(input_data(:)));
    BER2=BER2/(total*length(input_data(:)));%L0=2
    BER4=BER4/(total*length(input_data(:)));%L0=4
    BER8=BER8/(total*length(input_data(:)));%L0=8
    BER10=BER10/(total*length(input_data(:)));%L0=10
    BER25=BER25/(total*length(input_data(:)));%L0=25
    BERadam=BERadam/(total*length(input_data(:)));%L0=25
    pb=pb/total;
    Ladam=Ladam/total;
    pbrayleigh(n)=pb;
    BERrayleigh(n) = BER;
    BERrayleigh2(n) = BER2;
    BERrayleigh4(n) = BER4;
    BERrayleigh8(n) = BER8;
    BERrayleigh10(n) = BER10;
    BERrayleigh25(n) = BER25;
    BERrayleighadam(n) = BERadam;
    Ladamlist(n) = Ladam;
end
% Plot BER results.
semilogy(SNR,pbrayleigh,'-k*',SNR,BERrayleigh,'-r*',SNR,BERrayleigh2,'-*',SNR,BERrayleigh4,'-*',SNR,BERrayleigh8,'-*',SNR,BERrayleigh10,'-*',SNR,BERrayleigh25,'-*',SNR,BERrayleighadam,'-*');
legend('���ŵ�����BER','�������BER','�ض�����L=2','�ض�����L=4','�ض�����L=8','�ض�����L=10','�ض�����L=25','����Ӧ����');
xlabel('SNR (dB)'); 
ylabel('BER');
title('��ͬL0�ض�����������������');
function tmp=issame(input_data,decode)%�Ƚϱ���Σ����ڼ��������ʺ�У�麯�����ܵĹ��ߺ���
    tmp=0;
    for i=1:length(decode(:))
    if decode(i)==input_data(i)
        tmp=tmp+1;
    end
    end
end
function code_data=mycc2_1_2(input_data,D)%�������뺯��
    %code_data=zeros(1,2*length(input_data(:)));
    for i=1:length(input_data(:))
        y1=mod(input_data(i)+D(1)+D(2),2);
        y2=mod(input_data(i)+D(2),2);
        D(2)=D(1);
        D(1)=input_data(i);
        code_data(2*i-1)=y1;
        code_data(2*i)=y2;
    end
end
function m = viterbi_hard(x)   %viterbi�㷨���룬Ӳ�о�

% a 00
% b 10
% c 01
% d 11
%x = [ 1 1 0 1 0 1 0 0 0 1];
sizex = size(x);
s = sizex(2)/2;
% to record the value 
val_a = 0;
val_b = 0;
val_c = 0;
val_d = 0;

% ����ʽ��Ӧ��״��ͼ      aa0    ab1     bc0      bd1      ca0     cb1     dc0    dd1      
gra =      [ 0,0;    1,1;    1,0;    0,1;    1,1;    0,0;    0,1;    1,0  ];

% 4��·�� 
ma = zeros(1,s);
mb = zeros(1,s);
mc = zeros(1,s);
md = zeros(1,s);

%ǰ����ʱ�������·��
val_a = val_a + dis(gra(1,:), x(1:2));
ma(1)=0;
val_b = val_b + dis(gra(2,:), x(1:2));
mb(1)=1;

mc = mb;
md = mb;
val_c = val_b + dis(gra(3,:), x(3:4));
mc(2)=0;
val_d = val_b + dis(gra(4,:), x(3:4));
md(2)=1;

val_a = val_a + dis(gra(1,:), x(3:4));
ma(2)=0;
val_b = val_b + dis(gra(2,:), x(3:4));
mb(2)=1;

for i = 1:s
%     val_a_t =val_a;
%     val_b_t =val_b;
%     val_c_t =val_c;
%     val_d_t =val_d;
%     tempa = ma;
%     tempb = mb;
%     tempc = mc;
%     tempd = md;
    % for val_a
        if val_a + dis(gra(1,:), x(2*i-1:2*i)) >= val_c + dis(gra(5,:),x(2*i-1:2*i))%�ж���·��aa��ca��������a��a������c��a 
            tempa = mc; %ca�������̳�֮ǰ��C·����������״̬aΪ00��cΪ01
            val_a_t = val_c + dis(gra(5,:),x(2*i-1:2*i));%���´�ʱ��a״̬���ۼƾ��룬Ϊ֮ǰ��c�ľ���ʹ˿�c��a����֮�͡�
            tempa(i)=0;%��¼��ʱ�����루aa��ca��������0����������Ϊ������
        else
            val_a_t = val_a + dis(gra(1,:),x(2*i-1:2*i));
            tempa = ma;
            tempa(i)=0;
        end
      %for val_b
         if val_a + dis(gra(2,:), x(2*i-1:2*i)) >= val_c + dis(gra(6,:),x(2*i-1:2*i))
            tempb = mc; 
            val_b_t = val_c + dis(gra(6,:),x(2*i-1:2*i));
            tempb(i)=1;
        else
            val_b_t = val_a + dis(gra(2,:),x(2*i-1:2*i));
            tempb = ma;
            tempb(i)=1;         
         end
        
         %for val_c 
            if val_b + dis(gra(3,:), x(2*i-1:2*i)) >= val_d + dis(gra(7,:),x(2*i-1:2*i))
            tempc = md; 
            val_c_t = val_d + dis(gra(7,:),x(2*i-1:2*i));
            tempc(i)=0;
            else
            val_c_t = val_b + dis(gra(3,:),x(2*i-1:2*i));
            tempc = mb;
            tempc(i)=0;
            end
            
      %for val_c
            if val_b + dis(gra(4,:), x(2*i-1:2*i)) >= val_d + dis(gra(8,:),x(2*i-1:2*i))
            tempd = md; 
            val_d_t = val_d + dis(gra(8,:),x(2*i-1:2*i));
            tempd(i)=1;
        else
            val_d_t = val_b + dis(gra(4,:),x(2*i-1:2*i));
            tempd = mb;
            tempd(i)=1;
            end
    val_a =val_a_t;        
    val_b =val_b_t;
    val_c =val_c_t;
    val_d =val_d_t;
    ma = tempa;
    mb = tempb;
    mc = tempc;
    md = tempd;        
end

if val_a <= val_b
    m = ma;
    t = val_a;
else
    m = mb;
   t = val_b;
end

if val_c <= t
    m = mc;
    t =val_c;
end

if val_d <= t
    m = md;
    t = val_d;
end
end

function d = dis(x, y)   
  d = sum(xor(x,y));
end

function m = viterbi_short_hard(x,L)   %�ض��㷨��viterbi�㷨����,�ض����ΪL

% a 00
% b 10
% c 01
% d 11
%x = [ 1 1 0 1 0 1 0 0 0 1];
sizex = size(x);
s = sizex(2)/2;
% to record the value 
val_a = 0;
val_b = 0;
val_c = 0;
val_d = 0;
tmp=0;%���ڼ�¼·�����

% ����ʽ��Ӧ��״��ͼ      aa0    ab1     bc0      bd1      ca0     cb1     dc0    dd1      
gra =      [ 0,0;    1,1;    1,0;    0,1;    1,1;    0,0;    0,1;    1,0  ];

% 4��·�� 
ma = zeros(1,s);
mb = zeros(1,s);
mc = zeros(1,s);
md = zeros(1,s);

%ǰ����ʱ�������·��
val_a = val_a + dis(gra(1,:), x(1:2));
ma(1)=0;
val_b = val_b + dis(gra(2,:), x(1:2));
mb(1)=1;

mc = mb;
md = mb;
val_c = val_b + dis(gra(3,:), x(3:4));
mc(2)=0;
val_d = val_b + dis(gra(4,:), x(3:4));
md(2)=1;

val_a = val_a + dis(gra(1,:), x(3:4));
ma(2)=0;
val_b = val_b + dis(gra(2,:), x(3:4));
mb(2)=1;
tmp=2;
for i = 1:s
%     val_a_t =val_a;
%     val_b_t =val_b;
%     val_c_t =val_c;
%     val_d_t =val_d;
%     tempa = ma;
%     tempb = mb;
%     tempc = mc;
%     tempd = md;
    % for val_a
        if val_a + dis(gra(1,:), x(2*i-1:2*i)) >= val_c + dis(gra(5,:),x(2*i-1:2*i))%�ж���·��aa��ca��������a��a������c��a 
            tempa = mc; %ca�������̳�֮ǰ��C·����������״̬aΪ00��cΪ01
            val_a_t = val_c + dis(gra(5,:),x(2*i-1:2*i));%���´�ʱ��a״̬���ۼƾ��룬Ϊ֮ǰ��c�ľ���ʹ˿�c��a����֮�͡�
            tempa(i)=0;%��¼��ʱ�����루aa��ca��������0����������Ϊ������
        else
            val_a_t = val_a + dis(gra(1,:),x(2*i-1:2*i));
            tempa = ma;
            tempa(i)=0;
        end
      %for val_b
         if val_a + dis(gra(2,:), x(2*i-1:2*i)) >= val_c + dis(gra(6,:),x(2*i-1:2*i))
            tempb = mc; 
            val_b_t = val_c + dis(gra(6,:),x(2*i-1:2*i));
            tempb(i)=1;
        else
            val_b_t = val_a + dis(gra(2,:),x(2*i-1:2*i));
            tempb = ma;
            tempb(i)=1;         
         end
        
         %for val_c 
            if val_b + dis(gra(3,:), x(2*i-1:2*i)) >= val_d + dis(gra(7,:),x(2*i-1:2*i))
            tempc = md; 
            val_c_t = val_d + dis(gra(7,:),x(2*i-1:2*i));
            tempc(i)=0;
            else
            val_c_t = val_b + dis(gra(3,:),x(2*i-1:2*i));
            tempc = mb;
            tempc(i)=0;
            end
            
      %for val_c
            if val_b + dis(gra(4,:), x(2*i-1:2*i)) >= val_d + dis(gra(8,:),x(2*i-1:2*i))
            tempd = md; 
            val_d_t = val_d + dis(gra(8,:),x(2*i-1:2*i));
            tempd(i)=1;
        else
            val_d_t = val_b + dis(gra(4,:),x(2*i-1:2*i));
            tempd = mb;
            tempd(i)=1;
            end
    val_a =val_a_t;        
    val_b =val_b_t;
    val_c =val_c_t;
    val_d =val_d_t;
    ma = tempa;
    mb = tempb;
    mc = tempc;
    md = tempd;   
   %�ۼ���ȴﵽL������������������鲢�Ұ�����·��ɾ����ͬ��Ϊһ��·����
    if(mod(i,L)==0)
        bestval=min([val_a,val_b,val_c,val_d]);%ͨ������·�����ۼ�����о����·��
        if val_a==bestval%��һ·�������·��
            m = ma;%������飬�ɵ�һ·��ȷ��
            mb = ma;%��֦���������ڶ�������������·����ȡ���һ·����ͬ
            mc = ma;
            md = ma;
        end
         if val_b==bestval
            m = mb;
            ma = mb;
            mc = mb;
            md = mb;
         end
         if val_c==bestval
            m = mc;
            mb = mc;
            ma = mc;
            md = mc;
         end
         if val_d==bestval
            m = md;
            mb = md;
            mc = md;
            ma = md;
         end
        val_a=bestval;
        val_b=bestval;
        val_c=bestval;
        val_d=bestval;
        
    end
end
if val_a <= val_b
    m = ma;
    t = val_a;
else
    m = mb;
   t = val_b;
end

if val_c <= t
    m = mc;
    t =val_c;
end

if val_d <= t
    m = md;
    t = val_d;
end
end

function [m,L] = adam_viterbi(x)   %�ض��㷨��viterbi�㷨����,�ض����ΪL

% a 00
% b 10
% c 01
% d 11
%x = [ 1 1 0 1 0 1 0 0 0 1];
sizex = size(x);
s = sizex(2)/2;
% to record the value 
val_a = 0;
val_b = 0;
val_c = 0;
val_d = 0;
life=4;%������1,��������
thr=2;%������2,·������ֵ
survival=life;%��¼���·��
record=0;%����ͳ��ƽ�����
tmp=0;%���ڼ�¼·�����
times=0;
% ����ʽ��Ӧ��״��ͼ      aa0    ab1     bc0      bd1      ca0     cb1     dc0    dd1      
gra =      [ 0,0;    1,1;    1,0;    0,1;    1,1;    0,0;    0,1;    1,0  ];

% 4��·�� 
ma = zeros(1,s);
mb = zeros(1,s);
mc = zeros(1,s);
md = zeros(1,s);

%ǰ����ʱ�������·��
val_a = val_a + dis(gra(1,:), x(1:2));
ma(1)=0;
val_b = val_b + dis(gra(2,:), x(1:2));
mb(1)=1;

mc = mb;
md = mb;
val_c = val_b + dis(gra(3,:), x(3:4));
mc(2)=0;
val_d = val_b + dis(gra(4,:), x(3:4));
md(2)=1;

val_a = val_a + dis(gra(1,:), x(3:4));
ma(2)=0;
val_b = val_b + dis(gra(2,:), x(3:4));
mb(2)=1;
tmp=2;
for i = 1:s
%     val_a_t =val_a;
%     val_b_t =val_b;
%     val_c_t =val_c;
%     val_d_t =val_d;
%     tempa = ma;
%     tempb = mb;
%     tempc = mc;
%     tempd = md;
    % for val_a
        if val_a + dis(gra(1,:), x(2*i-1:2*i)) >= val_c + dis(gra(5,:),x(2*i-1:2*i))%�ж���·��aa��ca��������a��a������c��a 
            tempa = mc; %ca�������̳�֮ǰ��C·����������״̬aΪ00��cΪ01
            val_a_t = val_c + dis(gra(5,:),x(2*i-1:2*i));%���´�ʱ��a״̬���ۼƾ��룬Ϊ֮ǰ��c�ľ���ʹ˿�c��a����֮�͡�
            tempa(i)=0;%��¼��ʱ�����루aa��ca��������0����������Ϊ������
        else
            val_a_t = val_a + dis(gra(1,:),x(2*i-1:2*i));
            tempa = ma;
            tempa(i)=0;
        end
      %for val_b
         if val_a + dis(gra(2,:), x(2*i-1:2*i)) >= val_c + dis(gra(6,:),x(2*i-1:2*i))
            tempb = mc; 
            val_b_t = val_c + dis(gra(6,:),x(2*i-1:2*i));
            tempb(i)=1;
        else
            val_b_t = val_a + dis(gra(2,:),x(2*i-1:2*i));
            tempb = ma;
            tempb(i)=1;         
         end
        
         %for val_c 
            if val_b + dis(gra(3,:), x(2*i-1:2*i)) >= val_d + dis(gra(7,:),x(2*i-1:2*i))
            tempc = md; 
            val_c_t = val_d + dis(gra(7,:),x(2*i-1:2*i));
            tempc(i)=0;
            else
            val_c_t = val_b + dis(gra(3,:),x(2*i-1:2*i));
            tempc = mb;
            tempc(i)=0;
            end
            
      %for val_c
            if val_b + dis(gra(4,:), x(2*i-1:2*i)) >= val_d + dis(gra(8,:),x(2*i-1:2*i))
            tempd = md; 
            val_d_t = val_d + dis(gra(8,:),x(2*i-1:2*i));
            tempd(i)=1;
        else
            val_d_t = val_b + dis(gra(4,:),x(2*i-1:2*i));
            tempd = mb;
            tempd(i)=1;
            end
    val_a =val_a_t;        
    val_b =val_b_t;
    val_c =val_c_t;
    val_d =val_d_t;
    ma = tempa;
    mb = tempb;
    mc = tempc;
    md = tempd;   
   %�ۼ���ȴﵽL������������������鲢�Ұ�����·��ɾ����ͬ��Ϊһ��·����
    if(mod(i,1)==0)%�������ڣ����Ե��������ֶ��Ŵ�L�����ټ�������Ĭ��Ϊ1��ÿ�ζ�ִ�У�
        worstval=max([val_a,val_b,val_c,val_d]);%ͨ������·�����ۼ�����о����·��
        bestval=min([val_a,val_b,val_c,val_d]);
        if val_a==bestval%��¼���·��
            bestroute = ma;
        elseif val_b==bestval
            bestroute = mb;
        elseif val_c==bestval
            bestroute = mc;
        elseif val_d==bestval
             bestroute = md;
        end
        if(worstval-bestval)>=thr
            survival=survival-1;
        if val_a==worstval%��¼���·��
            ma=bestroute;
            val_a=bestval;
        elseif val_b==worstval
            mb=bestroute;
            val_b=bestval;
        elseif val_c==worstval
            mc=bestroute;
             val_c=bestval;
        elseif val_d==worstval
             md=bestroute;
             val_d=bestval;
        end
        else 
            if(worstval==bestval)
                survival=survival-1;end
            if(survival<=1)
               m=bestroute;
               survival=life;
               tmp=tmp+i-record;
               times=times+1;
               record=i;
            end   
        end
    end
end
if val_a <= val_b
    m = ma;
    t = val_a;
else
    m = mb;
   t = val_b;
end

if val_c <= t
    m = mc;
    t =val_c;
end

if val_d <= t
    m = md;
    t = val_d;
end
L=(tmp/times);
end
