function[x,fval,existFlag]=MyLPSolver2(f,A,b,Aeq,beq)
%xΪ���Ž⣬fvalΪ���ź���ֵ
%existFlag==1:Ψһ���Ž⣻existFlag==2:��������Ž⣻existFlag==3:�޽磻existFlag==4:�޿��н⣻existFlag==5:��ⳬʱ��
iter=0;%��ʼ����������
% b=[5;-2];
% f=[-3 1 -2];
% A=[1 3 1 ;-2 1 -1];
% Aeq=[4 3 -2];
% beq=[5];
b=[5;2];
f=[-3 1 -2];
A=[1 -3 1 ;-2 1 -1];
Aeq=[-4 -3 -2];
beq=[55];


% ����A Aeq���������Բ�һ�� �������0����
if length(A)>length(Aeq)
    l=length(A);
    A=[padarray(A,[l-length(A) 0],'post'); padarray(Aeq',[l-length(Aeq) 0],'post')'];
    
elseif length(A)<length(Aeq)
    l=length(Aeq) ;
    A=[padarray(A',[l-length(A) 0],'post')'; padarray(Aeq,[l-length(Aeq) 0],'post')];
else
    A=[A; Aeq];
end
b=[b;beq];

%���ɱ�׼��
[mA nA]=size(A);    %ȡA������mA,����nA
[mAeq nAeq]=size(Aeq);    %ȡAeq������mAeq,����nAeq
Arti1=mA+nA-mAeq+1; %ȷ���˹�������Χ Arti1Ϊ��ʼ�±�

slackV=eye(mA);        %ȡ��������ȵĵ�λ�󣬶�A:�������ɳڱ���  ��Aeq:�������˹�����
A1=[A slackV] ;        %��׼�������Ա���ϵ������


m=[(nA+1):(mA+nA)]; %��ʼ���������±�
mm=[1:nA];         %��ʼ�ǻ��������±ꣻ

tmp1=[zeros(1,(mA+nA-size(f,2)))];%�����ɳں��˹�������Ӧ�ļ�����
f1=[f tmp1];         %���������ļ�����,��Ϊf1

AA1=[A1 b]%�����α�
existFlag=0;
Arti1IsOut=false;
k=0;
k_l=Arti1;
r=1;
while existFlag==0%��ʼ����
    iter=iter+1;%����������1
    b=AA1(:,mA+nA+1) ;
    for i=1:mA
        if m(i)==k_l;
            m(i)=k;;
        end
    end
    for i=1:nA
        if mm(i)==k;
            mm(i)=k_l;
        end
    end
    slackV=AA1(:,m);%������
    %������ʼ�����α�
    fb=[];
    for i=1:mA
        fb=[fb f1(m(i))];%��������Ӧ������������
    end
    
    xb=(inv(slackV))*b ;%����x
    z0=fb*xb;%����Ŀ�꺯��
    B=zeros(1,mA+nA);%����������
    %�����������
    for j=1:nA
        cy=0;
        for i=1:mA
            cy=cy+f1(m(i))*AA1(i,mm(j));
        end
        z(mm(j))=cy;
        B(mm(j))=f1(mm(j))-z(mm(j));
    end
    tmp2=[B -z0];
    AA=[tmp2;AA1]    %���µ����α�
    
    if (iter)>1000 %��ⳬʱʱ��ֹ���
        existFlag=5
        disp('��ⳬʱ')
        break;
    end
    
    %���Ž����
    
    if (min(B))>=0 && Arti1IsOut==true%�˹�����������
        nbv=B(:,mm);%�ǻ�������Ӧ�ļ�����
        if max(m)>Arti1%������Ž⺬���˹���������α���Ž�->�޿��н�
            disp('�޿��н�')
            existFlag=4
            x(m)=xb
            fval=[]
            break;
            
            
        elseif min(nbv)==0%���ڷǻ������ļ�����Ϊ0
            x(m)=xb
            disp('����������Ž�')
            existFlag=2
            fval=z0
            break;
            
        else
            x(m)=xb
            disp('��Ψһ���Ž�')
            existFlag=1
            disp('���Ž�x=')
            x=x(:,1:max(nA,nAeq))
            fval=z0 %������Ψһ���Ž⣬��ʾ��������ֵfval�Լ�ȡֵx
            disp('��������=')
            disp(iter)
            break;
        end
    end
    
    if (min(B))<0         %���ڼ�����<0��ת��
        k=min(find(B==min(B(find(B<0)))));%ѡֵ��С����Ϊpivot
        pivot=AA1(:,k)
        if max(pivot)<=0   %���ڼ�����<0,�Ҷ�Ӧ�ġ�pivot���������Ԫ�ض�<=0�����޽��
            disp('���޽��')
            existFlag=3
            x=[]
            fval=[]
            break;
        end
   
        if   Arti1IsOut==false
    r=find(tmp3==min(tmp3(find(tmp3>0))));
    pri=AA1(r,k);            %�ҵ���Ԫ��pri������           
    leave=AA1(:,k_l)          %�����±�Ϊk_l %����leave
    
        if(k_l+1<=mA+nA)
            k_l=k_l+1
        else
            Arti1IsOut=true;
        end
        
        else
    %ȷ����������leave
    for i=1:mA
        tmp3(i)=(AA1(i,mA+nA+1))/(AA1(i,k));
    end
    r=find(tmp3==min(tmp3(find(tmp3>0))));
    pri=AA1(r,k);            %�ҵ���Ԫ��pri
    k_l=m(r);                %�����±�Ϊk_l
    leave=AA1(:,k_l)          %����leave
        end
     
    %����Ԫ�ĸ�˹��Ԫ�������Է����飬���µ����α�AA1
    for i=1:mA
        if i~=r
            AA1(i,:)=AA1(i,:)-AA1(r,:)*AA1(i,k)/AA1(r,k)
        else
            AA1(r,:)=AA1(r,:)/pri
        end
    end
    AA1
        end
   
end



