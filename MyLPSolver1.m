function[x,fval,existFlag]=MyLPSolver1(f,A,b)
%xΪ���Ž⣬fvalΪ���ź���ֵ
%existFlag==1:Ψһ���Ž⣻existFlag==2:��������Ž⣻existFlag==3:�޽磻existFlag==4:�޿��н⣻existFlag==5:��ⳬʱ��
iter=0;%��ʼ����������
% b=[2;6];
% f=[-2 -1];
% A=[1 -1;1 1];

% f = [-2 -4];
% A = [-1 2 ;
%     1 2;
%     1 -1];
% b = [4;10;2];

%���ɱ�׼��
[mA nA]=size(A);    %ȡA������mA,����nA
slackV=eye(mA);        %ȡ��������ȵĵ�λ�󣬼�ͳһ�����ɳڱ���  
A1=[A slackV]  ;       %��׼�������Ա���ϵ������


m=[(nA+1):(mA+nA)]; %��ʼ���������±�
mm=[1:nA];         %��ʼ�ǻ��������±ꣻ

tmp1=zeros(1,(mA+nA-size(f,2)));%�����ɳڱ�����Ӧ�ļ�����
f1=[f tmp1];         %���������ļ�����,��Ϊf1

AA1=[A1 b]%�����α�
existFlag=0;
k=0;
k_l=0;
r=1;
while existFlag==0%��ʼ����
    iter=iter+1;%����������1
    b=AA1(:,mA+nA+1);  %����b
    for i=1:mA
        if m(i)==k_l;
            m(i)=k;
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
    B=zeros(1,mA+nA);
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
    nbv=B(:,mm);%�ǻ�������Ӧ�ļ�����
    if (min(B))>=0
        if min(nbv)==0
            disp('����������Ž�')
            existFlag=2
            x(m)=xb
          fval=z0         
            break;
            
        else
                x(m)=xb 
              %  if min(x)>=0              
                disp('��Ψһ���Ž�')
                existFlag=1
                disp('���Ž�x=')
                 x=x(:,1:nA)
                fval=z0 %�ж��Ƿ���Ψһ���Ž⣬���������ʾ����ȡֵx�Լ�����ֵfval
                disp('��������=')
                 disp(iter)
                break;
           % end 
        end
    end
            
            %����ת��
            k=min(find(B==min(B(find(B<0)))));%ѡֵ��С����Ϊpivot
            pivot=AA1(:,k)
            if max(pivot)<=0   %���ڼ�����<0,�Ҷ�Ӧ�ġ�pivot���������Ԫ�ض�<=0����Ŀ�꺯�����޽�ģ������Ž⣻
                 disp('���޽��')
                 existFlag=3
                 x=[]
                fval=[]
                break;
            end
            %ȷ����������leave
            for i=1:mA
                tmp3(i)=(AA1(i,mA+nA+1))/(AA1(i,k));
            end
            r=find(tmp3==min(tmp3(find(tmp3>0))));
            r
            k
            pri=AA1(r,k)            %�ҵ���Ԫ��pri
            k_l=m(r);                %�����±�Ϊk_l
            leave=AA1(:,k_l)          %����leave
            
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
        
        
        
