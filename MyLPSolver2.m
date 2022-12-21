function[x,fval,existFlag]=MyLPSolver2(f,A,b,Aeq,beq)
%x为最优解，fval为最优函数值
%existFlag==1:唯一最优解；existFlag==2:无穷多最优解；existFlag==3:无界；existFlag==4:无可行解；existFlag==5:求解超时；
iter=0;%初始化迭代次数
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


% 假设A Aeq的列数可以不一样 因此作补0调整
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

%化成标准型
[mA nA]=size(A);    %取A的行数mA,列数nA
[mAeq nAeq]=size(Aeq);    %取Aeq的行数mAeq,列数nAeq
Arti1=mA+nA-mAeq+1; %确定人工变量范围 Arti1为起始下标

slackV=eye(mA);        %取与行数相等的单位阵，对A:即引入松弛变量  对Aeq:即引入人工变量
A1=[A slackV] ;        %标准单纯形自变量系数矩阵


m=[(nA+1):(mA+nA)]; %初始基变量的下标
mm=[1:nA];         %初始非基变量的下标；

tmp1=[zeros(1,(mA+nA-size(f,2)))];%构造松弛和人工变量对应的检验数
f1=[f tmp1];         %生成完整的检验数,记为f1

AA1=[A1 b]%单纯形表
existFlag=0;
Arti1IsOut=false;
k=0;
k_l=Arti1;
r=1;
while existFlag==0%开始迭代
    iter=iter+1;%迭代次数加1
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
    slackV=AA1(:,m);%基变量
    %建立初始单纯形表；
    fb=[];
    for i=1:mA
        fb=[fb f1(m(i))];%基变量对应的完整检验数
    end
    
    xb=(inv(slackV))*b ;%更新x
    z0=fb*xb;%计算目标函数
    B=zeros(1,mA+nA);%完整检验数
    %计算检验数；
    for j=1:nA
        cy=0;
        for i=1:mA
            cy=cy+f1(m(i))*AA1(i,mm(j));
        end
        z(mm(j))=cy;
        B(mm(j))=f1(mm(j))-z(mm(j));
    end
    tmp2=[B -z0];
    AA=[tmp2;AA1]    %更新单纯形表
    
    if (iter)>1000 %求解超时时终止求解
        existFlag=5
        disp('求解超时')
        break;
    end
    
    %最优解分析
    
    if (min(B))>=0 && Arti1IsOut==true%人工变量出过基
        nbv=B(:,mm);%非基变量对应的检验数
        if max(m)>Arti1%如果最优解含有人工变量，是伪最优解->无可行解
            disp('无可行解')
            existFlag=4
            x(m)=xb
            fval=[]
            break;
            
            
        elseif min(nbv)==0%存在非基变量的检验数为0
            x(m)=xb
            disp('有无穷多最优解')
            existFlag=2
            fval=z0
            break;
            
        else
            x(m)=xb
            disp('有唯一最优解')
            existFlag=1
            disp('最优解x=')
            x=x(:,1:max(nA,nAeq))
            fval=z0 %否则是唯一最优解，显示变量最优值fval以及取值x
            disp('迭代次数=')
            disp(iter)
            break;
        end
    end
    
    if (min(B))<0         %存在检验数<0，转轴
        k=min(find(B==min(B(find(B<0)))));%选值最小的作为pivot
        pivot=AA1(:,k)
        if max(pivot)<=0   %对于检验数<0,且对应的‘pivot’中任意的元素都<=0，是无界解
            disp('有无界解')
            existFlag=3
            x=[]
            fval=[]
            break;
        end
   
        if   Arti1IsOut==false
    r=find(tmp3==min(tmp3(find(tmp3>0))));
    pri=AA1(r,k);            %找到主元素pri所在行           
    leave=AA1(:,k_l)          %出基下标为k_l %出基leave
    
        if(k_l+1<=mA+nA)
            k_l=k_l+1
        else
            Arti1IsOut=true;
        end
        
        else
    %确定出基变量leave
    for i=1:mA
        tmp3(i)=(AA1(i,mA+nA+1))/(AA1(i,k));
    end
    r=find(tmp3==min(tmp3(find(tmp3>0))));
    pri=AA1(r,k);            %找到主元素pri
    k_l=m(r);                %出基下标为k_l
    leave=AA1(:,k_l)          %出基leave
        end
     
    %列主元的高斯消元法解线性方程组，更新单纯形表AA1
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



