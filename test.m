%part1 ���Ž����
%existFlag==1:Ψһ���Ž⣻existFlag==2:��������Ž⣻existFlag==3:�޽磻existFlag==4:�޿��н⣻existFlag==5:��ⳬʱ��
%�ֱ�������ֲ���������ʽ

%��ʽ1 MyLPSolver(f,A,b)
%Ψһ���Ž�
% b=[2;6];
% f=[-2 -1];
% A=[1 -1;1 1];
%��������Ž�
% f = [-2 -4];
% A = [-1 2 ;
%     1 2;
%     1 -1];
% b = [4;10;2];
%�޽��
% f = [-1 1 -1];
% A = [2 -1 -2 1 0 0 ;
%     -2 3 1 0 -1 0;
%     1 -1 -1 0 0 -1 ];
% b = [4;5;1];

%  [x,fval,existFlag]=MyLPSolver1(f,A,b)

%��ʽ2 MyLPSolver2(f,A,b,Aeq,beq)
%Ψһ���Ž�
% b=[5;-2];
% f=[-3 1 -2];
% A=[1 3 1 ;-2 1 -1];
% Aeq=[4 3 -2];
% beq=[5]; 
%��������Ž�
% f = [-2 -4];
% A = [-1 2 ;
%     1 2;
%     1 -1];
% b = [4;10;2];
% Aeq=[];
% beq=[];
%�޽��
% b=[5;2];
% f=[-3 1 -2];
% A=[1 -3 1 ;-2 1 -1];
% Aeq=[-4 -3 -2];
% beq=[55];

%  [x,fval,existFlag]=MyLPSolver2(f,A,b,Aeq,beq)
%  lb=[0;0;0]; ub=[100;100;100];
% [x,fval,existFlag]=linprog(f,A,b,Aeq,beq,lb,ub)

%��ʽ3 MyLPSolver3(f,A,b,Aeq,beq,lb,ub)
% b=[5;-2];
% f=[-3 1 -2];
% A=[1 3 1 ;-2 1 -1];
% Aeq=[4 3 -2];
% beq=[5];
% lb=[-1;-2;-3];
% ub=[11;12;13];
%  [x,fval,existFlag]=MyLPSolver3(f,A,b,Aeq,beq,lb,ub)


%part2 A ����Ĺ�ģ n ��С����ʱ��������������Թ滮��ģ������ʱ��ı仯 ����Ϊ���������㷨���ƣ�����ֻ�õ�һ�����ԣ�

    tic
   for n=2:100 %ע������ʱ��ȽϾ� 
 tic 
 b=randi( [-10,100],n,1) ;
 f=randi( [-10,100],1,n) ;
A=randi( [-10,100],n ,n) ;
[x,fval,existFlag]=MyLPSolver1(f,A,b)
dt = toc;
DT(n)=dt;
   end
   
  DT
 fid = fopen('D:\college\�����\�Ż�\��ĩ1\my\results.txt','wt');%������ʱ��洢��results�ļ���
fprintf(fid,'%g\n',DT);       
fclose(fid);
% % % save D:\college\�����\�Ż�\��ĩ1\my\results.txt DT -ascii ; %�洢������
% % % type('results.txt')

% min(DT,[],'all')
% max(DT,[],'all')
% mean(DT,'all')


%part3 ʵ�����Ч��
%����ʵ����������Ч��
%   problem = mpsread('adlittle.mps');
% %  problem = mpsread('afiro.mps')
%   tic;
%  [x,fval,existFlag] = MyLPSolver3(full(problem.f)',full(problem.Aineq),full(problem.bineq),full(problem.Aeq),full(problem.beq),full(problem.lb),full(problem.ub));
%   % [x,fval,existFlag] = MyLPSolver2(full(problem.f)',full(problem.Aineq),full(problem.bineq),full(problem.Aeq),full(problem.beq));
%  toc
