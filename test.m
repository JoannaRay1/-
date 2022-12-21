%part1 最优解测试
%existFlag==1:唯一最优解；existFlag==2:无穷多最优解；existFlag==3:无界；existFlag==4:无可行解；existFlag==5:求解超时；
%分别测试三种参数输入形式

%形式1 MyLPSolver(f,A,b)
%唯一最优解
% b=[2;6];
% f=[-2 -1];
% A=[1 -1;1 1];
%无穷多最优解
% f = [-2 -4];
% A = [-1 2 ;
%     1 2;
%     1 -1];
% b = [4;10;2];
%无界解
% f = [-1 1 -1];
% A = [2 -1 -2 1 0 0 ;
%     -2 3 1 0 -1 0;
%     1 -1 -1 0 0 -1 ];
% b = [4;5;1];

%  [x,fval,existFlag]=MyLPSolver1(f,A,b)

%形式2 MyLPSolver2(f,A,b,Aeq,beq)
%唯一最优解
% b=[5;-2];
% f=[-3 1 -2];
% A=[1 3 1 ;-2 1 -1];
% Aeq=[4 3 -2];
% beq=[5]; 
%无穷多最优解
% f = [-2 -4];
% A = [-1 2 ;
%     1 2;
%     1 -1];
% b = [4;10;2];
% Aeq=[];
% beq=[];
%无界解
% b=[5;2];
% f=[-3 1 -2];
% A=[1 -3 1 ;-2 1 -1];
% Aeq=[-4 -3 -2];
% beq=[55];

%  [x,fval,existFlag]=MyLPSolver2(f,A,b,Aeq,beq)
%  lb=[0;0;0]; ub=[100;100;100];
% [x,fval,existFlag]=linprog(f,A,b,Aeq,beq,lb,ub)

%形式3 MyLPSolver3(f,A,b,Aeq,beq,lb,ub)
% b=[5;-2];
% f=[-3 1 -2];
% A=[1 3 1 ;-2 1 -1];
% Aeq=[4 3 -2];
% beq=[5];
% lb=[-1;-2;-3];
% ub=[11;12;13];
%  [x,fval,existFlag]=MyLPSolver3(f,A,b,Aeq,beq,lb,ub)


%part2 A 矩阵的规模 n 从小到大时，求解器随着线性规划规模变大求解时间的变化 （因为三个函数算法类似，这里只用第一个测试）

    tic
   for n=2:100 %注意运行时间比较久 
 tic 
 b=randi( [-10,100],n,1) ;
 f=randi( [-10,100],1,n) ;
A=randi( [-10,100],n ,n) ;
[x,fval,existFlag]=MyLPSolver1(f,A,b)
dt = toc;
DT(n)=dt;
   end
   
  DT
 fid = fopen('D:\college\大二下\优化\期末1\my\results.txt','wt');%将运行时间存储在results文件里
fprintf(fid,'%g\n',DT);       
fclose(fid);
% % % save D:\college\大二下\优化\期末1\my\results.txt DT -ascii ; %存储方法二
% % % type('results.txt')

% min(DT,[],'all')
% max(DT,[],'all')
% mean(DT,'all')


%part3 实际求解效果
%测试实际问题的求解效果
%   problem = mpsread('adlittle.mps');
% %  problem = mpsread('afiro.mps')
%   tic;
%  [x,fval,existFlag] = MyLPSolver3(full(problem.f)',full(problem.Aineq),full(problem.bineq),full(problem.Aeq),full(problem.beq),full(problem.lb),full(problem.ub));
%   % [x,fval,existFlag] = MyLPSolver2(full(problem.f)',full(problem.Aineq),full(problem.bineq),full(problem.Aeq),full(problem.beq));
%  toc
