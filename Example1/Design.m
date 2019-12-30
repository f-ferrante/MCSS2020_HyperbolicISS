
%This code solves the optimization problem considered in the paper. 
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used.    

clear all;

%%%%%%%Plant Date Definition%%%%%%%%%%
Lambda =[1 0; 0 sqrt(2)];

H=[0 1.1;
   1 0]; 
B=eye(2); 
np=max(size(H));
nu=min(size(B));
%%%%%%%Optimization problem Definition%%%%%%%%%%
mu=0.1;
c=sdpvar(1,1,'full');
P=diag(sdpvar(np,1));
Q=sdpvar(np,np,'full');
Y=sdpvar(np,np,'full');
s=sqrt(exp(-mu)*min(eig(Lambda)));

 M=[Q+Q'+Lambda*P, -(Q'*H+Y), -Y;
    -(Q'*H+Y)', -exp(-mu)*P*Lambda, zeros(np,np);
    -Y', zeros(np,np),-eye(np)*s^2];

problem=[M<=-1e-8*eye(max(size(M))), P>=c*eye(np), c>=0, P>=1e-6*eye(np)];
    
options=sdpsettings('solver','sdpt3','verbose',0);

solution=solvesdp(problem, -c,options);

if(solution.problem==0)
P=double(P);
Y=double(Y); 
Q=double(Q);
K=inv(B)*inv(Q')*Y;
gamma=sqrt(1/min(eig(P)));
end
%%%%%%%%%%%%%%Definition of global variables for simulations%%%%%%%%%%%%
global H B K Lambda

