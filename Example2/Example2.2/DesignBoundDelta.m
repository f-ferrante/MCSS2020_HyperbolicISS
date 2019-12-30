
%This code solves the optimization problem considered in the paper for Example 2, Second Scenario.
%The codes requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used.    

clear all;

%%%%%%%Plant Data Definition%%%%%%%%%%

Lambda =diag([1, 1, 2]);
H=[0.2 0.4 0.2; 0.8, 0.2,0.1; 0.4 0 0.2]; 
B=[0 0 1; 1 0 0]'; 
np=max(size(H));
nu=min(size(B));

%%%%%%%Optimization problem Definition%%%%%%%%%
mu=0.3209*2;;
delta=1;
c=sdpvar(1,1,'full');
W=diag(sdpvar(np,1));
F=sdpvar(np,np,'full');
Y=sdpvar(nu,np,'full');
s=exp(-mu)*mu*min(eig(Lambda));

 M=[F+F', -(H*W+B*Y), -B*Y, F';
   -(H*W+B*Y)', -exp(-mu)*W*Lambda, zeros(np,np),zeros(np,np);
    (-B*Y)', zeros(np,np),eye(np)*s*(-2*delta*W+delta^2*eye(np)), zeros(np,np);
    F, zeros(np,np),zeros(np,np),-inv(Lambda)*W
    ];

problem=[M<=-1e-8*eye(max(size(M))), W<=c*eye(np), c>=0];
    
options=sdpsettings('solver','sdpt3','verbose',0);

solution=solvesdp(problem,c,options);
if solution.problem==0
W=double(W);
P=inv(W);
Y=double(Y);
K=Y*inv(W);
gamma=sqrt(1/min(eig(P)));
end
%%%%%%%%%%%%%%Definition of global variables for simulations%%%%%%%%%%%%
global H B K Lambda

