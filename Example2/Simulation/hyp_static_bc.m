
function [x1,t1, V1to, V2to, V3to]=hyp_static_bc
% Code originally written by Ying Tang to illustrate a CDC2014 paper
% and then adapted by Francesco Ferrante to illustrate a CDC2018 paper. The
% codes relies on the software for the solution to hyperbolic PDEs written
% by Shampine and freely available at http://faculty.smu.edu/shampine/current.html

global  K H lamb1 lamb2 lamb3 Thorizon Lambda
lamb1 = Lambda(1,1);
lamb2 = Lambda(2,2);
lamb3 = Lambda(3,3);
Thorizon =50;
method = 'LxF';
pdefun = @pdes;

% Mesh grid
Nx = 10; %number of points in the space domain   
x = linspace(0,1,Nx);
Tinit = 0;
t = Tinit;
% Initial conditions %%%% SHOULD SATISFY THE COMPATIBILITY CONDITION
V(1,:)=(cos(4*pi*x)-1);
V(2,:)=(cos(2*pi*x)-1);
V(3,:)=x.*(1-x);

sol = setup(1,pdefun,t,x,V,method,[],@bcfun);

% check CFL condition
dx = x(2)-x(1);
dt = 0.95*dx/4;

Nt=floor(Thorizon/dt);
howfar = Thorizon/Nt;

%to stack the solution
V1to=zeros(Nt,Nx);
V2to =zeros(Nt,Nx);
V3to =zeros(Nt,Nx);
V1to(1,:) = V(1,:);
V2to(1,:) = V(2,:);
V3to(1,:) = V(3,:);
for m = 1:Nt-1
    disp([num2str(m) '/' num2str(Nt)]) 
    %keyboard
    sol = hpde(sol,howfar,dt);  
    t = sol.t;
    V = sol.u;
    V1to(m+1,:) = V(1,:);
    V2to(m+1,:) = V(2,:);
    V3to(m+1,:) = V(3,:);

end

%%%% Generate the two vectors z and t
t = linspace(Tinit,Thorizon,Nt);
[x1,t1] = meshgrid(x,t);
end

    
%=========================================================================
% Subfunctions

    function F = pdes(t,x,V,V_x) % linear PDE
        global lamb1 lamb2 lamb3
        F = zeros(size(V));
        F(1,:) = -lamb1.*V_x(1,:);
        F(2,:) = -lamb2.*V_x(2,:);
        F(3,:) = -lamb3.*V_x(3,:);
    % end function pdes
    end
   
    function [XL,XR] = bcfun(t,XLex,XRex) 
    %boundary conditions
    global K B H
    XLtmp = XLex;
    Xbc=(H+B*K)*[XRex(1);XRex(2);XRex(3)]+0.2*B*K*[1;1;1]*sin(t);%boundary condition
    XLtmp(1) = [1 0 0]*Xbc;
    XLtmp(2) = [0 1 0]*Xbc;
    XLtmp(3) = [0 0 1]*Xbc;
    XRtmp = XRex;
    XL=XLtmp;
    XR=XRtmp;
   end