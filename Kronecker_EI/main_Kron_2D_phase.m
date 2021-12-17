%%%%%%%Script for computing the 2D Allen-Cahn equation with the %%%%%%%
%%%%%%%Exponential Euler method%%%%%%%

%Import the data of the problem
[eps,beta,u0,x1,x2,y1,y2,T,P]=data;
delta=0.01; %coefficient of the non-linear term 
pert=delta^2; %perturbation parameter

%Number of time steps, time step size and time 
steps=200;
tau=T/steps;
t=0:tau:T;

%Number of integration points (Lobatto quadrature)
nquad=2;

%The number of elements (2^r) in each space direction 
r=9;


%Meshes and Parameters
[xsol,ysol,nelx,nely] = mesh2D(r,x1,x2,y1,y2,P); 
[nx,ny,nel,nnode,coord,nodes]=parameters(xsol,ysol);

%Full sparse matrix
dir=BC(nx,ny,P); 
A=FEM_matrices_sparse(eps,beta,dir,nel,coord,nodes,nquad);
A=-tau*A;

%Store the size of the full matrix
size_A=size(A,1);

%1D matrices
[dirx,diry]=BC_1D(nx,ny,P); 
Ax=FEM_matrices_1D(eps,beta(1),dirx,xsol,nelx,nquad,pert);
Ay=FEM_matrices_1D(eps,beta(2),diry,ysol,nely,nquad,pert);
Ax=-tau*Ax;
Ay=-tau*Ay;

%Sanity check for A having Kronecker sum structure
%Akron=kron(eye(size(Ay)),Ax)+kron(Ay,eye(size(Ax)));
%norm(A-Akron)

%Initial condition vector
[xmat,ymat]=meshgrid(xsol,ysol);
U0=u0(xmat,ymat);
U0=reshape(U0',[],1);
U0(dir)=[];
dimx=size(U0,1);

%Reshape initial condition vector
[p,q]=size(Ax);
[m,n]=size(Ay);

%Solve full 2D system
%Initialize solution vector
uhat_2D=zeros(dimx,steps+1);
uhat_2D(:,1)=U0;
tic
%Loop through time steps
for i=1:steps
    %Action of Varphi_0 over the previous solution
    V0=expmv(1,A,uhat_2D(:,i));

    %Action of Varphi_1 over the non-linear term
    B=sparse(dimx+1,1,1,dimx+1,1);
    F=(uhat_2D(:,i)-uhat_2D(:,i).^3)/delta^2;
    V1=expmv(1,block_matrix(A,F,dimx,1),B);
    V1=V1(1:end-1);

    %Sum both contributions 
    uhat_2D(:,i+1)=V0+tau*V1;
end
times_full=toc;

%Solve the 1D system systems 
%Initialize solution vector
uhat_kron=zeros(dimx,steps+1);
uhat_kron(:,1)=U0;

tic
%Compute the Varphi matrices once
Phi0x=expm(Ax);
Phi0y=expm(Ay);
Phi1x=Ax*phipade(Ax,1);
Phi1y=Ay*phipade(Ay,1);
%Loop trough time steps 
for i=1:steps
    %Action of Varphi_0 over the previous solution
    U=reshape(uhat_kron(:,i),q,n);
    W0=Phi0x*U*Phi0y';
    W0kron=reshape(W0,m*p,1);

    %Action of Varphi_1 over the non-linear term
    F=reshape((uhat_kron(:,i)-uhat_kron(:,i).^3)/delta^2-pert*uhat_kron(:,i),q,n);
    W1=Phi1x*F*Phi1y'+F*Phi1y'+Phi1x*F;

    %Solve Sylvester equation
    W1=sylvester(Ax,Ay',W1);
    W1kron=reshape(W1,m*p,1);

    %Solve the full 2D system
    %W1kron=A\W1;

    %Sum both contributions 
    uhat_kron(:,i+1)=W0kron+tau*W1kron;
end
times_kron=toc;


%Compute the relative error of the solution between both approaches
error=zeros(steps+1,1);
for i=1:steps+1
    error(i)=norm(uhat_2D(:,i)-uhat_kron(:,i))/norm(uhat_2D(:,i));
end

%Export the times and errors to a .txt file 
Tab_error=table(t',error);
writetable(Tab_error,'RelativeErrorCahn_Euler_r9_2000.txt','Delimiter',' ');

saving=times_full./times_kron;
Tab_times=table(size_A,times_full,times_kron,saving);
writetable(Tab_times,'TimesCahn_Euler_r9_2000.txt','Delimiter',' ');

%Plot relative error 
figure
plot(t,error,'-b*')
xlabel({'$t$'},'interpreter','latex')
ylabel({'Relative Error'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)
    


    






