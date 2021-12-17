%%%%%%%Script for computing the computational times and errors of%%%%%%%
%%%%%%%the 2D heat equation and Eriksson-Johnson equation%%%%%%%

%Import the data of the problem
[eps,beta,u0,x1,x2,y1,y2,T,P]=data;

%Time step size and the order of the varphi-function
tau=1/16;
order=5;

%Number of integration points (Lobatto quadrature)
nquad=2;

%The initial and final number of elements (2^r) in each space direction 
r_in=3;
r_max=7;

%Initialize the vectors for the computational times, size of matrices and
%relative errors
times_full=zeros(r_max-r_in+1,1);
times_kron=zeros(r_max-r_in+1,1);
size_A=zeros(r_max-r_in+1,1);
error=zeros(r_max-r_in+1,1);

%Initialize the perturbartion parameter (only for the Laplacian with 
%homogeneous Neumann BC)
pert=0;

%Loop through the global space refinaments 
for r=r_in:r_max
    r

    %Meshes and Parameters
    [xsol,ysol,nelx,nely] = mesh2D(r,x1,x2,y1,y2,P);
    [nx,ny,nel,nnode,coord,nodes]=parameters(xsol,ysol);

    %Full sparse matrix
    dir=BC(nx,ny,P);
    A=FEM_matrices_sparse(eps,beta,dir,nel,coord,nodes,nquad);
    A=-tau*A;
    condest(A)

    %Store the size of the full matrix
    size_A(r-r_in+1)=size(A,1);

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
    U=reshape(U0,q,n);

    if order==0
        %Varphi order: p=0

        %Time for the 2D matrix
        tic
        V0=expmv(1,A,U0);
        times_full(r-r_in+1)=toc;

        %Time for the Kronecker product
        tic
        W=expm(Ax)*U*expm(Ay');
        W0kron=reshape(W,m*p,1);
        times_kron(r-r_in+1)=toc;

        %Relative error
        error(r-r_in+1)=norm(V0-W0kron)/norm(V0);

    else
        %Varphi order: p>=1

        %Time for the 2D matrix
        tic
        B=sparse(dimx+order,1,1,dimx+order,1);
        V0=expmv(1,block_matrix(A,U0,dimx,order),B);
        V0=V0(1:end-order);
        times_full(r-r_in+1)=toc;

        %Compute the actions for p>=1
        if order==1

        tic
        Phi1x=Ax*phipade(Ax,1); %It is not necessary to store this matrices!
        Phi1y=Ay*phipade(Ay,1);
        W=Phi1x*U*Phi1y'+U*Phi1y'+Phi1x*U;

        elseif order==2

        tic
        Phi1x=phipade(Ax,1);
        Phi1y=phipade(Ay,1);
        Phi2x=phipade(Ax,2);
        Phi2y=phipade(Ay,2);
        W=Ax*Phi1x*U*(Ay*Phi1y)'+action(Ax,Phi2x*U,2)+action(Ay,Phi2y*U',2)';

        elseif order==3

        tic
        Phi2x=phipade(Ax,2);
        Phi2y=phipade(Ay,2);
        Phi3x=phipade(Ax,3);
        Phi3y=phipade(Ay,3);
        W=action(Ay,Phi2y*(action(Ax,Phi2x*U,2))',2)'+action(Ay,Phi2y*U'*Ax',2)'...
          +action(Ax,Phi2x*U*Ay',2)+action(Ax,Phi3x*U,3)+action(Ay,Phi3y*U',3)';

        elseif order==4

        tic
        Phi2x=phipade(Ax,2);
        Phi2y=phipade(Ay,2);
        Phi3x=phipade(Ax,3);
        Phi3y=phipade(Ay,3);
        Phi4x=phipade(Ax,4);
        Phi4y=phipade(Ay,4);
        W=action(Ay,Phi2y*(action(Ax,Phi2x*U,2))',2)'+action(Ay,Phi3y*U'*Ax',3)'...
          +action(Ax,Phi3x*U*Ay',3)+action(Ax,Phi4x*U,4)+action(Ay,Phi4y*U',4)';

        elseif order==5
        tic
        Phi2x=phipade(Ax,2);
        Phi2y=phipade(Ay,2);
        Phi3x=phipade(Ax,3);
        Phi3y=phipade(Ay,3);
        Phi4x=phipade(Ax,4);
        Phi4y=phipade(Ay,4);
        Phi5x=phipade(Ax,5);
        Phi5y=phipade(Ay,5);
        W=action(Ay,Phi3y*(action(Ax,Phi3x*U,3))',3)'+(1/2)*action(Ay,Phi3y*action(Ax,U,2)',3)'...
          +(1/2)*action(Ay,action(Ax,Phi3x*U,3)',2)'+action(Ay,Phi4y*U'*Ax',4)'...
          +action(Ax,Phi4x*U*Ay',4)+action(Ax,Phi5x*U,5)+action(Ay,Phi5y*U',5)';
        end


        %Solve the sylvester equations
%         for j=1:order
%             W=sylvester(Ax,Ay',W);
%         end
%         W0kron=reshape(W,m*p,1);

        %Solve the full 2D system
        W0kron=reshape(W,m*p,1);
        for j=1:order
           W0kron=A\W0kron;
        end

        %Time for the Kronecker product
        times_kron(r-r_in+1)=toc;

        %Relative error
        error(r-r_in+1)=norm(V0-W0kron)/norm(V0);

    end

end

%Export the times and errors to a .txt file 
saving=times_full./times_kron;
Tab_times=table(size_A,times_full,times_kron,saving)
writetable(Tab_times,'TimesLaplace_p1_16_diag.txt','Delimiter',' ');

Tab_error=table(size_A,error)
writetable(Tab_error,'RelativeErrorLaplace_p1_16_diag.txt','Delimiter',' ');

    






