%%%%Script for plotting the solution and the errors of the%%%
%%%Allen-Cahn equation%%%

%%%Plot a particular time step%%%
i=1;
Umat=reshape(uhat_kron(:,i),[nelx+1,nely+1])';
figure
h=pcolor(xsol,ysol,Umat);
set(h, 'EdgeColor', 'none')
%mesh(xmat,ymat,Umat);
colorbar, set(colorbar,'TickLabelInterpreter','latex')
xlabel({'$x$'},'interpreter','latex')
ylabel({'$y$'},'interpreter','latex')
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
title({'$t=$'},'interpreter','latex')
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)
caxis([-1 1])


%%%Movie of the whole solution%%%
figure
hold on
title('Allen-Cahn equation','interpreter','latex')
xlabel({'$x$'},'interpreter','latex')
ylabel({'$y$'},'interpreter','latex')
xticks([0 0.25 0.5 0.75 1]);
yticks([0 0.25 0.5 0.75 1]);
set(gca,'TickLabelInterpreter','latex'),set(gca,'fontsize',16)
v=[1:50 51:20:steps+1];
for j=1:size(v,2)
    i=v(j);
    Umat=reshape(uhat_2D(:,i),[nelx+1,nely+1])';
    
    h=pcolor(xsol,ysol,Umat);
    set(h, 'EdgeColor', 'none')
    %mesh(xmat,ymat,Umat);
    colorbar, set(colorbar,'TickLabelInterpreter','latex')
    caxis([-1 1])
    pause(0.01)
end

