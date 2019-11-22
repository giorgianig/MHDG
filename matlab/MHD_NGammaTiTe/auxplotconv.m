%% aux plot convergence 
close all
% load 'Results/Auxconv_Nc3_Np3Testcase9_Wx1_Wy1.mat'
load 'Results/Auxconv.mat'

% c = {'b', 'r', 'k', 'm','g'};
c = {'k', 'k', 'k', 'k','k'};
m = {'o-','>-','s-','<-','*-'};
figure
for i = 1:size(L2err_U1,2)
    slope = diff(log10(L2err_U1(:,i))) ./ diff(log10(elsize(:,i)));
    loglog(elsize(:,i),L2err_U1(:,i),[c{i} m{i}])
    leg{i} = ['p=' num2str(i)];
    % print slopes
    for j=1:numel(slope)
        text(elsize(j+1,i),L2err_U1(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
            'fontname','times new roman','color',c{i})
    end
    hold on
end
% ylim([1e-12 1e0])
legend(leg,'location','southeast','fontname','times new roman','interpreter','latex')
title('$n$','fontname','times new roman','interpreter','latex')
xlabel('Element size','fontname','times new roman','interpreter','latex')
ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'ConvDens')

% 
% 
figure
for i = 1:size(L2err_U2,2)
    slope = diff(log10(L2err_U2(:,i))) ./ diff(log10(elsize(:,i)));
    loglog(elsize(:,i),L2err_U2(:,i),[c{i} m{i}])
    leg{i} = ['p=' num2str(i)];
    % print slopes
    for j=1:numel(slope)
        text(elsize(j+1,i),L2err_U2(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
            'fontname','times new roman','color',c{i})
    end
    hold on
end
% ylim([1e-12 1e0])
legend(leg,'location','southeast','fontname','times new roman','interpreter','latex')
title('$nu$','fontname','times new roman','interpreter','latex')
xlabel('Element size','fontname','times new roman','interpreter','latex')
ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'ConvMome')

figure
for i = 1:size(L2err_U3,2)
    slope = diff(log10(L2err_U3(:,i))) ./ diff(log10(elsize(:,i)));
    loglog(elsize(:,i),L2err_U3(:,i),[c{i} m{i}])
    leg{i} = ['p=' num2str(i)];
    % print slopes
    for j=1:numel(slope)
        text(elsize(j+1,i),L2err_U3(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
            'fontname','times new roman','color',c{i})
    end
    hold on
end
% ylim([1e-12 1e0])
legend(leg,'location','southeast','fontname','times new roman','interpreter','latex')
title('$nE_i$','fontname','times new roman','interpreter','latex')
xlabel('Element size','fontname','times new roman','interpreter','latex')
ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'ConvEneI')




figure
for i = 1:size(L2err_U4,2)
    slope = diff(log10(L2err_U4(:,i))) ./ diff(log10(elsize(:,i)));
    loglog(elsize(:,i),L2err_U4(:,i),[c{i} m{i}])
    leg{i} = ['p=' num2str(i)];
    % print slopes
    for j=1:numel(slope)
        text(elsize(j+1,i),L2err_U4(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
            'fontname','times new roman','color',c{i})
    end
    hold on
end
% ylim([1e-12 1e0])
legend(leg,'location','southeast','fontname','times new roman','interpreter','latex')
title('$nE_e$','fontname','times new roman','interpreter','latex')
xlabel('Element size','fontname','times new roman','interpreter','latex')
ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'ConvEneE')



%% Convergence p
% readyforprintnew([8 6],24,[],[],1,[],[],'Convergence_pvel')
% 
% figure
% semilogy(1:size(L2err_U1,1),L2err_U1(5,:),'k*-')
% hold on
% semilogy(1:size(L2err_U2,1),L2err_U2(5,:),'ko-')
% xlabel('Polynomial degree','fontname','times new roman')
% ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
% legend('Density','Parallel velocity')
% grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'Pconvergence_Circle')
