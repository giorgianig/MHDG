%% aux plot convergence 
close all
% load 'Results/Conv_p_Np5Testcase25_Wx8_Wy8.mat'
load 'Results/Auxconv.mat'

% c = {'b', 'r', 'k', 'm','g'};
c = {'k', 'k', 'k', 'k','k'};
m = {'o-','>-','s-','<-','*-'};
% figure
% for i = 1:size(L2err_dens,2)
%     slope = diff(log10(L2err_dens(:,i))) ./ diff(log10(elsize(:,i)));
%     loglog(elsize(:,i),L2err_dens(:,i),[c{i} m{i}])
%     leg{i} = ['p=' num2str(i)];
%     % print slopes
%     for j=1:numel(slope)
%         text(elsize(j+1,i),L2err_dens(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
%             'fontname','times new roman','color',c{i})
%     end
%     hold on
% end
% ylim([1e-12 1e0])
% legend(leg,'location','southeast','fontname','times new roman')
% title('$n$','fontname','times new roman','interpreter','latex')
% xlabel('Element size','fontname','times new roman')
% ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
% grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'Convergence_density')

% 
% 
% figure
% for i = 1:size(L2err_pvel,2)
%     slope = diff(log10(L2err_pvel(:,i))) ./ diff(log10(elsize(:,i)));
%     loglog(elsize(:,i),L2err_pvel(:,i),[c{i} m{i}])
%     leg{i} = ['p=' num2str(i)];
%     % print slopes
%     for j=1:numel(slope)
%         text(elsize(j+1,i),L2err_pvel(j+1,i)*(1-0.5),num2str(slope(j),'%1.1f'),...
%             'fontname','times new roman','color',c{i})
%     end
%     hold on
% end
% ylim([1e-12 1e0])
% legend(leg,'location','southeast','fontname','times new roman')
% title('$\Gamma$','fontname','times new roman','interpreter','latex')
% xlabel('Element size','fontname','times new roman')
% ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
% grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'Convergence_Gamma')
% 
figure
semilogy(1:size(L2err_dens,1),L2err_dens(5,:),'k*-')
hold on
semilogy(1:size(L2err_pvel,1),L2err_pvel(5,:),'ko-')
xlabel('$p$','fontname','times new roman','interpreter','latex')
ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
legend({'$n$','$\Gamma$'},'location','northeast','fontname','times new roman','interpreter','latex')
grid on
readyforprintnew([8 6],24,[],[],1,[],[],'Pconvergence_Circle')



% p convergence density
% np = size(L2err_dens,1);
% semilogy(1:np,L2err_dens(:,1),[c{1} '*-'])
% hold on 
% semilogy(1:np,L2err_pvel(:,1),[c{2} 'o-'])
% legend({'$n$','$\Gamma$'},'location','northeast','fontname','times new roman','interpreter','latex')
% % title('Density','fontname','times new roman')
% xlabel('$p$','fontname','times new roman','interpreter','latex')
% ylabel('$\mathcal{L}^2$ error','fontname','times new roman','interpreter','latex')
% set(gca,'xtick',1:np ,'fontname','times new roman')
% grid on
% readyforprintnew([8 6],24,[],[],1,[],[],'PconvWest')
