n = 1000;
ys=linspace(-3,3,n+1)'; xs=7.7*ones(size(ys)); 

load XHDG_soln
q1=q(1:2:end); q2=q(2:2:end);
Q = [q1,q2];
QlineXHDG=evalDGapproximationAtPoints([xs,ys],Q,X,T,referenceElement);

Qline = QlineXHDG;
% figure(1), clf, plot(ys,Qline,'LineWidth',2), legend('q1','q2')
% title('X-HDG velocity components along the line x=8','FontSize',14)
I1_XHDG=trapz(ys,Qline(:,1)), I2_XHDG=trapz(ys,Qline(:,2))

load HDGCoarse_soln
q1=q(1:2:end); q2=q(2:2:end);
Q = [q1,q2];
QlineHDG=evalDGapproximationAtPoints([xs,ys],Q,X,T,referenceElement);

Qline = QlineHDG;
% figure(2), clf, plot(ys,Qline,'LineWidth',2), legend('q1','q2')
% title('HDG velocity components along the line x=8','FontSize',14)
I1_HDG=trapz(ys,Qline(:,1)), I2_HDG=trapz(ys,Qline(:,2))

load HDGFine_soln
q1=q(1:2:end); q2=q(2:2:end);
Q = [q1,q2];
QlineHDGFine=evalDGapproximationAtPoints([xs,ys],Q,X,T,referenceElement);

Qline = QlineHDGFine;
% figure(3), clf, plot(ys,Qline,'LineWidth',2), legend('q1','q2')
% title('HDG (fine) velocity components along the line x=8','FontSize',14)
I1_HDGFine=trapz(ys,Qline(:,1)), I2_HDGFine=trapz(ys,Qline(:,2))

figure(1), clf
plot(ys,QlineHDGFine(:,1),'k-',ys,QlineHDG(:,1),'r-',ys,QlineXHDG(:,1),'b-')
xlabel('y'), ylabel('Velocity - HORIZONTAL direction'), legend('reference','HDG','XHDG')

figure(2), clf
plot(ys,QlineHDGFine(:,2),'k-',ys,QlineHDG(:,2),'r-',ys,QlineXHDG(:,2),'b-')
xlabel('y'), ylabel('Velocity - VERTIVCAL direction'), legend('reference','HDG','XHDG')

figure(3), clf
plot(ys,QlineHDG(:,1)-QlineHDGFine(:,1),'r-',ys,QlineXHDG(:,1)-QlineHDGFine(:,1),'b-')
xlabel('y'), ylabel('Velocity ERROR - HORIZONTAL direction'), legend('HDG','XHDG')

figure(4), clf
plot(ys,QlineHDG(:,2)-QlineHDGFine(:,2),'r-',ys,QlineXHDG(:,2)-QlineHDGFine(:,2),'b-')
xlabel('y'), ylabel('Velocity ERROR - VERTICAL direction'), legend('HDG','XHDG')

Error_normalFlux_HDG =  I1_HDG-I1_HDGFine
Error_normalFlux_HDG =  I1_XHDG-I1_HDGFine
