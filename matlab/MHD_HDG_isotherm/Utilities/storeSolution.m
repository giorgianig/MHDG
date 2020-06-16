function storeSolution(X,T,Xp,Tp,L,u,p,refElv,boundary)

% to plot vectors and pressure
% plotVelocity(X,T,u,5)
% plotPressure(Xp,Tp,p)
% figure(1), clf,
% plotDiscontinuousSolution(X,T,sqrt(u(1:2:end-1).^2 + u(2:2:end).^2),refElv);
% caxis([0 1])
% 
% figure(2), clf,
% plotDiscontinuousSolution(X,T,L(2:4:end-2) - L(3:4:end-1),refElv);
% caxis([-15 15])
% to plot the streamlines
u_cont = createContinuousSolution_vec(X,T,u);
plotStreamlinesWithStreamFunction(X,T,u_cont,refElv,boundary,100)