function up = cons2phys(uc)
% uc = np x neq
% up = np x n phys var

global Mref prolongateExponents neq

% mi = 3.35e-27;       % Ionic mass (kg)
% kb = 1.38e-23;       % Boltzmann constant
if neq ==2
    up = zeros(size(uc,1),2);
    up(:,1) = uc(:,1);                                           % density
    up(:,2) = uc(:,2)./(uc(:,1)*sqrt(2*Mref));                   % Mach
else
    up = zeros(size(uc,1),10);
    % Adimensional physical variables
    up(:,1) = uc(:,1);                                           % density
    up(:,2) = uc(:,2)./uc(:,1);                                  % u parallel
    up(:,3) = uc(:,3)./uc(:,1);                                  % total energy for ions
    up(:,4) = uc(:,4)./uc(:,1);                                  % total energy for electrons
    up(:,5) = (2/(3*Mref)*(uc(:,3)-(0.5*uc(:,2).^2)./uc(:,1)));    % pressure for ions
    up(:,6) = (2/(3*Mref)*uc(:,4));                              % pressure for electrons
    up(:,7) = up(:,5)./up(:,1);                                  % temperature of ions
    up(:,8) = up(:,6)./up(:,1);                                  % temperature of electrons
    % if ~prolongateExponents
    up(:,9)   = sqrt((abs(up(:,7))+abs(up(:,8)))*Mref);          % sound speed
    % else
    %     up(:,9)   = prolExp((up(:,7)+up(:,8))*Mref, 0.5);      % sound speed
    % end
    up(:,10) = up(:,2)./up(:,9);                                 % Mach
end


% % conservative to physical variables
% function up = cons2phys(uc)
% % uc = np x neq
% % up = np x n phys var
% 
% global Mref prolongateExponents
% 
% % mi = 3.35e-27;       % Ionic mass (kg)
% % kb = 1.38e-23;       % Boltzmann constant
% up = zeros(size(uc,1),7);
% % Adimensional physical variables
% up(:,1) = uc(:,1);                                                  % density
% up(:,2) = uc(:,2)./uc(:,1);                                      % u parallel
% up(:,3) = uc(:,3)./uc(:,1);                                      % total energy for ions
% up(:,4) = uc(:,4)./uc(:,1);                                      % total energy for electrons
% up(:,5) = (2/(3*Mref)*(uc(:,3)-0.5*uc(:,2).^2./uc(:,1)));     % pressure for ions
% up(:,6) = (2/(3*Mref)*uc(:,4));                                % pressure for electrons
% up(:,7) = up(:,5)./up(:,1);                              % temperature of ions
% up(:,8) = up(:,6)./up(:,1);                              % temperature of electrons
% % if ~prolongateExponents
%     up(:,9)   = sqrt((abs(up(:,7))+abs(up(:,8)))*Mref);                 % sound speed
% % else
% %     up(:,9)   = prolExp((up(:,7)+up(:,8))*Mref, 0.5);                 % sound speed
% % end
% up(:,10) = up(:,2)./up(:,9);                                     % Mach
% % conservative to physical variables
