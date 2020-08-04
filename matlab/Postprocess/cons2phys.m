% conservative to physical variables
function up = cons2phys(uc,simulation_parameters)
% uc = np x neq
% up = np x n phys var



if strcmpi(simulation_parameters.model,'N-Gamma')
    a = simulation_parameters.physics.a;
    up = zeros(size(uc,1),2);
    if simulation_parameters.switches.logrho
        up(:,1) = exp(uc(:,1));
    else
        up(:,1) = uc(:,1);                                        % density
    end
    up(:,2) = uc(:,2)./up(:,1)/sqrt(a);                       % Mach
elseif strcmpi(simulation_parameters.model,'N-Gamma-Ti-Te')
    Mref = simulation_parameters.physics.Mref;
    up = zeros(size(uc,1),10);
    if simulation_parameters.switches.logrho
        up(:,1) = exp(uc(:,1));
    else
        up(:,1) = uc(:,1);                                        % density
    end                                       % density
    up(:,2) = uc(:,2)./up(:,1)/Mref;                          % Mach
    up(:,3) = uc(:,3)./up(:,1);                               % total energy for ions
    up(:,4) = uc(:,4)./up(:,1);                               % total energy for electrons
    up(:,5) = (2/(3*Mref)*(uc(:,3)-0.5*uc(:,2).^2./up(:,1))); % pressure for ions
    up(:,6) = (2/(3*Mref)*uc(:,4));                           % pressure for electrons
    up(:,7) = up(:,5)./up(:,1);                               % temperature of ions
    up(:,8) = up(:,6)./up(:,1);                               % temperature of electrons
    up(:,9)   = sqrt((abs(up(:,7))+abs(up(:,8)))*Mref);       % sound speed
    up(:,10) = up(:,2)./up(:,9);                              % Mach
elseif strcmpi(simulation_parameters.model,'N-Gamma-Vorticity')
    a = simulation_parameters.physics.a;
    up = zeros(size(uc,1),4);
    if simulation_parameters.switches.logrho
        up(:,1) = exp(uc(:,1));
    else
        up(:,1) = uc(:,1);                                        % density
    end                                        % density
    up(:,2) = uc(:,2)./up(:,1)/sqrt(a);                       % Mach
    up(:,3) = uc(:,3) ;                                       % total energy for ions
    up(:,4) = uc(:,4);                                        % total energy for electrons
end
% Adimensional physical variables



