% save data & solution for Navier Stokes problem

% mesh
data.mesh.X = X;
data.mesh.T = T;
data.mesh.referenceElement_v = refElv;
data.mesh.Xp = Xp;
data.mesh.Tp = Tp;
data.mesh.referenceElement_p = refElp;
data.mesh.F = F;
data.mesh.flipFace = flipFace;
data.mesh.infoFaces = infoFaces;

% simulation parameters
data.simul_param.Re = Re;
data.simul_param.dt0 = dt0;
data.simul_param.baseIncrDt = baseIncrDt;

% solution
data.solution.velocity = u;
data.solution.pressure = p;
data.solution.time = time;

save(['Saves/data_Re' num2str(Re)] ,'data')
