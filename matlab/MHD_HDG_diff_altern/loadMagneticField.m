function loadMagneticField(X,T,F,refEl)
%% Load magnetic field from external file and store it in each Gauss point
%% in the elements and in the faces

global Magnetic axisym Mesh driftvel testcase

%% Load file
magfieldwest = load('/home/giorgio/Matlab/Meshes/Meshes_2D/WEST/WEST_far_465.mat');

magfieldwest.r2D = magfieldwest.r2D/Mesh.lscale;
magfieldwest.z2D = magfieldwest.z2D/Mesh.lscale;

%% Store in the elements Gauss points
xgauss = refEl.N*reshape(X(T',1),size(T'));
ygauss = refEl.N*reshape(X(T',2),size(T'));
bx = -magfieldwest.Br2D./sqrt(magfieldwest.Br2D.^2+magfieldwest.Bz2D.^2+magfieldwest.Bphi2D.^2);
by = -magfieldwest.Bz2D./sqrt(magfieldwest.Br2D.^2+magfieldwest.Bz2D.^2+magfieldwest.Bphi2D.^2);

Magnetic.bxnodes = interp2(magfieldwest.r2D,magfieldwest.z2D,bx,reshape(X(T',1),size(T')),reshape(X(T',2),size(T')));
Magnetic.bynodes = interp2(magfieldwest.r2D,magfieldwest.z2D,by,reshape(X(T',1),size(T')),reshape(X(T',2),size(T')));
Magnetic.bx = interp2(magfieldwest.r2D,magfieldwest.z2D,bx,xgauss,ygauss);
Magnetic.by = interp2(magfieldwest.r2D,magfieldwest.z2D,by,xgauss,ygauss);
if testcase.n>50
    Magnetic.flux2D = interp2(magfieldwest.r2D,magfieldwest.z2D,magfieldwest.flux2D,xgauss,ygauss);
end
div = zeros(size(bx));

if axisym
    div(2:end-1,2:end-1) = 1./magfieldwest.r2D(2:end-1,2:end-1).* (magfieldwest.r2D(2:end-1,3:end).*bx(2:end-1,3:end)- ...
        magfieldwest.r2D(2:end-1,1:end-2).*bx(2:end-1,1:end-2))./...
        (magfieldwest.r2D(2:end-1,3:end)-magfieldwest.r2D(2:end-1,1:end-2))+...
        (by(3:end,2:end-1)-by(1:end-2,2:end-1))./...
        (magfieldwest.z2D(3:end,2:end-1)-magfieldwest.z2D(1:end-2,2:end-1));
    
else
    div(2:end-1,2:end-1) = (bx(2:end-1,3:end)-bx(2:end-1,1:end-2))./...
        (magfieldwest.r2D(2:end-1,3:end)-magfieldwest.r2D(2:end-1,1:end-2))+...
        (by(3:end,2:end-1)-by(1:end-2,2:end-1))./...
        (magfieldwest.z2D(3:end,2:end-1)-magfieldwest.z2D(1:end-2,2:end-1));
end
Magnetic.div = interp2(magfieldwest.r2D,magfieldwest.z2D,div,xgauss,ygauss);



if driftvel
    dvx =  zeros(size(bx));
    dvy =  zeros(size(bx));

    Bmod = sqrt(magfieldwest.Br2D.^2+magfieldwest.Bz2D.^2+magfieldwest.Bphi2D.^2);
    Bmodx = (Bmod(2:end-1,3:end)-Bmod(2:end-1,1:end-2))./...
        (magfieldwest.r2D(2:end-1,3:end)-magfieldwest.r2D(2:end-1,1:end-2));
    Bmody = (Bmod(3:end,2:end-1)-Bmod(1:end-2,2:end-1))./...
        (magfieldwest.z2D(3:end,2:end-1)-magfieldwest.z2D(1:end-2,2:end-1));
    dvx(2:end-1,2:end-1) =  -magfieldwest.Bphi2D(2:end-1,2:end-1).*Bmody./Bmod(2:end-1,2:end-1).^3;
    dvy(2:end-1,2:end-1) =   magfieldwest.Bphi2D(2:end-1,2:end-1).*Bmodx./Bmod(2:end-1,2:end-1).^3;
    Magnetic.dvxnodes = interp2(magfieldwest.r2D,magfieldwest.z2D,dvx,reshape(X(T',1),size(T')),reshape(X(T',2),size(T')));
    Magnetic.dvynodes = interp2(magfieldwest.r2D,magfieldwest.z2D,dvy,reshape(X(T',1),size(T')),reshape(X(T',2),size(T')));
    Magnetic.dvx = interp2(magfieldwest.r2D,magfieldwest.z2D,dvx,xgauss,ygauss);
    Magnetic.dvy = interp2(magfieldwest.r2D,magfieldwest.z2D,dvy,xgauss,ygauss);
end
%% Store in the faces Gauss points
nFac = max(max(F));
nf = size(refEl.faceNodes,1);
check = false(nFac,1);
Magnetic.bxfaces = zeros(numel(refEl.IPweights1d),nFac);
Magnetic.byfaces = zeros(numel(refEl.IPweights1d),nFac);
for iel = 1:size(F,1)
    
    for ifac = 1:nf
        Face = F(iel,ifac);
        if check(Face)
            continue
        end
        Tf = T(iel,refEl.faceNodes(ifac,:));
        xgauss = refEl.N1d*X(Tf,1);
        ygauss = refEl.N1d*X(Tf,2);
        Magnetic.bxfaces(:,Face) = interp2(magfieldwest.r2D,magfieldwest.z2D,bx,xgauss,ygauss);
        Magnetic.byfaces(:,Face) = interp2(magfieldwest.r2D,magfieldwest.z2D,by,xgauss,ygauss);
        check(Face) = true;
    end
end

