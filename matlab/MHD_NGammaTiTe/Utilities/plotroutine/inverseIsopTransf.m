function xieta = inverseIsopTransf(x,Xe,refEl)

tol = 1.e-10;
maxit = 10;
npoints = size(x,1);
nNod = size(Xe,1);

% First trial assuming straight-sided element
if refEl.elemType==1
    xieta0 = inverseLinearTransformation(x,Xe);
    x0 = isoparametricTransformationHighOrder(xieta0,Xe,refEl.degree,refEl.NodesCoord);
else
    xieta0 = inverseLinearTransformationQua(x,Xe);
    x0 = isoparametricTransformationHighOrderQua(xieta0,Xe,refEl.degree,refEl.NodesCoord1d);    
end
ind = sqrt((x(:,1)-x0(:,1)).^2+(x(:,2)-x0(:,2)).^2)>tol*sqrt(x(:,1).^2+x(:,2).^2)+1.e-14;
if any(ind)
    x = x(ind,:);
    x0 = x0(ind,:);
    xieta = xieta0(ind,:);
    
    for iter = 1:maxit
        
        if all(sqrt((x(:,1)-x0(:,1)).^2+(x(:,2)-x0(:,2)).^2)<tol*sqrt(x(:,1).^2+x(:,2).^2)+1.e-14)
            break
        end
        
        N = zeros(npoints,nNod); Nx = N; Ny = N;
        if refEl.elemType==1
            
            % Compute shape functions at interpolation points
            V = Vandermonde_LP(refEl.degree,refEl.NodesCoord);
            invV = V'\eye(size(V'));
            [~,p_xi,p_eta] = orthopoly2D_deriv_xieta(xieta,refEl.degree);
            Nx = (invV*p_xi)';
            Ny = (invV*p_eta)';
        else
            sp=computeShapeFunctionsAtPointsQuads(refEl.degree,refEl.NodesCoord1d,xieta);
            Nx = sp(:,:,2);
            Ny = sp(:,:,3);
        end
        Jxx = Nx*Xe(:,1);
        Jxy = Ny*Xe(:,1);
        Jyx = Nx*Xe(:,2);
        Jyy = Ny*Xe(:,2);
        detJ = Jxx.*Jyy-Jxy.*Jyx;
        rhs = x-x0;
        xieta(:,1)=xieta(:,1)+( rhs(:,1).*Jyy-rhs(:,2).*Jxy)./detJ;
        xieta(:,2)=xieta(:,2)+( rhs(:,2).*Jxx-rhs(:,1).*Jyx)./detJ;
        
        %         J = [Nx(i,:)*Xe(:,1)  Ny(i,:)*Xe(:,1)
        %              Nx(i,:)*Xe(:,2)  Ny(i,:)*Xe(:,2)];
        %         xieta(i,:) = xieta(i,:) + (x(i,:)-x0(i,:))/J';
        if refEl.elemType==1
            x0 = isoparametricTransformationHighOrder(xieta,Xe,refEl.degree,refEl.NodesCoord);
        else
            x0 = isoparametricTransformationHighOrderQua(xieta,Xe,refEl.degree,refEl.NodesCoord1d);
        end
        
    end
    if any(sqrt((x(:,1)-x0(:,1)).^2+(x(:,2)-x0(:,2)).^2)>tol*sqrt(x(:,1).^2+x(:,2).^2)+1.e-14)
        error('Not converging')
    end
    aux_xieta = xieta0;
    aux_xieta(ind,:) = xieta;
    xieta = aux_xieta;
else
    xieta = xieta0;
end