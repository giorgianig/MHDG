function plotSolutionOnALine(X,T,u,refEl,func,h)

nel = size(T,1);

for iel = 1:nel
    
   Te = T(iel,:);
   Xe = X(Te,:);
   x = X(Te,1);y = X(Te,2);
  if checkcrossing(func(x,y))
    
      
      [xp,yp,zp] = computeElement(Xe,refEl,func,h);
      
      
  end
end


function [xp,yp,up] = computeElement(Xe,refEl,func,h)

N1d = refEl.N1d;
nf = size(refEl.faceNodes,1); 

% number of crossing for each face
ncf = zeros(nf,1);
xycross = zeros(6,2);

% find the crossing in each face
for ifac = 1:nf 
    
    xg = [Xe(refEl.faceNodes(ifac,end),1);...
        N1d*Xe(refEl.faceNodes(ifac,:),1);...
        Xe(refEl.faceNodes(ifac,1),1)];
    yg = [Xe(refEl.faceNodes(ifac,end),2);...
        N1d*Xe(refEl.faceNodes(ifac,:),2);...
        Xe(refEl.faceNodes(ifac,1),2)];
    
    if checkcrossing(func(xg,yg))
        % there is at least one crossing in this face
        
        % check the number of crossing
        sn = func(xg(1),yg(1));
        for ig = 2:numel(xg)
           if sn*func(xg(ig),yg(ig))<0
               
               xy = findzero([xg(ig-1),yg(ig-1)],[xg(ig),yg(ig)],func);
               ncf(ifac) = ncf(ifac)+1;
               xycross(sum(ncf),:) = xy;
               sn = func(xg(ig),yg(ig));
           end
        end
    end
end

% Only one crossing per face is allowed (so far)
if any(ncf>1),error('More then one crossing per face!'), end

xycross = xycross(1:ncf,:);

% find points
for ifac = 1:nf-1
   
    if ~ncf(ifac)
        continue
    end
    xy1 = xycross(ifac,:);
    
    % distance 
     for ifacn = ifac:nf
        
         if ~ncf(ifacn)
             continue
         end
         xy2 = xycross(ifacn,:);  
     end
     
     dist = norm(xy2-xy1);
     
     np = ceil(dist/h);
     
     if abs(xy2(1)-xy1(1))>abs(xy2(2)-xy1(2))
         % use x
         xp = linspace(xy1(1),xy2(1),np);
         ypst = linspace(xy1(2),xy2(2),np);
         yp = zeros(size(ypst));
         del = abs(ypst(2)-ypst(1));
         for ip = 1:np
             yp(ip) = findcoord(func,xp(ip),ypst(ip),1,del);
         end
     else
         % use y
         xpst = linspace(xy1(1),xy2(1),np);
         yp = linspace(xy1(2),xy2(2),np);
         xp = zeros(size(xpst));
         del = abs(xpst(2)-xpst(1));
         for ip = 1:np
             xp(ip) = findcoord(func,xpst(ip),yp(ip),2,del);
         end
     end
     
     up = computeValues(xp,yp,Xe,refEl);

end






function res = checkcrossing(f)

res = all(any([f>0,f<0],1));


function xy = findzero(xy1,xy2,func)

tol = 1e-8;
maxit = 1000;
check = 1;
% little check
if ~func(xy1(1),xy1(2))*func(xy2(1),xy2(2))<0
    error('Something wrong')
end

% make sure f(xy1)>0
if func(xy1(1),xy1(2))<0
    temp = xy1;
    xy1 = xy2;
    xy2 = temp;
end

for i = 1:maxit 
   xy = 0.5*(xy1+xy2);
   
   err = func(xy(1),xy(2));
   if abs(err)<tol
       check = 0;
       break
   else 
       if err>0
           xy1 = xy;
       else
           xy2 = xy;
       end
   end
end
if check, error('Not converging'), end



function res = findcoord(f,x,y,swit,del)
    
maxit = 1000;
tol = 1e-8;

if swit==1
    %% x fixed, search for y
    
    % first check
    if abs(f(x,y))<tol
        res = y;
        return
    else
        % search next y
        yn = nan;
        for i = 1:maxit
            if f(x,y)*f(x,y+i*del)<0
                yn = y+i*del;
                break
            elseif f(x,y)*f(x,y-i*del)<0
                yn = y-i*del;
                break
            end
        end
        if isnan(yn),error('Next y not found'),end
        if f(x,y)>0
            y1 = y;
            y2 = yn;
        else
            y1 = yn;
            y2 = y;
        end        
    end
    check = 1;
    for i = 1:maxit
       
        y = 0.5*(y1+y2);
        err = abs(f(x,y));
        if err<tol
            res = y;
            check = 0;
            break
        else
            if err>0
                y1 = y;
            else
                y2 = y;
            end
        end
    end
    if check, error('Not converging'), end
elseif swit ==2
%% y fixed, search for x
        % first check
    if abs(f(x,y))<tol
        res = x;
        return
    else
        % search next x
        xn = nan;
        for i = 1:maxit
            if f(x,y)*f(x+i*del,y)<0
                xn = x+i*del;
                break
            elseif f(x,y)*f(x-i*del,y)<0
                xn = x-i*del;
                break
            end
        end
        if isnan(xn),error('Next x not found'),end
        if f(x,y)>0
            x1 = x;
            x2 = xn;
        else
            x1 = xn;
            x2 = x;
        end           
    end
    
    check = 1;
    for i = 1:maxit
       
        x = 0.5*(x1+x2);
        err = abs(f(x,y));
        if err<tol
            res = x;
            check = 0;
            break
        else
            if err>0
                x1 = x;
            else
                x2 = x;
            end
        end
    end 
    if check, error('Not converging'), end
    
else
   error('Wrong switch') 
end
    

function u = computeValues(x,y,X,refEl)








    
    