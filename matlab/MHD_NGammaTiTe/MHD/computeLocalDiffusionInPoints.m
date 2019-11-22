function [diff,xy_corn] = computeLocalDiffusionInPoints(X)

global localDiffPoints

switch localDiffPoints
    case 1
        % Small limiter
        xcorn =  [1.712987425726;2.287012574274;1000;1000];
        ycorn = [-0.692909649383;-0.692909649383;1000;1000];
        xy_corn = [xcorn';ycorn'];
        
        x = X(:,1);
        y = X(:,2);
        
        r1 = sqrt( (x-xcorn(1)).^2+(y-ycorn(1)).^2 );
        r2 = sqrt( (x-xcorn(2)).^2+(y-ycorn(2)).^2 );
        
        d = 1e-2;
        lref = 1e-3;
%         diff_n = zeros(size(x));

        diff = d* (exp(-(r1/lref).^2) + exp( -(r2/lref).^2) );
        
    case 2
        % Infinitely small limiter
         xcorn =  3.4;
        ycorn = -0.75;
        xy_corn = [xcorn';ycorn'];
        
        x = X(:,1);
        y = X(:,2);
        
        r = sqrt( (x-xcorn).^2+(y-ycorn).^2 );
        
        d = 1e1;
        lref = 1e-2;
        diff = d* exp(-(r/lref).^2);
%         diff_n = zeros(size(x));
%         diff_u = d* exp(-(r/lref).^2);
        
end