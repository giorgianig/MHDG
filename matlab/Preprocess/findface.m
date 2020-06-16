function f = findface(s)
% find face number

if numel(s)==3
    if isequal(s, [1 1 0])
        f = 1;
    elseif isequal(s,[1 0 1])
        f = 2;
    elseif isequal(s, [0 1 1])
        f = 3;
    else
        error('Something wrong in the input')
    end
elseif numel(s)==4
    if isequal(s,[1 1 0 0])
        f = 1;
    elseif isequal(s, [0 1 1 0])
        f = 2;
    elseif isequal(s,[0 0 1 1])
        f =3;
    elseif isequal(s,[1 0 0 1])
        f = 4;
    else
        error('Something wrong in the input')
    end
else
    error('Case not implemented')
end
