function res = faceNodes_aux(nOfElementNodes)
switch nOfElementNodes
    case 3 %P1
        res = [1 2; 2 3; 3 1];
    case 6 %P2
        res = [1 4 2; 2 5 3; 3 6 1];
    case 10 %P3
        res = [1 4 5 2; 2 6 7 3; 3 8 9 1];
    case 15 %P4
        res = [1 4 5 6 2; 2 7 8 9 3; 3 10 11 12 1];
    case 21 %P5
        res = [1 4:7 2; 2 8:11 3; 3 12:15 1];
    case 28 %P6
        res = [1 4:8 2; 2 9:13 3; 3 14:18 1];
    case 36 %P7
        res = [1 4:9 2; 2 10:15 3; 3 16:21 1];
    case 45 %P8
        res = [1 4:10 2; 2 11:17 3; 3 18:24 1];
    case 55 %P9
        res = [1 4:11 2; 2 12:19 3; 3 20:27 1];
    case 66 %P10
        res = [1 4:12 2; 2 13:21 3; 3 22:30 1];
    case 78 %P11
        res = [1 4:13 2; 2 14:23 3; 3 24:33 1];
end
end