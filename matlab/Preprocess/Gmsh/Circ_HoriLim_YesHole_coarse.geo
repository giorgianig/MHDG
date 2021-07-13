// Mesh refinement
lc0 = 1.0;
lc1 = 0.1;
lc2 = 0.033;

// Center
Point(1)  = { 0                , 0                , 0, lc0};

// Interior circle (core)
Point(2)  = { 0.5              , 0                , 0, lc1};
Point(3)  = { 0                , 0.5              , 0, lc1};
Point(4)  = {-0.5              , 0                , 0, lc1};
Point(5)  = { 0                ,-0.5              , 0, lc1};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Separatrix (just for though, result in bad distortion)
//Point(6)  = { 0.75             , 0                , 0, lc2};
//Point(7)  = { 0                , 0.75             , 0, lc2};
//Point(8)  = {-0.75             , 0                , 0, lc2};
//Point(9)  = { 0                ,-0.75             , 0, lc2};
//Circle(5) = {6, 1, 7};
//Circle(6) = {7, 1, 8};
//Circle(7) = {8, 1, 9};
//Circle(8) = {9, 1, 6};

// Exterior circle (Wall)
Point(10) = { 0.6614           ,-0.75             , 0, lc2};
Point(11) = { 1                , 0                , 0, lc2};
Point(12) = { 0                , 1                , 0, lc2};
Point(13) = {-1                , 0                , 0, lc2};
Point(14) = {-0.6614           ,-0.75             , 0, lc2};
Circle(9) = {10, 1, 11};
Circle(10)= {11, 1, 12};
Circle(11)= {12, 1, 13};
Circle(12)= {13, 1, 14};
BSpline(13)= {14,10};

// Surface (egde)
//Curve Loop(14) = {1, 2, 3, 4};
//Curve Loop(15) = {5, 6, 7, 8};
//Plane Surface(16) = {14, 15};

// Surface (SOL)
//Curve Loop(17) = {5, 6, 7, 8};
//Curve Loop(18) = {9, 10, 11, 12, 13};
//Plane Surface(19) = {17, 18};

// Surface (all) (better distortion) Check orientation (normal toward you)
Curve Loop(14) = {1, 2, 3, 4};
Curve Loop(15) = {9, 10, 11, 12, 13};
Plane Surface(16) = {14, -15};

Physical Curve("IN", 1) = {1,2,3,4};
Physical Curve("OUT", 2) = {9,10,11,12};
Physical Curve("LIM", 3) = {13};

Physical Surface("DOM", 4) = {16};
