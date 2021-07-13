// Gmsh project created on Fri Jun 14 10:02:26 2019
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {-1, 0, 0, 1.0};
//+
Point(5) = {-0.0001, -1., 0, 1.0};
//+
Point(6) = {0.0001, -1., 0, 1.0};
//+
Point(7) = {0, 0.75, 0, 1.0};
//+
Point(8) = {0, -0.75, 0, 1.0};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {4, 1, 5};
//+
Circle(4) = {6, 1, 3};
//+
Circle(5) = {3, 1, 2};
//+
Point(9) = {-0.75, 0, 0, 1.0};
//+
Point(10) = {0.75, 0, 0, 1.0};
//+
Circle(6) = {7, 1, 9};
//+
Circle(7) = {9, 1, 8};
//+
Circle(8) = {8, 1, 10};
//+
Circle(9) = {10, 1, 7};
//+
Line(10) = {5, 8};
//+
Line(11) = {8, 6};
//+
Line(12) = {2, 7};
//+
Curve Loop(1) = {1, 2, 10, -7, -6, -12};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-9, -8, 11, 4, 5, 12};
//+
Plane Surface(2) = {2};
//+
Point(11) = {0.5, 0, 0, 1.0};
//+
Point(12) = {-0.5, 0, 0, 1.0};
//+
Point(13) = {0, 0.5, 0, 1.0};
//+
Point(14) = {0, -0.5, 0, 1.0};
//+
Circle(13) = {13, 1, 12};
//+
Circle(14) = {12, 1, 14};
//+
Circle(15) = {14, 1, 11};
//+
Circle(16) = {11, 1, 13};
//+
Line(17) = {7, 13};
//+
Line(18) = {8, 14};
//+
Curve Loop(3) = {6, 7, 18, -14, -13, -17};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, 17, -16, -15, -18};
//+
Plane Surface(4) = {4};
//+
Physical Curve("IN", 1) = {16, 15, 14, 13};
//+
Physical Curve("OUT", 2) = {1, 2, 4, 5};
//+
Physical Curve("LIM", 3) = {10, 11};
//+
Physical Surface("DOM", 4) = {1, 3, 4, 2};
//+
Characteristic Length {8} = 1.;
//+
Characteristic Length {5, 6} = 1.;
//+
Characteristic Length {7, 9, 10} = 1;
//+
Characteristic Length {11, 12, 13, 14} = 1;
//+
//Recombine Surface {3, 1, 2, 4};
//+
