// Gmsh project created on Fri Jun 14 10:02:26 2019
//+
Point(1) = {0, 0, -0, 1.0};
//+
Point(2) = {0, 1, -0, 1.0};
//+
Point(3) = {1, 0, -0, 1.0};
//+
Point(4) = {-1, 0, -0, 1.0};
//+
Point(5) = { -0.707106781186547, -0.707106781186547, 0, 1.0};
//+
Point(6) = { 0.707106781186547, -0.707106781186547, 0, 1.0};
//+
Point(15) = { -0.530330085889911, -0.530330085889911, 0, 1.0};
//+
Point(16) = { 0.530330085889911, -0.530330085889911, 0, 1.0};
//+
Point(7) = {0, 0.75, -0, 1.0};
//+
Point(8) = {0, -0.75, -0, 1.0};
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
Circle(7) = {9, 1, 15};
//+
Circle(19) = {15, 1, 8};
//+
Circle(20) = {8, 1, 16};
//+
Circle(8) = {16, 1, 10};
//+
Circle(9) = {10, 1, 7};
//+
Line(10) = {5, 15};
//+
Line(11) = {16, 6};
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
Curve Loop(3) = {6, 7, 19, 18, -14, -13, -17};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {8, 9, 17, -16, -15, -18, 20};
//+
Plane Surface(4) = {4};
//+
Physical Curve("IN", 1) = {16, 15, 14, 13};
//+
Physical Curve("OUT", 2) = {1, 2, 4, 5};
//+
Physical Curve("LIM", 3) = {10, 19, 20, 11};
//+
Physical Surface("DOM", 4) = {1, 3, 4, 2};
//+
Characteristic Length {8, 3, 2, 4} = 1.;
//+
Characteristic Length {5, 15, 6, 16, 10, 7, 8, 9} = 1.;
//+

//+
Transfinite Surface {2};
//+
Transfinite Surface {2} = {6, 2, 7, 16};
//+
Transfinite Surface {2} = {16, 6, 2, 7};
//+
Transfinite Curve {9, 8, 9} = 20 Using Progression 1;
//+
Transfinite Curve {5, 4} = 20 Using Progression 1;
//+
Transfinite Curve {12, 12} = 4 Using Progression 1;
//+
Transfinite Curve {11} = 4 Using Progression 1;
