// Gmsh project created on Sat Oct 22 18:32:58 2022
w1 = 100;
l1 = 30;
w2 = 10;
l2 = 40;
in = 20;

//Square
Point(0) = {-w1/2,0,0,1.0};  
Point(1) = {w1/2,0,0,1.0};
Point(2) = {w1/2,l1,0,1.0};
Point(3) = {w2/2,l1,0,1};
Point(4) = {-w2/2,l1,0,1};
Point(5) = {-w1/2,l1,0,1};

//Bar
Point(6) = {w2/2,l1-in,0,1};  
Point(7) = {w2/2,l1+l2-in,0,1};
Point(8) = {-w2/2,l1+l2-in,0,1};
Point(9) = {-w2/2,l1-in,0,1};

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 9};
Line(6) = {9, 4};
Line(7) = {4, 5};
Line(8) = {5, 0};


Line(9) = {3, 7};
Line(10) = {7, 8};
Line(11) = {8, 4};


Curve Loop(1) = {1, 2, 3, 4,5,6,7,8};
Curve Loop(2) = {-4,9,10,11,-6,-5};

Plane Surface(1) = {1};
Physical Surface("1") = {1};
Plane Surface(2) = {2};
Physical Surface("2") = {2};

Physical Curve(1) = {2,8};
Physical Curve(2) = {1};
Physical Curve(3) = {10};


Mesh.Algorithm = 9;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 5;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 1;
Show "*";



