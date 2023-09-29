// Gmsh project created on Sat Oct 22 18:32:58 2022
len=1; width=0.5;
c1=0.4; c2=0.6;  
r1=0.25; r2=0.25;


Point(0) = {0,0,0,1.0};  
Point(1) = {c1-r1,0,0,1.0};
Point(2) = {c1,0,0,1.0};
Point(3) = {c1+r1,0,0,1.0};
Point(4) = {len,0,0,1.0};
Point(5) = {len,width,0,1.0};
Point(6) = {c2+r2,width,0,1.0};
Point(7) = {c2,width,0,1.0};
Point(8) = {c2-r2,width,0,1.0};
Point(9) = {0,width,0,1.0};


Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 0};


Circle(11) = {3, 2, 1};
Circle(12) = {8, 7, 6};

Curve Loop(1) = {1, -11, 4, 5, 6, -12, 9,10};
Plane Surface(1) = {1};
Physical Surface("1") = {1};


Curve Loop(2) = {2, 3, 11};
Plane Surface(2) = {2};
Physical Surface("2") = {2};

Curve Loop(3) = {12, 7, 8};
Plane Surface(3) = {3};
Physical Surface("3") = {3};


Physical Curve(0) = {9, 8, 6, 7, 4, 3, 2, 1};
Physical Curve(1) = {10};
Physical Curve(2) = {5};

Mesh.Algorithm = 9;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 10;
Show "*";


