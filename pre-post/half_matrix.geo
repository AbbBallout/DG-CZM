// Gmsh project created on Sat Oct 22 18:32:58 2022
len=0.5;  
r=len/4;


Point(0) = {0,0,0,1.0};  
Point(1) = {r,0,0,1.0};
Point(2) = {len,0,0,1.0};
Point(3) = {len,len,0,1.0};
Point(4) = {0,len,0,1.0};
Point(5) = {0,r,0,1.0};


Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 0};

Circle(7) = {1, 0, 5};


Curve Loop(1) = {2, 3, 4, 5,-7};
Curve Loop(2) = {1,7,6};

Plane Surface(1) = {1};
Physical Surface("1") = {1};
Plane Surface(2) = {2};
Physical Surface("2") = {2};


Physical Curve(0) = {5,6};
Physical Curve(1) = {3};
Physical Curve(2) = {1,2};
Physical Curve(3) = {4};

Mesh.Algorithm = 9;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 10;
Show "*";



