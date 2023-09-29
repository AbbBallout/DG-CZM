// Gmsh project created on Sat Oct 22 18:32:58 2022
len=1;
c=len/2;  
r=len/4;


Point(0) = {0,0,0,10.0};  
Point(1) = {len,0,0,10.0};
Point(2) = {len,len,0,10.0};
Point(3) = {0,len,0,10.0};
Point(4) = {c,c,0,1.0,1.0};
Point(5) = {c,r,0,1.0,1.0};
Point(6) = {c,c+r,0,1.0,1};

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 0};


Circle(5) = {5, 4, 6};
Circle(6) = {6, 4, 5};

Curve Loop(1) = {1, 2, 3, 4};

Curve Loop(2) = {5,6};
Plane Surface(1) = {1,2};
Physical Surface("1") = {1};
Plane Surface(2) = {2};
Physical Surface("2") = {2};


Physical Curve(0) = {4};
Physical Curve(1) = {2};
Physical Curve(2) = {1};
Physical Curve(3) = {3};

Mesh.Algorithm = 9;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 0.05;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 1;
Show "*";



