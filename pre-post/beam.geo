// Gmsh project created on Sat Oct 22 18:32:58 2022
w1=1;
w2=0.7;

l1=0.25;
l2=1.5;
epsi=0.05;

Point(0) = {0,0,0,1.0};
Point(1) = {w1,0,0,1.0};
Point(2) = {w1,l1,0,1.0};
Point(3) = {w1-w2/2,l1+epsi,0,1.0};
Point(4) = {w1-w2/2,l1+l2-epsi,0,1.0};
Point(5) = {w1,l1+l2,0,1.0};
Point(6) = {w1,2*l1+l2,0,1.0};
Point(7) = {0,2*l1+l2,0,1.0};
Point(8) = {0,l1+l2,0,1.0};
Point(9) = {w2/2,l1+l2-epsi,0,1.0};
Point(10) = {w2/2,l1+epsi,0,1.0};
Point(11) = {0,l1,0,1.0};

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 0};


Curve Loop(1) = {1, 2, 3, 4,5,6,7,8,9,10,11,12};
Plane Surface(1) = {1};
Physical Surface("1") = {1};

Physical Curve(0) = {8,9,10,11,12};
Physical Curve(1) = {2,3,4,5,6};
Physical Curve(2) = {1};
Physical Curve(3) = {7};



Mesh.Algorithm = 9;
Mesh.RecombineAll = 1;
Mesh.CharacteristicLengthFactor = 1;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 10;
Show "*";




