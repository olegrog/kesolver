Kn = 0.1;   // Knudsen number
N = 50;     // Total number of cells

L = 1./Kn;  // The length is measured in terms of the mean free path

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};

Line(1) = {1, 2};           // Create a line connecting points 1 and 2
Transfinite Curve{1} = N;   // Uniformly place N nodes on this line

// Convert 1D geometry to 2D
Extrude {0, 1, 0} {
	Line{1};
	Layers{1};
	Recombine;
}

// Convert 2D geometry to 3D
Extrude {0, 0, 1} {
	Surface{5}; // Hereinafter, the number of the surface can be found in GMSH GUI
	Layers{1};
	Recombine;
}

// We have to name some physical surfaces to specify the BC
Physical Surface("left") = {26};
Physical Surface("right") = {18};
// We have to name the whole volume to specify the IC
Physical Volume("volume") = {1};
