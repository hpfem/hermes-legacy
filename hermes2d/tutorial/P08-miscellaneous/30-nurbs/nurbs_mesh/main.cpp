#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example illustrates how to use full-featured NURBS. Simplified
// format is enabled for circular arcs (see example 03-poisson). 
//
// Domain: Rectangle (0, 2) x (0, 1) where the upper edge is a NURBS
//         (see the end of the mesh file for details).
//
// Choose one of the following mesh files:

const int INIT_REF_NUM = 2;

//const char* mesh_file = "domain-1.mesh";            // One control point.
//const char* mesh_file = "domain-2.mesh";          // Two control points.
const char* mesh_file = "domain-4.mesh";          // Three control points.

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load(mesh_file, &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  MeshView mview("Nurbs", new WinGeom(0, 0, 350, 350));
  mview.show(&mesh);
	
  // Wait for the view to be closed.
  View::wait();
  return 0;
}

