
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"

#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"


#include "TPZVTKGeoMesh.h"
#include "pzbndcond.h"


// ******************* (Geometry linear description: Zero and One D element) *****************

void ZeroDElements();
void OneDElements();

// ******************* (Geometry linear description: Two D element) **************************

void TwoDTriElements();
void TwoDQuadElements();

// ******************* (Geometry linear description: Three D element) ************************

void ThreeDTetraElements();
void ThreeDPyraElements();
void ThreeDPrisElements();
void ThreeDHexaElements();

// ******************* (Create linear meshes: 1D) *******************************************

TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el);

TPZGeoMesh *CreateOneDNLGMesh(long num_el, REAL size_el);



// ******************* (Create linear meshes: 2D) *******************************************

TPZGeoMesh *CreateTwoDSimpGMesh(int num_div, REAL Lx, REAL Ly);

TPZGeoMesh *CreateTwoDTriGMesh(int num_div, REAL Lx, REAL Ly);

TPZGeoMesh *CreateTwoDQuadGMesh(int num_div, REAL Lx, REAL Ly);







// ******************* (Create linear meshes: 3D) *******************************************




int main() {

    ZeroDElements();
    OneDElements();
    
    TwoDTriElements();
    TwoDQuadElements();
    
    ThreeDTetraElements();
    ThreeDPyraElements();
    ThreeDPrisElements();
    ThreeDHexaElements();
    
    
    // ******************* (Create linear meshes: 1D) *******************************************
    
    REAL domain = 1.;
    long num_el = 10;
    REAL size_el = domain/num_el;
    
    TPZGeoMesh *gmesh_OneDL = CreateOneDLGMesh(num_el, size_el); // function to create the 1D geometric mesh
    
    std::ofstream outgmeshOneDL("geomesh_OneDL.txt");
    gmesh_OneDL->Print(outgmeshOneDL);
    
    std::ofstream vtkgmeshOneDL("geomesh_OneDL.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDL, vtkgmeshOneDL);

    // -----------------------------------------
    
    
    TPZGeoMesh *gmesh_OneDNL = CreateOneDNLGMesh(num_el, size_el); // function to create the 1D geometric mesh
    
    std::ofstream outgmeshOneDNL("geomesh_OneDNL.txt");
    gmesh_OneDNL->Print(outgmeshOneDNL);
    
    std::ofstream vtkgmeshOneDNL("geomesh_OneDNL.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDNL, vtkgmeshOneDNL);
    
    // ******************* (Create linear meshes: 2D) *******************************************

    long num_divsi=2; // number of divition
    REAL Lx=1.; // length of domain in x direction
    REAL Ly=1.; // length of domain in y direction
    
    TPZGeoMesh *gmesh_TwoDSimp = CreateTwoDSimpGMesh(num_divsi, Lx, Ly); // function to create the 2D geometric mesh
    
    std::ofstream outgmeshTwoDSimp("geomesh_TwoDSimp.txt");
    gmesh_TwoDSimp->Print(outgmeshTwoDSimp);
    
    std::ofstream vtkgmeshTwoDSimp("geomesh_TwoDSimp.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDSimp, vtkgmeshTwoDSimp);
    
    // --------------------------------------------------------
    
    long num_divtri=1; // number of divition

    TPZGeoMesh *gmesh_TwoDTri = CreateTwoDTriGMesh(num_divtri, Lx, Ly); // function to create the 2D geometric mesh
    
    std::ofstream outgmeshTwoDTri("geomesh_TwoDTri.txt");
    gmesh_TwoDTri->Print(outgmeshTwoDTri);
    
    std::ofstream vtkgmeshTwoDTri("geomesh_TwoDTri.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDTri, vtkgmeshTwoDTri);
    
    // --------------------------------------------------------
    
    long num_divtqu=1; // number of divition
    
    TPZGeoMesh *gmesh_TwoDQuad = CreateTwoDQuadGMesh(num_divtqu, Lx, Ly); // function to create the 2D geometric mesh
    
    std::ofstream outgmeshTwoDQuad("geomesh_TwoDQuad.txt");
    gmesh_TwoDQuad->Print(outgmeshTwoDQuad);
    
    std::ofstream vtkgmeshTwoDQuad("geomesh_TwoDQuad.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDQuad, vtkgmeshTwoDQuad);


    
    return 0;
}





void ZeroDElements(){
    
    TPZGeoMesh geometry_1DP; // Create the objet that will describe the geometry (Point).
    int n_nodes = 3; // number of nodes
    int geometry_dim = 1; // geometry dimension
    std::string name("geometry 1DP"); // geometry name
    
    // setting the object
    geometry_1DP.SetName(name);
    geometry_1DP.SetDimension(geometry_dim);
    geometry_1DP.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> point_topology(1);
    TPZVec<REAL> x(3,0.0);
    
    // 1st node is located at x ={-100,PI,0.5}
    {
        x[0] = -100.0; // x coordinate
        x[1] = M_PI; // Y coordinate
        x[2] = 0.5; // Z coordinate
        geometry_1DP.NodeVec()[0].SetNodeId(node_id);
        geometry_1DP.NodeVec()[0].SetCoord(x);
        point_topology[0] = node_id;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
    }
    
    // 2nd node is located at x ={0.0,PI,0.5}
    node_id++;
    element_id++;
    {
        x[0] = 0.0; // x coordinate
        x[1] = M_PI; // Y coordinate
        x[2] = 0.5; // Z coordinate
        geometry_1DP.NodeVec()[1].SetNodeId(node_id);
        geometry_1DP.NodeVec()[1].SetCoord(x);
        point_topology[0] = node_id;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
    }
    
    // 3nd node is located at x ={100,PI,0.5}
    node_id++;
    element_id++;
    {
        x[0] = 100.0; // x coordinate
        x[1] = M_PI; // Y coordinate
        x[2] = 0.5; // Z coordinate
        geometry_1DP.NodeVec()[2].SetNodeId(node_id);
        geometry_1DP.NodeVec()[2].SetCoord(x);
        point_topology[0] = node_id;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
    }
    
    geometry_1DP.BuildConnectivity();
    
    
    std::ofstream file("geometry_1DP.txt");
    geometry_1DP.Print(file);
    
}

void OneDElements(){
    
    TPZGeoMesh geometry_1DL; // Create the objet that will describe the geometry (Line).
    int n_nodes = 2; // number of nodes
    int geometry_dim = 1; // geometry dimension
    std::string name("geometry 1DL"); // geometry name
    
    // setting the object
    geometry_1DL.SetName(name);
    geometry_1DL.SetDimension(geometry_dim);
    geometry_1DL.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Linear_topology(2);
    TPZVec<REAL> coord(3,0.0);
    
    {
    // 1st node of line element is located at x ={-10,PI,0.5}
    
        coord[0] = -10.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 0.5; // Z coordinate
        geometry_1DL.NodeVec()[0].SetNodeId(node_id);
        geometry_1DL.NodeVec()[0].SetCoord(coord);
        Linear_topology[0] = node_id;
     
    // 2nd node of line element is located at x ={10,PI,0.5}
    node_id++;
    
        coord[0] = 10.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 0.5; // Z coordinate
        geometry_1DL.NodeVec()[1].SetNodeId(node_id);
        geometry_1DL.NodeVec()[1].SetCoord(coord);
        Linear_topology[1] = node_id;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, geometry_1DL);
    }

    geometry_1DL.BuildConnectivity();
    
    
    std::ofstream file("geometry_1DL.txt");
    geometry_1DL.Print(file);
    
}


void TwoDTriElements() {
    
    TPZGeoMesh geometry_2DTri; // Create the objet that will describe the geometry (2D Triangle).
    int n_nodes = 3; // number of nodes
    int geometry_dim = 2; // geometry dimension
    std::string name("geometry 2DTri"); // geometry name
    
    // setting the object
    geometry_2DTri.SetName(name);
    geometry_2DTri.SetDimension(geometry_dim);
    geometry_2DTri.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Triangle_topology(3);
    TPZVec<REAL> coord(3,0.0);
    
    {
    
    // 1st node of triangle element is located at x ={PI,0.0,0.0}
    
        coord[0] = M_PI; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_2DTri.NodeVec()[0].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[0].SetCoord(coord);
        Triangle_topology[0] = node_id;
    
    // 2nd node of triangle element is located at x ={0.0,PI,0.0}
        node_id++;
    
        coord[0] = 0.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_2DTri.NodeVec()[1].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[1].SetCoord(coord);
        Triangle_topology[1] = node_id;
    
    // 3rd node of triangle element is located at x ={0.0,0.0,PI}
        node_id++;

        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_2DTri.NodeVec()[2].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[2].SetCoord(coord);
        Triangle_topology[2] = node_id;
    
        
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (element_id, Triangle_topology, physical_id, geometry_2DTri);
    }
    
    geometry_2DTri.BuildConnectivity();
    
    
    std::ofstream file("geometry_2DTri.txt");
    geometry_2DTri.Print(file);
    
}



void TwoDQuadElements() {
    
    TPZGeoMesh geometry_2DQuad; // Create the objet that will describe the geometry (2D Quadrilateral).
    int n_nodes = 4; // number of nodes
    int geometry_dim = 2; // geometry dimension
    std::string name("geometry 2DQuad"); // geometry name
    
    // setting the object
    geometry_2DQuad.SetName(name);
    geometry_2DQuad.SetDimension(geometry_dim);
    geometry_2DQuad.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Quadrilateral_topology(4);
    TPZVec<REAL> coord(3,0.0);
    
    {
        
        // 1st node of quadrilateral element is located at x ={-1.0,PI,-1.0}
        
        coord[0] = -1.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = -1.0; // Z coordinate
        geometry_2DQuad.NodeVec()[0].SetNodeId(node_id);
        geometry_2DQuad.NodeVec()[0].SetCoord(coord);
        Quadrilateral_topology[0] = node_id;
        
        // 2nd node of quadrilateral element is located at x ={1.0,PI,-1.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = -1.0; // Z coordinate
        geometry_2DQuad.NodeVec()[1].SetNodeId(node_id);
        geometry_2DQuad.NodeVec()[1].SetCoord(coord);
        Quadrilateral_topology[1] = node_id;
        
        // 3rd node of quadrilateral element is located at x ={1.0,PI,1.0}
        node_id++;

        coord[0] = 1.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 1.0; // Z coordinate
        geometry_2DQuad.NodeVec()[2].SetNodeId(node_id);
        geometry_2DQuad.NodeVec()[2].SetCoord(coord);
        Quadrilateral_topology[2] = node_id;
        
        // 4th node of quadrilateral element is located at x ={-1.0,PI,1.0}
        node_id++;

        coord[0] = -1.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 1.0; // Z coordinate
        geometry_2DQuad.NodeVec()[3].SetNodeId(node_id);
        geometry_2DQuad.NodeVec()[3].SetCoord(coord);
        Quadrilateral_topology[3] = node_id;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (element_id, Quadrilateral_topology, physical_id, geometry_2DQuad);
    }
    
    geometry_2DQuad.BuildConnectivity();
    
    
    std::ofstream file("geometry_2DQuad.txt");
    geometry_2DQuad.Print(file);
    
}



void ThreeDTetraElements() {
    
    TPZGeoMesh geometry_3DTetra; // Create the objet that will describe the geometry (3D Tetrahedron).
    int n_nodes = 4; // number of nodes
    int geometry_dim = 3; // geometry dimension
    std::string name("geometry 3DTetra"); // geometry name
    
    // setting the object
    geometry_3DTetra.SetName(name);
    geometry_3DTetra.SetDimension(geometry_dim);
    geometry_3DTetra.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Tetrahedron_topology(4);
    TPZVec<REAL> coord(3,0.0);
    
    {
        
        // 1st node of tetrahedron element is located at x ={1.0,1.0,1.0}
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 1.0; // Z coordinate
        geometry_3DTetra.NodeVec()[0].SetNodeId(node_id);
        geometry_3DTetra.NodeVec()[0].SetCoord(coord);
        Tetrahedron_topology[0] = node_id;
        
        // 2nd node of tetrahedron element is located at x ={1.0,-1.0,-1.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = -1.0; // Z coordinate
        geometry_3DTetra.NodeVec()[1].SetNodeId(node_id);
        geometry_3DTetra.NodeVec()[1].SetCoord(coord);
        Tetrahedron_topology[1] = node_id;
        
        // 3rd node of tetrahedron element is located at x ={-1.0,1.0,-1.0}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = -1.0; // Z coordinate
        geometry_3DTetra.NodeVec()[2].SetNodeId(node_id);
        geometry_3DTetra.NodeVec()[2].SetCoord(coord);
        Tetrahedron_topology[2] = node_id;
        
        // 4th node of tetrahedron element is located at x ={-1.0,-1.0,1.0}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = 1.0; // Z coordinate
        geometry_3DTetra.NodeVec()[3].SetNodeId(node_id);
        geometry_3DTetra.NodeVec()[3].SetCoord(coord);
        Tetrahedron_topology[3] = node_id;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (element_id, Tetrahedron_topology, physical_id, geometry_3DTetra);
    }
    
    geometry_3DTetra.BuildConnectivity();
    
    
    std::ofstream file("geometry_3DTetra.txt");
    geometry_3DTetra.Print(file);
    
}


void ThreeDPyraElements() {
    
    TPZGeoMesh geometry_3DPyra; // Create the objet that will describe the geometry (3D Pyramid).
    int n_nodes = 5; // number of nodes
    int geometry_dim = 3; // geometry dimension
    std::string name("geometry 3DPyra"); // geometry name
    
    // setting the object
    geometry_3DPyra.SetName(name);
    geometry_3DPyra.SetDimension(geometry_dim);
    geometry_3DPyra.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Pyramid_topology(5);
    TPZVec<REAL> coord(3,0.0);
    
    {
        
        // 1st node of pyramid element is located at x ={-1.0,-1.0,0.0}
        
        coord[0] = -1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPyra.NodeVec()[0].SetNodeId(node_id);
        geometry_3DPyra.NodeVec()[0].SetCoord(coord);
        Pyramid_topology[0] = node_id;
        
        // 2nd node of pyramid element is located at x ={1.0,-1.0,0.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPyra.NodeVec()[1].SetNodeId(node_id);
        geometry_3DPyra.NodeVec()[1].SetCoord(coord);
        Pyramid_topology[1] = node_id;
        
        // 3rd node of pyramid element is located at x ={1.0,1.0,0.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPyra.NodeVec()[2].SetNodeId(node_id);
        geometry_3DPyra.NodeVec()[2].SetCoord(coord);
        Pyramid_topology[2] = node_id;
        
        // 4th node of pyramid element is located at x ={-1.0,1.0,0.0}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPyra.NodeVec()[3].SetNodeId(node_id);
        geometry_3DPyra.NodeVec()[3].SetCoord(coord);
        Pyramid_topology[3] = node_id;
        
        // 5th node of pyramid element is located at x ={0.0,0.0,1.0}
        node_id++;
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 1.0; // Z coordinate
        geometry_3DPyra.NodeVec()[4].SetNodeId(node_id);
        geometry_3DPyra.NodeVec()[4].SetCoord(coord);
        Pyramid_topology[4] = node_id;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (element_id, Pyramid_topology, physical_id, geometry_3DPyra);
    }
    
    geometry_3DPyra.BuildConnectivity();
    
    
    std::ofstream file("geometry_3DPyra.txt");
    geometry_3DPyra.Print(file);
    
}

void ThreeDPrisElements() {
    
    TPZGeoMesh geometry_3DPris; // Create the objet that will describe the geometry (3D Prism).
    int n_nodes = 6; // number of nodes
    int geometry_dim = 3; // geometry dimension
    std::string name("geometry 3DPris"); // geometry name
    
    // setting the object
    geometry_3DPris.SetName(name);
    geometry_3DPris.SetDimension(geometry_dim);
    geometry_3DPris.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Prism_topology(6);
    TPZVec<REAL> coord(3,0.0);
    
    {
        
        // 1st node of prism element is located at x ={0.0,0.0,0.0}
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPris.NodeVec()[0].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[0].SetCoord(coord);
        Prism_topology[0] = node_id;
        
        // 2nd node of prism element is located at x ={1.0,0.0,0.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPris.NodeVec()[1].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[1].SetCoord(coord);
        Prism_topology[1] = node_id;
        
        // 3rd node of prism element is located at x ={0.0,1.0,0.0}
        node_id++;
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DPris.NodeVec()[2].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[2].SetCoord(coord);
        Prism_topology[2] = node_id;
        
        // 4th node of prism element is located at x ={0.0,0.0,PI}
        node_id++;
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DPris.NodeVec()[3].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[3].SetCoord(coord);
        Prism_topology[3] = node_id;
        
        // 5th node of prism element is located at x ={1.0,0.0,PI}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DPris.NodeVec()[4].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[4].SetCoord(coord);
        Prism_topology[4] = node_id;
        
        // 6th node of prism element is located at x ={0.0,1.0,PI}
        node_id++;
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DPris.NodeVec()[5].SetNodeId(node_id);
        geometry_3DPris.NodeVec()[5].SetCoord(coord);
        Prism_topology[5] = node_id;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (element_id, Prism_topology, physical_id, geometry_3DPris);
    }
    
    geometry_3DPris.BuildConnectivity();
    
    
    std::ofstream file("geometry_3DPris.txt");
    geometry_3DPris.Print(file);
    
}


void ThreeDHexaElements() {
    
    TPZGeoMesh geometry_3DHexa; // Create the objet that will describe the geometry (3D Hexahedron).
    int n_nodes = 8; // number of nodes
    int geometry_dim = 3; // geometry dimension
    std::string name("geometry 3DHexa"); // geometry name
    
    // setting the object
    geometry_3DHexa.SetName(name);
    geometry_3DHexa.SetDimension(geometry_dim);
    geometry_3DHexa.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> Hexahedron_topology(8);
    TPZVec<REAL> coord(3,0.0);
    
    {
        
        // 1st node of hexahedron element is located at x ={-1.0,-1.0,0.0}
        
        coord[0] = -1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DHexa.NodeVec()[0].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[0].SetCoord(coord);
        Hexahedron_topology[0] = node_id;
        
        // 2nd node of hexahedron element is located at x ={1.0,-1.0,0.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DHexa.NodeVec()[1].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[1].SetCoord(coord);
        Hexahedron_topology[1] = node_id;
        
        // 3rd node of hexahedron element is located at x ={1.0,1.0,0.0}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_3DHexa.NodeVec()[2].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[2].SetCoord(coord);
        Hexahedron_topology[2] = node_id;
        
        // 4th node of hexahedron element is located at x ={-1.0,1.0,0.0}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = 0; // Z coordinate
        geometry_3DHexa.NodeVec()[3].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[3].SetCoord(coord);
        Hexahedron_topology[3] = node_id;
        
        // 5th node of hexahedron element is located at x ={-1.0,-1.0,PI}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DHexa.NodeVec()[4].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[4].SetCoord(coord);
        Hexahedron_topology[4] = node_id;
        
        // 6th node of hexahedron element is located at x ={1.0,-1.0,PI}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = -1.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DHexa.NodeVec()[5].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[5].SetCoord(coord);
        Hexahedron_topology[5] = node_id;
        
        // 7th node of hexahedron element is located at x ={1.0,1.0,PI}
        node_id++;
        
        coord[0] = 1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DHexa.NodeVec()[6].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[6].SetCoord(coord);
        Hexahedron_topology[6] = node_id;
        
        // 8th node of hexahedron element is located at x ={-1.0,1.0,PI}
        node_id++;
        
        coord[0] = -1.0; // x coordinate
        coord[1] = 1.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_3DHexa.NodeVec()[7].SetNodeId(node_id);
        geometry_3DHexa.NodeVec()[7].SetCoord(coord);
        Hexahedron_topology[7] = node_id;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (element_id, Hexahedron_topology, physical_id, geometry_3DHexa);
    }
    
    geometry_3DHexa.BuildConnectivity();
    
    
    std::ofstream file("geometry_3DHexa.txt");
    geometry_3DHexa.Print(file);
    
}

// ************************************** Create 1D linear meshes ***************************************

TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el)
{
    TPZGeoMesh * gmesh_OneDL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("geomesh OneDL"); // geometry name
    
    gmesh_OneDL->SetName(name);
    gmesh_OneDL->SetDimension(geometry_dim);
    
    long num_nodes = num_el + 1; // Number of the nodes
    gmesh_OneDL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
    
    
    const int physical_id = 1; // Define id for material
    const int bc0 = -1; // Define id for left boundary condition
    const int bc1 = -2; // Define id for right boundary condition
    
    for (long i = 0 ; i < num_nodes; i++)
    {
        const REAL valElem = i * size_el;
        TPZVec <REAL> coord(3,0.);
        coord[0] = valElem;
        gmesh_OneDL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
        gmesh_OneDL->NodeVec()[i].SetNodeId(i); // The id identification
    }
    
    // Creating linear element and  zero-dimensional boundary element
    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
    long element_id = 0;

    
    // Elements
    
    for (long iel = 0; iel < num_el; iel++)
    {
        const long inod_l = iel;
        const long inod_r = iel + 1;
        Linear_topology[0] = inod_l;
        Linear_topology[1] = inod_r;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDL);

    }
    
    element_id++;
    
    
    // Left boundary condition
    point_topology[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDL);
    element_id++;
    
    
    // Right boundary condition
    point_topology[0] = num_nodes-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDL);
    
    
    gmesh_OneDL->BuildConnectivity(); // Construct mesh neighbor connectivity
    
    return gmesh_OneDL;
}



// ************************************** Create 1D non linear meshes ***************************************

TPZGeoMesh *CreateOneDNLGMesh(long num_el, REAL size_el)
{
    TPZGeoMesh * gmesh_OneDNL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("geomesh OneDNL"); // geometry name
    
    gmesh_OneDNL->SetName(name);
    gmesh_OneDNL->SetDimension(geometry_dim);
    
    long num_nodes = num_el + 1; // Number of the nodes
    gmesh_OneDNL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
    
    
    const int physical_id = 1; // Define id for material
    const int bc0 = -1; // Define id for left boundary condition
    const int bc1 = -2; // Define id for right boundary condition
    
    for (long i = 0 ; i < num_nodes; i++)
    {
        const REAL valElem = i * size_el;
        TPZVec <REAL> coord(3,0.);
        coord[0] = valElem;
        coord[1] = valElem;
        coord[2] = pow (valElem, 3);

        gmesh_OneDNL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
        gmesh_OneDNL->NodeVec()[i].SetNodeId(i); // The id identification
    }
    
    // Creating linear element and  zero-dimensional boundary element
    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
    long element_id = 0;
    
    
    // Elements
    
    for (long iel = 0; iel < num_el; iel++)
    {
        const long inod_l = iel;
        const long inod_r = iel + 1;
        Linear_topology[0] = inod_l;
        Linear_topology[1] = inod_r;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDNL);
    }
    
    element_id++;

    
    // Left boundary condition
    point_topology[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDNL);
    element_id++;
    
    
    // Right boundary condition
    point_topology[0] = num_nodes-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDNL);
    
    
    gmesh_OneDNL->BuildConnectivity(); // Construct mesh neighbor connectivity
    
    return gmesh_OneDNL;
}



// ************************************** Create 2D simple meshes ***************************************

TPZGeoMesh *CreateTwoDSimpGMesh(int num_divsi, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDSimp = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    long num_Quadnodes = 4; // Number of the nodes

    
    std::string name("geomesh TwoDSimp"); // geometry name
    gmesh_TwoDSimp->SetName(name);
    gmesh_TwoDSimp->SetDimension(geometry_dim);
    
    
    gmesh_TwoDSimp->NodeVec().Resize(num_Quadnodes); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(num_Quadnodes);

    
    
    int physical_id = 1; // Define id for material
    long id = 0;
    long x0 = 0;
    long y0 = 0;
    REAL valx;
    
    
    // Setting node coordantes for quadrilateral element
    for(long i = 0; i < num_Quadnodes/2; i++)
    {
        valx = (i * Lx) + x0;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);       //coord X
        Node[id].SetCoord(1,y0);         //coord Y
        gmesh_TwoDSimp->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(long i = 0; i < num_Quadnodes/2; i++)
    {
        valx = (Lx - i * Lx) + x0;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);      //coord X
        Node[id].SetCoord(1,(Ly + y0));   //coord Y
        gmesh_TwoDSimp->NodeVec()[id] = Node[id];
        id++;
    }
    
    // Index of element
    long elementid = 0;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    
    // Creating quadrilateral element and  one-dimensional boundary element
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec <long> Linear_topology(2);
    
    // Internal elements
    Quadrilateral_topology[0] = 1;
    Quadrilateral_topology[1] = 2;
    Quadrilateral_topology[2] = 3;
    Quadrilateral_topology[3] = 0;
    
    
    // Elements
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDSimp); // create quadrilateral element
    elementid++;
    
    // Boundray elements
    
    Linear_topology[0] = 0;
    Linear_topology[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDSimp); // create boundary element; bottom; bc0
    elementid++;
    
    Linear_topology[0] = 1;
    Linear_topology[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDSimp); // create boundary element; right; bc1
    elementid++;
    
    Linear_topology[0] = 2;
    Linear_topology[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDSimp); // create boundary element; top; bc2
    elementid++;
    
    Linear_topology[0] = 3;
    Linear_topology[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDSimp); // create boundary element; left; bc3
    elementid++;
    
    
    // Build the mesh
    gmesh_TwoDSimp->BuildConnectivity();
    
    // Uniform refinement
    
    
    for( int ref = 0; ref < num_divsi; ref++)
    {
        TPZVec<TPZGeoEl *> children;
        long nel = gmesh_TwoDSimp->NElements();
        for ( long i = 0; i < nel ; i++)
        {
            TPZGeoEl * gel = gmesh_TwoDSimp->ElementVec()[i];
            gel->Divide (children);
        }
    }
    
    return gmesh_TwoDSimp;

}


// ************************************** Create 2D triangle meshes ***************************************

TPZGeoMesh *CreateTwoDTriGMesh(int num_divtri, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDTri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    long numnod_Shape = 9; // Number of the nodes
    
    
    std::string name("geomesh TwoDTri"); // geometry name
    gmesh_TwoDTri->SetName(name);
    gmesh_TwoDTri->SetDimension(geometry_dim);
    
    
    gmesh_TwoDTri->NodeVec().Resize(numnod_Shape); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(numnod_Shape);
    
    
    
    int physical_id = 1; // Define id for material
    long id = 0;
    long x0 = 0;
    long y0 = 0;
    REAL valx;
    REAL valy;
    

    for(long i = 0; i < numnod_Shape/3; i++)
    {
        valx = (i * Lx/2) + x0;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);       //coord X
        
        
        for(long j = 0; j < numnod_Shape/3; j++)
        {
            valy = (j * Ly/2) + y0;

            Node[id].SetCoord(1,valy);   //coord Y
            gmesh_TwoDTri->NodeVec()[id] = Node[id];
            id++;
        }
        
    }
    
    
    // Index of element
    long elementid = 0;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    
    // Creating quadrilateral element and  one-dimensional boundary element
    TPZVec <long> Triangle_topology(3);
    TPZVec <long> Linear_topology(2);
    
    // Internal elements
    Triangle_topology[0] = 0;
    Triangle_topology[1] = 1;
    Triangle_topology[2] = 2;
    
    
    // Elements
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create quadrilateral element
    elementid++;
    
    // Boundray elements
    
    Linear_topology[0] = 0;
    Linear_topology[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom; bc0
    elementid++;
    
    Linear_topology[0] = 1;
    Linear_topology[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; right; bc1
    elementid++;
    
    Linear_topology[0] = 2;
    Linear_topology[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDTri); // create boundary element; top; bc2
    elementid++;
    
    Linear_topology[0] = 3;
    Linear_topology[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; left; bc3
    elementid++;
    
    
    // Build the mesh
    gmesh_TwoDTri->BuildConnectivity();
    
    
    return gmesh_TwoDTri;
    
}




// ************************************** Create 2D quadrilateral meshes ***************************************

TPZGeoMesh *CreateTwoDQuadGMesh(int num_divqu, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDQuad = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    long num_Quadnodes = 4; // Number of the nodes
    
    
    std::string name("geomesh TwoDQuad"); // geometry name
    gmesh_TwoDQuad->SetName(name);
    gmesh_TwoDQuad->SetDimension(geometry_dim);
    
    
    gmesh_TwoDQuad->NodeVec().Resize(num_Quadnodes); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(num_Quadnodes);
    
    
    
    int physical_id = 1; // Define id for material
    long id = 0;
    long x0 = 0;
    long y0 = 0;
    REAL valx;
    
    
    // Setting node coordantes for quadrilateral element
    for(long i = 0; i < num_Quadnodes/2; i++)
    {
        valx = (i * Lx) + x0;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);       //coord X
        Node[id].SetCoord(1,y0);         //coord Y
        gmesh_TwoDQuad->NodeVec()[id] = Node[id];
        id++;
    }
    
    for(long i = 0; i < num_Quadnodes/2; i++)
    {
        valx = (Lx - i * Lx) + x0;
        Node[id].SetNodeId(id);
        Node[id].SetCoord(0,valx);      //coord X
        Node[id].SetCoord(1,(Ly + y0));   //coord Y
        gmesh_TwoDQuad->NodeVec()[id] = Node[id];
        id++;
    }
    
    // Index of element
    long elementid = 0;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    
    // Creating quadrilateral element and  one-dimensional boundary element
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec <long> Linear_topology(2);
    
    // Internal elements
    Quadrilateral_topology[0] = 1;
    Quadrilateral_topology[1] = 2;
    Quadrilateral_topology[2] = 3;
    Quadrilateral_topology[3] = 0;
    
    
    // Elements
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
    elementid++;
    
    // Boundray elements
    
    Linear_topology[0] = 0;
    Linear_topology[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom; bc0
    elementid++;
    
    Linear_topology[0] = 1;
    Linear_topology[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; right; bc1
    elementid++;
    
    Linear_topology[0] = 2;
    Linear_topology[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; top; bc2
    elementid++;
    
    Linear_topology[0] = 3;
    Linear_topology[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; left; bc3
    elementid++;
    
    
    // Build the mesh
    gmesh_TwoDQuad->BuildConnectivity();
    
    // Uniform refinement
    for( int ref = 0; ref < num_divqu; ref++)
    {
        TPZVec<TPZGeoEl *> children;
        long nel = gmesh_TwoDQuad->NElements();
        for ( long i = 0; i < nel-1 ; i++)
        {
            TPZGeoEl * gel = gmesh_TwoDQuad->ElementVec()[i];
            gel->Divide (children);
        }
    }
    
    return gmesh_TwoDQuad;
    
}

