
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


// ******************* Zero and One D element *****************

void ZeroDElements();
void OneDElements();

// ******************* Two D element **************************

void TwoDTriElements();
void TwoDQuadElements();

// ******************* Three D element ************************

void ThreeDTetraElements();
void ThreeDPyraElements();
void ThreeDPrisElements();
void ThreeDHexaElements();


int main() {

    ZeroDElements();
    OneDElements();
    
    TwoDTriElements();
    TwoDQuadElements();
    
    ThreeDTetraElements();
    ThreeDPyraElements();
    ThreeDPrisElements();
    ThreeDHexaElements();


    
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


