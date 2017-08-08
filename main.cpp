
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"

int main() {

    TPZGeoMesh geometry_1D; // Create the objet that will describe the geometry.
    int n_nodes = 2; // number of nodes
    int geometry_dim = 1; // geometry dimension
    std::string name("geometry 1D"); // geometry name

    // setting the object
    geometry_1D.SetName(name);
    geometry_1D.SetDimension(geometry_dim);
    geometry_1D.NodeVec().Resize(n_nodes);
    
    
    int node_id = 0;
    int element_id = 0;
    int physical_id = 1;
    TPZVec<long> point_topology(1);
    TPZVec<REAL> x(3,0.0);

    // 1st node is located at x ={100,PI,0.5}
    {
        x[0] = 100.0; // x coordinate
        x[1] = M_PI; // Y coordinate
        x[2] = 0.5; // Z coordinate
        geometry_1D.NodeVec()[0].SetNodeId(node_id);
        geometry_1D.NodeVec()[0].SetCoord(x);
        point_topology[0] = node_id;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1D);
    }
    
    // 2nd node is located at x ={-100,PI,0.5}
    node_id++;
    element_id++;
    {
        x[0] = -100.0; // x coordinate
        x[1] = M_PI; // Y coordinate
        x[2] = 0.5; // Z coordinate
        geometry_1D.NodeVec()[1].SetNodeId(node_id);
        geometry_1D.NodeVec()[1].SetCoord(x);
        point_topology[0] = node_id;
        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1D);
    }
    
    
    geometry_1D.BuildConnectivity();
    
    
    std::ofstream file("geometry_1D.txt");
    geometry_1D.Print(file);

    return 0;
}
