
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"


void ZeroDElements();
void OneDElements();
void TwoDTriElements();



int main() {

    ZeroDElements();
    OneDElements();
    TwoDTriElements();
    
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
    // 1st node of element is located at x ={-10,PI,0.5}
    
        coord[0] = -10.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 0.5; // Z coordinate
        geometry_1DL.NodeVec()[0].SetNodeId(node_id);
        geometry_1DL.NodeVec()[0].SetCoord(coord);
        Linear_topology[0] = node_id;
    
    
     
    // 2nd node of element is located at x ={10,PI,0.5}
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
    
    // 1st node of element is located at x ={PI,0.0,0.0}
    
        coord[0] = M_PI; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_2DTri.NodeVec()[0].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[0].SetCoord(coord);
        Triangle_topology[0] = node_id;
   
    
    
    // 2nd node of element is located at x ={0.0,PI,0.0}
    node_id++;
    
        coord[0] = 0.0; // x coordinate
        coord[1] = M_PI; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        geometry_2DTri.NodeVec()[1].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[1].SetCoord(coord);
        Triangle_topology[1] = node_id;
    
    // 3rd node of element is located at x ={0.0,0.0,PI}
   
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = M_PI; // Z coordinate
        geometry_2DTri.NodeVec()[0].SetNodeId(node_id);
        geometry_2DTri.NodeVec()[0].SetCoord(coord);
        Triangle_topology[0] = node_id;
    
        
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (element_id, Triangle_topology, physical_id, geometry_2DTri);
    }
    
    geometry_2DTri.BuildConnectivity();
    
    
    std::ofstream file("geometry_2DTri.txt");
    geometry_2DTri.Print(file);
    
}

