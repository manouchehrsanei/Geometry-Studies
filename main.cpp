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


#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticprism.h"

#include "tpzarc3d.h"
#include "tpzellipse3d.h"
#include "tpzgeoblend.h"

#include "pzvec.h"
#include "pzquad.h"
#include "tpzintrulep3d.h"


#include "Topology.h"



// ************************************* (Geometry linear description: Zero and One D element) ******************

void ZeroDElements();
void OneDElements();

// ************************************* (Geometry linear description: Two D element) ***************************

void TwoDTriElements();
void TwoDQuadElements();

// ************************************* (Geometry linear description: Three D element) *************************

void ThreeDTetraElements();
void ThreeDPyraElements();
void ThreeDPrisElements();
void ThreeDHexaElements();

// ************************************** (Create linear meshes: 1D) ********************************************

TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el);

TPZGeoMesh *CreateOneDNLGMesh(long num_el, REAL size_el);



// ************************************** (Create linear meshes: 2D) ********************************************

TPZGeoMesh *CreateTwoDSimpGMesh(int num_div, REAL Lx, REAL Ly);

TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly);

TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly);


// ************************************* (Create linear meshes: 3D) ********************************************

TPZGeoMesh *CreateThreeDHexPriGMesh(long nnodesthr, REAL Lx, REAL Ly, REAL Lz);

TPZGeoMesh *CreateThreeDHexPrytetGMesh(long nnodesthrhpt, REAL Lx, REAL Ly, REAL Lz);


// ********************************** (Geometry nonlinear description: One D element) **************************

void NonOneDElements();

// ************************************* (Geometry linear description: Two D element) **************************

void NonTwoDTriElements();
void NonTwoDQuadElements();

// ************************************* (Geometry linear description: Three D element) *************************

void NonThreeDTetraElements();
void NonThreeDPyraElements();
void NonThreeDPrisElements();
void NonThreeDHexaElements();


// ************************************* (Create nonlinear meshes: 3D) ********************************************

TPZGeoMesh *CreateNonThreeDHexaGMesh(long nnodesthrehex, REAL Lx, REAL Ly, REAL Lz);

TPZGeoMesh *CreateThreeDarcGMesh(long nnodesarc, REAL Rad);

// --------------------------------------

void ThreeEllipsMeshGenerate();


// ************************************* (Finding nodes and elements) ********************************************

TPZGeoMesh *CreateOneDFindNoElGMesh(long num_elem, REAL size_elem);



// ************************************* (Integral) ********************************************









// ******************************************** (main of program) ***********************************************


int main() {
    
    CheckTopology();


    
    return 0;
    
}

    // ********************************* (linear elements) ******************************************************
//
//    ZeroDElements();
//    OneDElements();
//    
//    TwoDTriElements();
//    TwoDQuadElements();
//    
//    ThreeDTetraElements();
//    ThreeDPyraElements();
//    ThreeDPrisElements();
//    ThreeDHexaElements();
//    
//    
//    // ******************* (Create linear meshes: 1D) ***********************************************************
//    
//    REAL domain = 1.0;
//    long num_el = 10;
//    REAL size_el = domain/num_el;
//    
//    TPZGeoMesh *gmesh_OneDL = CreateOneDLGMesh(num_el, size_el); // function to create the 1D geometric mesh
//    
//
//
//    // ---------------------------------------------------------------------------------------
//    
//    
//    TPZGeoMesh *gmesh_OneDNL = CreateOneDNLGMesh(num_el, size_el); // function to create the 1D geometric mesh
//    
//
//    
//    // ********************************* (Create linear meshes: 2D) *******************************************
//
//    long num_divsi = 2; // number of divition
//    REAL Lx = 1.0; // length of domain in x direction
//    REAL Ly = 1.0; // length of domain in y direction
//    
//    TPZGeoMesh *gmesh_TwoDSimp = CreateTwoDSimpGMesh(num_divsi, Lx, Ly); // function to create the 2D geometric mesh
//
//    
//    
//    // ---------------------------------------------------------------------------------------
//    
//    long nnodes = 9; // Number of the nodes
//    
//    TPZGeoMesh *gmesh_TwoDTri = CreateTwoDTriGMesh(nnodes, Lx, Ly); // function to create the 2D geometric mesh
//    
//
//    
//    
//    // ---------------------------------------------------------------------------------------
//    
//    long nnodesqu = 10; // number of divition
//    
//    TPZGeoMesh *gmesh_TwoDQuad = CreateTwoDQuadGMesh(nnodesqu, Lx, Ly); // function to create the 2D geometric mesh
//    
//
//    
//    // ****************************** (Create linear meshes: 3D) ********************************************
//    
//    REAL Lz = 1.0;
//    long nnodesthr = 25; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDHexPri = CreateThreeDHexPriGMesh(nnodesthr, Lx, Ly, Lz); // function to create the 3D geometric mesh
//
//    
//    // ----------------------------------------------------------------------------------------
//    
//    long nnodesthrhpt = 12; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDHexPrytet = CreateThreeDHexPrytetGMesh(nnodesthrhpt, Lx, Ly, Lz); // function to create the 3D geometric mesh
//    
//
//    
//    
//    // *************************** (nonlinear elements) *****************************************************
//    
//    NonOneDElements();
//    
//    NonTwoDTriElements();
//    NonTwoDQuadElements();
//    
//    NonThreeDTetraElements();
//    NonThreeDPyraElements();
//    NonThreeDPrisElements();
//    NonThreeDHexaElements();
//    
//    
//    // ----------------------------------------------------------------------------------------
//    
//    long nnodesthrehex = 20; // number of divition
//    
//    TPZGeoMesh *gmesh_NonThreeDHexa = CreateNonThreeDHexaGMesh(nnodesthrehex, Lx, Ly, Lz); // function to create the 3D geometric mesh
//
//    
//    
//    // ----------------------------------------------------------------------------------------
//
//    REAL Rad = 1.0;
//    long nnodesarc = 9; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDarc = CreateThreeDarcGMesh(nnodesarc, Rad); // function to create the 3D geometric mesh
//    
//
//
//
////     ---------------------------------------  check shavad
////    
////    
////    ThreeEllipsMeshGenerate();
//
//    
//    // ************************************* (Finding nodes and elements) ********************************************
//    
//    REAL OneDdomain = 1.0;
//    long num_elem = 20;
//    REAL size_elem = OneDdomain/num_elem;
//    
//    TPZGeoMesh *gmesh_OneDFindNoEl = CreateOneDFindNoElGMesh(num_elem, size_elem); // function to create the 1D geometric mesh
//    
//    
//    // ------------------------------ finding the boundary element anc change the id
////    long target_id = 100;
////    
////    long nel = gmesh_OneDFindNoEl->NElements();
////    
////    for (long iel = 0; iel < nel; iel++)
////    {
////        TPZGeoEl * element = gmesh_OneDFindNoEl->Element(iel);
////        
////        int dim = element->Dimension();
////        if (dim == 0)
////        { //Filtering boundary elements
////            std::cout << "before = " <<std::endl;
////            element->Print();
////    
////            element->SetMaterialId(target_id);
////            std::cout << "after = " <<std::endl;
////            element->Print();
////        }
////    }
//    
//    // --------------------------- finding the espeical element and change the coordiantes
//    
////    long target_index = 10;
////    REAL z_modified = 1.0;
////    
////    bool one_way = false;
////    
////    TPZVec<REAL> x(3,0.0);
////    
////    TPZGeoEl * element = gmesh_OneDFindNoEl->Element(target_index);
////    int n_nodes = element->NNodes();
////    for (int i = 0; i < n_nodes; i++) {
////        TPZGeoNode & inode = element->Node(i);
////        
////        if (one_way) {
////            inode.GetCoordinates(x);
////            x[2] = z_modified;
////            inode.SetCoord(x);
////        }
////        else{
////            inode.SetCoord(2, z_modified);
////        }
////        
////    }
////    
////    std::ofstream outgmeshOneDFindNoEl("geomesh_OneDFindNoElModified.txt");
////    gmesh_OneDFindNoEl->Print(outgmeshOneDFindNoEl);
////    
////    std::ofstream vtkgmeshOneDFindNoEl("geomesh_OneDFindNoElModified.vtk");
////    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDFindNoEl, vtkgmeshOneDFindNoEl);
//    
//    // --------------------------- finding the jacobian of espeical element
//    
//    long target_index = 10;
//    REAL z_modified = 1.0;
//    TPZVec<REAL> x(3,0.0);
//
//    TPZVec<REAL> par_coordinate(3,0.0);
//    TPZFMatrix<REAL> jacobian(3,3);
//    TPZFMatrix<REAL> Axes(3,3);
//    REAL detJacobian;
//    REAL detJacobianSub;
//    REAL detJacobianPlu;
//
//
//    TPZFMatrix<REAL> InvJac(3,3);
//    
//    // -------------------- the element before the target index ------------
//    
//
//    TPZGeoEl * elementTarget = gmesh_OneDFindNoEl->Element(target_index);
//    int numnodes = elementTarget->NNodes();
//    for (int i = 0; i < numnodes; i++)
//    {
//        TPZGeoNode & inode = elementTarget->Node(i);
//        
//        inode.GetCoordinates(x);
//        x[2] = z_modified;
//        inode.SetCoord(x);
//
//    }
//    
//    TPZManVector<REAL,3> x_1(3,0.0);
//    elementTarget->X(par_coordinate, x_1);
//    elementTarget->Jacobian(par_coordinate, jacobian, Axes, detJacobian, InvJac);
//    
//    TPZManVector<REAL,3> x_2(3,0.0);
//    TPZGeoEl * elementTarSub = gmesh_OneDFindNoEl->Element(target_index-1);
//    elementTarSub->X(par_coordinate, x_2);
//    elementTarSub->Jacobian(par_coordinate, jacobian, Axes, detJacobianSub, InvJac);
//    
//    
//    TPZManVector<REAL,3> x_3(3,0.0);
//    TPZGeoEl * elementTarPlu = gmesh_OneDFindNoEl->Element(target_index+1);
//    elementTarPlu->X(par_coordinate, x_3);
//    elementTarPlu->Jacobian(par_coordinate, jacobian, Axes, detJacobianPlu, InvJac);
//    
//    
//    
//    
//    
//    
//    return 0;
//}
//
//
//
//
//// ++++++++++++++++++++++++++++++++++++++++++++ functions +++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
//
//// ********************************* (linear elements) ******************************************************
//
//
//
//void ZeroDElements(){
//    
//    TPZGeoMesh geometry_1DP; // Create the objet that will describe the geometry (Point).
//    int n_nodes = 3; // number of nodes
//    int geometry_dim = 1; // geometry dimension
//    std::string name("geometry 1DP"); // geometry name
//    
//    // setting the object
//    geometry_1DP.SetName(name);
//    geometry_1DP.SetDimension(geometry_dim);
//    geometry_1DP.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> point_topology(1);
//    TPZVec<REAL> x(3,0.0);
//    
//    // 1st node is located at x ={-100,PI,0.5}
//    {
//        x[0] = -100.0; // x coordinate
//        x[1] = M_PI; // Y coordinate
//        x[2] = 0.5; // Z coordinate
//        geometry_1DP.NodeVec()[0].SetNodeId(node_id);
//        geometry_1DP.NodeVec()[0].SetCoord(x);
//        point_topology[0] = node_id;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
//    }
//    
//    // 2nd node is located at x ={0.0,PI,0.5}
//    node_id++;
//    element_id++;
//    {
//        x[0] = 0.0; // x coordinate
//        x[1] = M_PI; // Y coordinate
//        x[2] = 0.5; // Z coordinate
//        geometry_1DP.NodeVec()[1].SetNodeId(node_id);
//        geometry_1DP.NodeVec()[1].SetCoord(x);
//        point_topology[0] = node_id;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
//    }
//    
//    // 3nd node is located at x ={100,PI,0.5}
//    node_id++;
//    element_id++;
//    {
//        x[0] = 100.0; // x coordinate
//        x[1] = M_PI; // Y coordinate
//        x[2] = 0.5; // Z coordinate
//        geometry_1DP.NodeVec()[2].SetNodeId(node_id);
//        geometry_1DP.NodeVec()[2].SetCoord(x);
//        point_topology[0] = node_id;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (element_id, point_topology, physical_id, geometry_1DP);
//    }
//    
//    geometry_1DP.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_1DP.txt");
//    geometry_1DP.Print(file);
//    
//}
//// ------------------------------------------------------------------------------------------
//
//void OneDElements(){
//    
//    TPZGeoMesh geometry_1DL; // Create the objet that will describe the geometry (Line).
//    int n_nodes = 2; // number of nodes
//    int geometry_dim = 1; // geometry dimension
//    std::string name("geometry 1DL"); // geometry name
//    
//    // setting the object
//    geometry_1DL.SetName(name);
//    geometry_1DL.SetDimension(geometry_dim);
//    geometry_1DL.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Linear_topology(2);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//    // 1st node of line element is located at x ={-10,PI,0.5}
//    
//        coord[0] = -10.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_1DL.NodeVec()[0].SetNodeId(node_id);
//        geometry_1DL.NodeVec()[0].SetCoord(coord);
//        Linear_topology[0] = node_id;
//     
//    // 2nd node of line element is located at x ={10,PI,0.5}
//    node_id++;
//    
//        coord[0] = 10.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_1DL.NodeVec()[1].SetNodeId(node_id);
//        geometry_1DL.NodeVec()[1].SetCoord(coord);
//        Linear_topology[1] = node_id;
//    
////    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, geometry_1DL);
//    }
//
//    geometry_1DL.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_1DL.txt");
//    geometry_1DL.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//void TwoDTriElements() {
//    
//    TPZGeoMesh geometry_2DTri; // Create the objet that will describe the geometry (2D Triangle).
//    int n_nodes = 3; // number of nodes
//    int geometry_dim = 2; // geometry dimension
//    std::string name("geometry 2DTri"); // geometry name
//    
//    // setting the object
//    geometry_2DTri.SetName(name);
//    geometry_2DTri.SetDimension(geometry_dim);
//    geometry_2DTri.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Triangle_topology(3);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//    
//    // 1st node of triangle element is located at x ={PI,0.0,0.0}
//    
//        coord[0] = M_PI; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_2DTri.NodeVec()[0].SetNodeId(node_id);
//        geometry_2DTri.NodeVec()[0].SetCoord(coord);
//        Triangle_topology[0] = node_id;
//    
//    // 2nd node of triangle element is located at x ={0.0,PI,0.0}
//        node_id++;
//    
//        coord[0] = 0.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_2DTri.NodeVec()[1].SetNodeId(node_id);
//        geometry_2DTri.NodeVec()[1].SetCoord(coord);
//        Triangle_topology[1] = node_id;
//    
//    // 3rd node of triangle element is located at x ={0.0,0.0,PI}
//        node_id++;
//
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_2DTri.NodeVec()[2].SetNodeId(node_id);
//        geometry_2DTri.NodeVec()[2].SetCoord(coord);
//        Triangle_topology[2] = node_id;
//    
//        
////    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (element_id, Triangle_topology, physical_id, geometry_2DTri);
//    }
//    
//    geometry_2DTri.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_2DTri.txt");
//    geometry_2DTri.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void TwoDQuadElements() {
//    
//    TPZGeoMesh geometry_2DQuad; // Create the objet that will describe the geometry (2D Quadrilateral).
//    int n_nodes = 4; // number of nodes
//    int geometry_dim = 2; // geometry dimension
//    std::string name("geometry 2DQuad"); // geometry name
//    
//    // setting the object
//    geometry_2DQuad.SetName(name);
//    geometry_2DQuad.SetDimension(geometry_dim);
//    geometry_2DQuad.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Quadrilateral_topology(4);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of quadrilateral element is located at x ={-1.0,PI,-1.0}
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_2DQuad.NodeVec()[0].SetNodeId(node_id);
//        geometry_2DQuad.NodeVec()[0].SetCoord(coord);
//        Quadrilateral_topology[0] = node_id;
//        
//        // 2nd node of quadrilateral element is located at x ={1.0,PI,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_2DQuad.NodeVec()[1].SetNodeId(node_id);
//        geometry_2DQuad.NodeVec()[1].SetCoord(coord);
//        Quadrilateral_topology[1] = node_id;
//        
//        // 3rd node of quadrilateral element is located at x ={1.0,PI,1.0}
//        node_id++;
//
//        coord[0] = 1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_2DQuad.NodeVec()[2].SetNodeId(node_id);
//        geometry_2DQuad.NodeVec()[2].SetCoord(coord);
//        Quadrilateral_topology[2] = node_id;
//        
//        // 4th node of quadrilateral element is located at x ={-1.0,PI,1.0}
//        node_id++;
//
//        coord[0] = -1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_2DQuad.NodeVec()[3].SetNodeId(node_id);
//        geometry_2DQuad.NodeVec()[3].SetCoord(coord);
//        Quadrilateral_topology[3] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (element_id, Quadrilateral_topology, physical_id, geometry_2DQuad);
//    }
//    
//    geometry_2DQuad.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_2DQuad.txt");
//    geometry_2DQuad.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void ThreeDTetraElements() {
//    
//    TPZGeoMesh geometry_3DTetra; // Create the objet that will describe the geometry (3D Tetrahedron).
//    int n_nodes = 4; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry 3DTetra"); // geometry name
//    
//    // setting the object
//    geometry_3DTetra.SetName(name);
//    geometry_3DTetra.SetDimension(geometry_dim);
//    geometry_3DTetra.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Tetrahedron_topology(4);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of tetrahedron element is located at x ={1.0,1.0,1.0}
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_3DTetra.NodeVec()[0].SetNodeId(node_id);
//        geometry_3DTetra.NodeVec()[0].SetCoord(coord);
//        Tetrahedron_topology[0] = node_id;
//        
//        // 2nd node of tetrahedron element is located at x ={1.0,-1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_3DTetra.NodeVec()[1].SetNodeId(node_id);
//        geometry_3DTetra.NodeVec()[1].SetCoord(coord);
//        Tetrahedron_topology[1] = node_id;
//        
//        // 3rd node of tetrahedron element is located at x ={-1.0,1.0,-1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_3DTetra.NodeVec()[2].SetNodeId(node_id);
//        geometry_3DTetra.NodeVec()[2].SetCoord(coord);
//        Tetrahedron_topology[2] = node_id;
//        
//        // 4th node of tetrahedron element is located at x ={-1.0,-1.0,1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_3DTetra.NodeVec()[3].SetNodeId(node_id);
//        geometry_3DTetra.NodeVec()[3].SetCoord(coord);
//        Tetrahedron_topology[3] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (element_id, Tetrahedron_topology, physical_id, geometry_3DTetra);
//    }
//    
//    geometry_3DTetra.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_3DTetra.txt");
//    geometry_3DTetra.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void ThreeDPyraElements() {
//    
//    TPZGeoMesh geometry_3DPyra; // Create the objet that will describe the geometry (3D Pyramid).
//    int n_nodes = 5; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry 3DPyra"); // geometry name
//    
//    // setting the object
//    geometry_3DPyra.SetName(name);
//    geometry_3DPyra.SetDimension(geometry_dim);
//    geometry_3DPyra.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Pyramid_topology(5);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of pyramid element is located at x ={-1.0,-1.0,0.0}
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPyra.NodeVec()[0].SetNodeId(node_id);
//        geometry_3DPyra.NodeVec()[0].SetCoord(coord);
//        Pyramid_topology[0] = node_id;
//        
//        // 2nd node of pyramid element is located at x ={1.0,-1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPyra.NodeVec()[1].SetNodeId(node_id);
//        geometry_3DPyra.NodeVec()[1].SetCoord(coord);
//        Pyramid_topology[1] = node_id;
//        
//        // 3rd node of pyramid element is located at x ={1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPyra.NodeVec()[2].SetNodeId(node_id);
//        geometry_3DPyra.NodeVec()[2].SetCoord(coord);
//        Pyramid_topology[2] = node_id;
//        
//        // 4th node of pyramid element is located at x ={-1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPyra.NodeVec()[3].SetNodeId(node_id);
//        geometry_3DPyra.NodeVec()[3].SetCoord(coord);
//        Pyramid_topology[3] = node_id;
//        
//        // 5th node of pyramid element is located at x ={0.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_3DPyra.NodeVec()[4].SetNodeId(node_id);
//        geometry_3DPyra.NodeVec()[4].SetCoord(coord);
//        Pyramid_topology[4] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (element_id, Pyramid_topology, physical_id, geometry_3DPyra);
//    }
//    
//    geometry_3DPyra.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_3DPyra.txt");
//    geometry_3DPyra.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void ThreeDPrisElements() {
//    
//    TPZGeoMesh geometry_3DPris; // Create the objet that will describe the geometry (3D Prism).
//    int n_nodes = 6; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry 3DPris"); // geometry name
//    
//    // setting the object
//    geometry_3DPris.SetName(name);
//    geometry_3DPris.SetDimension(geometry_dim);
//    geometry_3DPris.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Prism_topology(6);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of prism element is located at x ={0.0,0.0,0.0}
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPris.NodeVec()[0].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[0].SetCoord(coord);
//        Prism_topology[0] = node_id;
//        
//        // 2nd node of prism element is located at x ={1.0,0.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPris.NodeVec()[1].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[1].SetCoord(coord);
//        Prism_topology[1] = node_id;
//        
//        // 3rd node of prism element is located at x ={0.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DPris.NodeVec()[2].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[2].SetCoord(coord);
//        Prism_topology[2] = node_id;
//        
//        // 4th node of prism element is located at x ={0.0,0.0,PI}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DPris.NodeVec()[3].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[3].SetCoord(coord);
//        Prism_topology[3] = node_id;
//        
//        // 5th node of prism element is located at x ={1.0,0.0,PI}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DPris.NodeVec()[4].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[4].SetCoord(coord);
//        Prism_topology[4] = node_id;
//        
//        // 6th node of prism element is located at x ={0.0,1.0,PI}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DPris.NodeVec()[5].SetNodeId(node_id);
//        geometry_3DPris.NodeVec()[5].SetCoord(coord);
//        Prism_topology[5] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (element_id, Prism_topology, physical_id, geometry_3DPris);
//    }
//    
//    geometry_3DPris.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_3DPris.txt");
//    geometry_3DPris.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void ThreeDHexaElements() {
//    
//    TPZGeoMesh geometry_3DHexa; // Create the objet that will describe the geometry (3D Hexahedron).
//    int n_nodes = 8; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry 3DHexa"); // geometry name
//    
//    // setting the object
//    geometry_3DHexa.SetName(name);
//    geometry_3DHexa.SetDimension(geometry_dim);
//    geometry_3DHexa.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> Hexahedron_topology(8);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of hexahedron element is located at x ={-1.0,-1.0,0.0}
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DHexa.NodeVec()[0].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[0].SetCoord(coord);
//        Hexahedron_topology[0] = node_id;
//        
//        // 2nd node of hexahedron element is located at x ={1.0,-1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DHexa.NodeVec()[1].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[1].SetCoord(coord);
//        Hexahedron_topology[1] = node_id;
//        
//        // 3rd node of hexahedron element is located at x ={1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_3DHexa.NodeVec()[2].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[2].SetCoord(coord);
//        Hexahedron_topology[2] = node_id;
//        
//        // 4th node of hexahedron element is located at x ={-1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0; // Z coordinate
//        geometry_3DHexa.NodeVec()[3].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[3].SetCoord(coord);
//        Hexahedron_topology[3] = node_id;
//        
//        // 5th node of hexahedron element is located at x ={-1.0,-1.0,PI}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DHexa.NodeVec()[4].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[4].SetCoord(coord);
//        Hexahedron_topology[4] = node_id;
//        
//        // 6th node of hexahedron element is located at x ={1.0,-1.0,PI}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DHexa.NodeVec()[5].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[5].SetCoord(coord);
//        Hexahedron_topology[5] = node_id;
//        
//        // 7th node of hexahedron element is located at x ={1.0,1.0,PI}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DHexa.NodeVec()[6].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[6].SetCoord(coord);
//        Hexahedron_topology[6] = node_id;
//        
//        // 8th node of hexahedron element is located at x ={-1.0,1.0,PI}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_3DHexa.NodeVec()[7].SetNodeId(node_id);
//        geometry_3DHexa.NodeVec()[7].SetCoord(coord);
//        Hexahedron_topology[7] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (element_id, Hexahedron_topology, physical_id, geometry_3DHexa);
//    }
//    
//    geometry_3DHexa.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_3DHexa.txt");
//    geometry_3DHexa.Print(file);
//    
//}
//
//// ************************************** Create 1D linear meshes ***************************************
//
//TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el)
//{
//    TPZGeoMesh * gmesh_OneDL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 1; // geometry dimension
//    std::string name("geomesh OneDL"); // geometry name
//    
//    gmesh_OneDL->SetName(name);
//    gmesh_OneDL->SetDimension(geometry_dim);
//    
//    long num_nodes = num_el + 1; // Number of the nodes
//    gmesh_OneDL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
//    
//    
//    const int physical_id = 1; // Define id for material
//    const int bc0 = -1; // Define id for left boundary condition
//    const int bc1 = -2; // Define id for right boundary condition
//    
//    for (long i = 0 ; i < num_nodes; i++)
//    {
//        const REAL valElem = i * size_el;
//        TPZVec <REAL> coord(3,0.);
//        coord[0] = valElem;
//        gmesh_OneDL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
//        gmesh_OneDL->NodeVec()[i].SetNodeId(i); // The id identification
//    }
//    
//    // Creating linear element and  zero-dimensional boundary element
//    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
//    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
//    long element_id = 0;
//
//    
//    // Elements
//    
//    for (long iel = 0; iel < num_el; iel++)
//    {
//        const long inod_l = iel;
//        const long inod_r = iel + 1;
//        Linear_topology[0] = inod_l;
//        Linear_topology[1] = inod_r;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDL);
//
//    }
//    
//    element_id++;
//    
//    
//    // Left boundary condition
//    point_topology[0] = 0;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDL);
//    element_id++;
//    
//    
//    // Right boundary condition
//    point_topology[0] = num_nodes-1;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDL);
//    
//    
//    gmesh_OneDL->BuildConnectivity(); // Construct mesh neighbor connectivity
//    
//    
//    std::ofstream outgmeshOneDL("geomesh_OneDL.txt");
//    gmesh_OneDL->Print(outgmeshOneDL);
//    
//    std::ofstream vtkgmeshOneDL("geomesh_OneDL.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDL, vtkgmeshOneDL);
//    
//    
//    return gmesh_OneDL;
//}
//
//
//
//// ************************************** Create 1D non linear meshes ***************************************
//
//TPZGeoMesh *CreateOneDNLGMesh(long num_el, REAL size_el)
//{
//    TPZGeoMesh * gmesh_OneDNL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 1; // geometry dimension
//    std::string name("geomesh OneDNL"); // geometry name
//    
//    gmesh_OneDNL->SetName(name);
//    gmesh_OneDNL->SetDimension(geometry_dim);
//    
//    long num_nodes = num_el + 1; // Number of the nodes
//    gmesh_OneDNL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
//    
//    
//    const int physical_id = 1; // Define id for material
//    const int bc0 = -1; // Define id for left boundary condition
//    const int bc1 = -2; // Define id for right boundary condition
//    
//    for (long i = 0 ; i < num_nodes; i++)
//    {
//        const REAL valElem = i * size_el;
//        TPZVec <REAL> coord(3,0.);
//        coord[0] = valElem;
//        coord[1] = valElem;
//        coord[2] = pow (valElem, 3);
//
//        gmesh_OneDNL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
//        gmesh_OneDNL->NodeVec()[i].SetNodeId(i); // The id identification
//    }
//    
//    // Creating linear element and  zero-dimensional boundary element
//    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
//    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
//    long element_id = 0;
//    
//    
//    // Elements
//    
//    for (long iel = 0; iel < num_el; iel++)
//    {
//        const long inod_l = iel;
//        const long inod_r = iel + 1;
//        Linear_topology[0] = inod_l;
//        Linear_topology[1] = inod_r;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDNL);
//    }
//    
//    element_id++;
//
//    
//    // Left boundary condition
//    point_topology[0] = 0;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDNL);
//    element_id++;
//    
//    
//    // Right boundary condition
//    point_topology[0] = num_nodes-1;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDNL);
//    
//    
//    gmesh_OneDNL->BuildConnectivity(); // Construct mesh neighbor connectivity
//    
//    std::ofstream outgmeshOneDNL("geomesh_OneDNL.txt");
//    gmesh_OneDNL->Print(outgmeshOneDNL);
//    
//    std::ofstream vtkgmeshOneDNL("geomesh_OneDNL.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDNL, vtkgmeshOneDNL);
//    
//    return gmesh_OneDNL;
//}
//
//
//
//// ************************************** Create 2D simple meshes ***************************************
//
//TPZGeoMesh *CreateTwoDSimpGMesh(int num_divsi, REAL Lx, REAL Ly)
//{
//    TPZGeoMesh * gmesh_TwoDSimp = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 2; // geometry dimension
//    long num_Quadnodes = 4; // Number of the nodes
//
//    
//    std::string name("geomesh TwoDSimp"); // geometry name
//    gmesh_TwoDSimp->SetName(name);
//    gmesh_TwoDSimp->SetDimension(geometry_dim);
//    
//    
//    gmesh_TwoDSimp->NodeVec().Resize(num_Quadnodes); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(num_Quadnodes);
//
//    
//    
//    int physical_id = 1; // Define id for material
//    long id = 0;
//    long x0 = 0;
//    long y0 = 0;
//    REAL valx;
//    
//    
//    // Setting node coordantes for quadrilateral element
//    for(long i = 0; i < num_Quadnodes/2; i++)
//    {
//        valx = (i * Lx) + x0;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0,valx);       //coord X
//        Node[id].SetCoord(1,y0);         //coord Y
//        gmesh_TwoDSimp->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    for(long i = 0; i < num_Quadnodes/2; i++)
//    {
//        valx = (Lx - i * Lx) + x0;
//        Node[id].SetNodeId(id);
//        Node[id].SetCoord(0,valx);      //coord X
//        Node[id].SetCoord(1,(Ly + y0));   //coord Y
//        gmesh_TwoDSimp->NodeVec()[id] = Node[id];
//        id++;
//    }
//    
//    // Index of element
//    long elementid = 0;
//    
//    // Index of boundary element
//    const int bc_bottom = -1; // define id for a material (border bottom)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_top = -3; // define id for a material (border top)
//    const int bc_left = -4; // define id for a material (border left)
//    
//    
//    // Creating quadrilateral element and  one-dimensional boundary element
//    TPZVec <long> Quadrilateral_topology(4);
//    TPZVec <long> Linear_topology(2);
//    
//    // Internal elements
//    Quadrilateral_topology[0] = 1;
//    Quadrilateral_topology[1] = 2;
//    Quadrilateral_topology[2] = 3;
//    Quadrilateral_topology[3] = 0;
//    
//    
//    // Elements
//    
////    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDSimp); // create quadrilateral element
//    elementid++;
//    
//    // Boundray elements
//    
//    Linear_topology[0] = 0;
//    Linear_topology[1] = 1;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDSimp); // create boundary element; bottom; bc0
//    elementid++;
//    
//    Linear_topology[0] = 1;
//    Linear_topology[1] = 2;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDSimp); // create boundary element; right; bc1
//    elementid++;
//    
//    Linear_topology[0] = 2;
//    Linear_topology[1] = 3;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDSimp); // create boundary element; top; bc2
//    elementid++;
//    
//    Linear_topology[0] = 3;
//    Linear_topology[1] = 0;
////    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDSimp); // create boundary element; left; bc3
//    elementid++;
//    
//    
//    // Build the mesh
//    gmesh_TwoDSimp->BuildConnectivity();
//    
//    // Uniform refinement
//    
//    
//    for( int ref = 0; ref < num_divsi; ref++)
//    {
//        TPZVec<TPZGeoEl *> children;
//        long nel = gmesh_TwoDSimp->NElements();
//        for ( long i = 0; i < nel ; i++)
//        {
//            TPZGeoEl * gel = gmesh_TwoDSimp->ElementVec()[i];
//            gel->Divide (children);
//        }
//    }
//    
//    
//    std::ofstream outgmeshTwoDSimp("geomesh_TwoDSimp.txt");
//    gmesh_TwoDSimp->Print(outgmeshTwoDSimp);
//    
//    std::ofstream vtkgmeshTwoDSimp("geomesh_TwoDSimp.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDSimp, vtkgmeshTwoDSimp);
//    
//    
//    return gmesh_TwoDSimp;
//
//}
//
//
//// ************************************** Create 2D triangle meshes ***************************************
//
//TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly)
//{
//    TPZGeoMesh * gmesh_TwoDTri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 2; // geometry dimension
//    
//    std::string name("geomesh TwoDTri"); // geometry name
//    gmesh_TwoDTri->SetName(name);
//    gmesh_TwoDTri->SetDimension(geometry_dim);
//    
//    
//    gmesh_TwoDTri->NodeVec().Resize(nnodes); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodes);
//    
//    
//    TPZVec<long> Triangle_topology(3);
//    TPZVec <long> Linear_topology(2);
//
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//    long elementid = 0;
//    int physical_id = 1;
//
//    // Index of boundary element
//    const int bc_bottom = -1; // define id for a material (border bottom)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_top = -3; // define id for a material (border top)
//    const int bc_left = -4; // define id for a material (border left)
//    
//    {
//        
//        // 0th element
//        {
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = Ly/2; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[0].SetNodeId(3);
//            gmesh_TwoDTri->NodeVec()[0].SetCoord(coord);
//            Triangle_topology[0] = 0; // index
//            
//            coord[0] = Lx/2; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[1].SetNodeId(7);
//            gmesh_TwoDTri->NodeVec()[1].SetCoord(coord);
//            Triangle_topology[1] = 1;
//            
//            coord[0] = Lx/2; // x coordinate
//            coord[1] = Ly/2; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[2].SetNodeId(4);
//            gmesh_TwoDTri->NodeVec()[2].SetCoord(coord);
//            Triangle_topology[2] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        
//        // 1st element
//        
//        {
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_TwoDTri->NodeVec()[3].SetNodeId(0);
//        gmesh_TwoDTri->NodeVec()[3].SetCoord(coord);
//        Triangle_topology[0] = 3;
//            
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_TwoDTri->NodeVec()[4].SetNodeId(1);
//        gmesh_TwoDTri->NodeVec()[4].SetCoord(coord);
//        Triangle_topology[1] = 4;
//            
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_TwoDTri->NodeVec()[5].SetNodeId(5);
//        gmesh_TwoDTri->NodeVec()[5].SetCoord(coord);
//        Triangle_topology[2] = 5;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//        
//        }
//        
//        // 2nd element
//        
//        {
//
//            Triangle_topology[0] = 4;
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = 0.0; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[6].SetNodeId(2);
//            gmesh_TwoDTri->NodeVec()[6].SetCoord(coord);
//            Triangle_topology[1] = 6;
//            
//            Triangle_topology[2] = 0;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        
//        // 3rd element
//
//        {
//
//            Triangle_topology[0] = 4;
//            
//            Triangle_topology[1] = 0;
//            
//            Triangle_topology[2] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        
//        // 4th element
//        
//        {
//
//            Triangle_topology[0] = 5;
//            
//            Triangle_topology[1] = 2;
//            
//            Triangle_topology[2] = 1;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        
//        
//        // 5th element
//        
//        {
//
//            Triangle_topology[0] = 5;
//            
//            Triangle_topology[1] = 1;
//            
//            coord[0] = 0.0; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[7].SetNodeId(8);
//            gmesh_TwoDTri->NodeVec()[7].SetCoord(coord);
//            Triangle_topology[2] = 7;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        
//        
//        // 6th element
//        
//        {
// 
//            Triangle_topology[0] = 4;
//  
//            Triangle_topology[1] = 2;
//
//            Triangle_topology[2] = 5;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//        // 7th element
//        
//        {
// 
//            Triangle_topology[0] = 0;
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDTri->NodeVec()[8].SetNodeId(6);
//            gmesh_TwoDTri->NodeVec()[8].SetCoord(coord);
//            Triangle_topology[1] = 8;
//
//            Triangle_topology[2] = 1;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
//            elementid++;
//            
//        }
//    }
//   
//    // bottom
//    {
//        {
//            
//            Linear_topology[0] = 3;
//
//            Linear_topology[1] = 4;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom
//            elementid++;
//        }
//    
//        {
//
//            Linear_topology[0] = 4;
//            
//            Linear_topology[1] = 6;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom
//            elementid++;
//        }
//    }
//       // right
//    
//    {
//        {
//
//            Linear_topology[0] = 6;
//            
//            Linear_topology[1] = 0;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; right
//            elementid++;
//        }
//        
//        {
//
//            Linear_topology[0] = 0;
//
//            Linear_topology[1] = 8;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; right
//            elementid++;
//        }
//    }
//    // top
//    {
//        {
//
//            Linear_topology[0] = 8;
//            
//            Linear_topology[1] = 1;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDTri); // create boundary element; top
//            elementid++;
//        }
//        
//        {
//            
//            Linear_topology[0] = 1;
//            
//            Linear_topology[1] = 7;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDTri); // create boundary element; top
//            elementid++;
//        }
//    }
//      // left
//    {
//        {
//
//            Linear_topology[0] = 7;
//
//            Linear_topology[1] = 5;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; left
//            elementid++;
//        }
//        
//        {
//
//            Linear_topology[0] = 5;
//            
//            Linear_topology[1] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; left
//            elementid++;
//        }
//    }
//    
//    // Build the mesh
//    gmesh_TwoDTri->BuildConnectivity();
//    
//    std::ofstream outgmeshTwoDTri("geomesh_TwoDTri.txt");
//    gmesh_TwoDTri->Print(outgmeshTwoDTri);
//    
//    std::ofstream vtkgmeshTwoDTri("geomesh_TwoDTri.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDTri, vtkgmeshTwoDTri);
//    return gmesh_TwoDTri;
//    
//}
//
//
//
//// ************************************** Create 2D quadrilateral meshes ***************************************
//
//TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly)
//{
//    TPZGeoMesh * gmesh_TwoDQuad = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 2; // geometry dimension
//    
//    std::string name("geomesh TwoDQuad"); // geometry name
//    gmesh_TwoDQuad->SetName(name);
//    gmesh_TwoDQuad->SetDimension(geometry_dim);
//    
//    
//    gmesh_TwoDQuad->NodeVec().Resize(nnodesqu); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodesqu);
//    
//    TPZVec <long> Quadrilateral_topology(4);
//    TPZVec<long> Triangle_topology(3);
//    TPZVec <long> Linear_topology(2);
//    TPZVec<long> point_topology(1);
//
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//
//    long elementid = 0;
//    int physical_id = 1;
//    
//    // Index of boundary element
//    const int bc_bottom = -1; // define id for a material (border bottom)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_top = -3; // define id for a material (border top)
//    const int bc_left = -4; // define id for a material (border left)
//    
//    
//    {
//        // 0th element
//        {
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = Ly/2; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[0].SetNodeId(3);
//            gmesh_TwoDQuad->NodeVec()[0].SetCoord(coord);
//            Quadrilateral_topology[0] = 0; // index
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[1].SetNodeId(9);
//            gmesh_TwoDQuad->NodeVec()[1].SetCoord(coord);
//            Quadrilateral_topology[1] = 1;
//            
//            coord[0] = Lx/2; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[2].SetNodeId(8);
//            gmesh_TwoDQuad->NodeVec()[2].SetCoord(coord);
//            Quadrilateral_topology[2] = 2;
//            
//            coord[0] = Lx/2; // x coordinate
//            coord[1] = Ly/2; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[3].SetNodeId(4);
//            gmesh_TwoDQuad->NodeVec()[3].SetCoord(coord);
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
//            elementid++;
//            
//        }
//    
//        // 1st element
//        {
//            
//            coord[0] = Lx/2; // x coordinate
//            coord[1] = 0.0; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[4].SetNodeId(1);
//            gmesh_TwoDQuad->NodeVec()[4].SetCoord(coord);
//            Quadrilateral_topology[0] = 4; // index
//            
//
//            Quadrilateral_topology[1] = 3;
//            
//            coord[0] = 0.0; // x coordinate
//            coord[1] = Ly/2; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[5].SetNodeId(5);
//            gmesh_TwoDQuad->NodeVec()[5].SetCoord(coord);
//            Quadrilateral_topology[2] = 5;
//            
//            coord[0] = 0.0; // x coordinate
//            coord[1] = 0.0; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[6].SetNodeId(0);
//            gmesh_TwoDQuad->NodeVec()[6].SetCoord(coord);
//            Quadrilateral_topology[3] = 6;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
//            elementid++;
//            
//        }
//  
//        // 2nd element
//        {
//            
//            coord[0] = Lx; // x coordinate
//            coord[1] = 0.0; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[7].SetNodeId(2);
//            gmesh_TwoDQuad->NodeVec()[7].SetCoord(coord);
//            Quadrilateral_topology[0] = 7; // index
//            
//            
//            Quadrilateral_topology[1] = 0;
//            
//            Quadrilateral_topology[2] = 3;
//
//            Quadrilateral_topology[3] = 4;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
//            elementid++;
//            
//        }
//        // 3rd element
//        
//        {
//
//            Triangle_topology[0] = 3;
//            
//            Triangle_topology[1] = 2;
//            
//            coord[0] = Lx/4; // x coordinate
//            coord[1] =0.75*Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[8].SetNodeId(7);
//            gmesh_TwoDQuad->NodeVec()[8].SetCoord(coord);
//            Triangle_topology[2] = 8;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
//            elementid++;
//        }
//        
//        // 4th element
//        
//        {
//            
//            Triangle_topology[0] = 5;
//            
//            Triangle_topology[1] = 8;
//            
//            coord[0] = 0.0; // x coordinate
//            coord[1] = Ly; // Y coordinate
//            coord[2] = 0.0; // Z coordinate
//            gmesh_TwoDQuad->NodeVec()[9].SetNodeId(6);
//            gmesh_TwoDQuad->NodeVec()[9].SetCoord(coord);
//            Triangle_topology[2] = 9;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
//            elementid++;
//        }
//        
//        // 5th element
//        
//        {
//            
//            Triangle_topology[0] = 3;
//            
//            Triangle_topology[1] = 8;
//
//            Triangle_topology[2] = 5;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
//            elementid++;
//        }
//        
//        // 6th element
//        
//        {
//            
//            Triangle_topology[0] = 8;
//            
//            Triangle_topology[1] = 9;
//            
//            Triangle_topology[2] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDQuad); // create triangle element
//            elementid++;
//        }
//        
//        // 7th element point
//        
//        {
//
//            point_topology[0] = 3;
////            new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (elementid, point_topology, physical_id, *gmesh_TwoDQuad); //create point element
//            elementid++;
//
//        }
//
//    }
//        
//        // bottom
//        {
//            {
//                
//                Linear_topology[0] = 6;
//                
//                Linear_topology[1] = 4;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom
//                elementid++;
//            }
//            
//            {
//                
//                Linear_topology[0] = 4;
//                
//                Linear_topology[1] = 7;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom
//                elementid++;
//            }
//        }
//        // right
//        
//        {
//            {
//                
//                Linear_topology[0] = 7;
//                
//                Linear_topology[1] = 0;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; right
//                elementid++;
//            }
//            
//            {
//                
//                Linear_topology[0] = 0;
//                
//                Linear_topology[1] = 1;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; right
//                elementid++;
//            }
//        }
//        // top
//        {
//            {
//                
//                Linear_topology[0] = 1;
//                
//                Linear_topology[1] = 2;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; top
//                elementid++;
//            }
//            
//            {
//                
//                Linear_topology[0] = 2;
//                
//                Linear_topology[1] = 9;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; top
//                elementid++;
//            }
//        }
//        // left
//        {
//            {
//                
//                Linear_topology[0] = 9;
//                
//                Linear_topology[1] = 5;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; left
//                elementid++;
//            }
//            
//            {
//                
//                Linear_topology[0] = 5;
//                
//                Linear_topology[1] = 6;
//                
////                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; left
//                elementid++;
//            }
//        }
//
//        
//        // Build the mesh
//        gmesh_TwoDQuad->BuildConnectivity();
//        
//    std::ofstream outgmeshTwoDQuad("geomesh_TwoDQuad.txt");
//    gmesh_TwoDQuad->Print(outgmeshTwoDQuad);
//    
//    std::ofstream vtkgmeshTwoDQuad("geomesh_TwoDQuad.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDQuad, vtkgmeshTwoDQuad);
//
//    
//        return gmesh_TwoDQuad;
//        
//    }
//   
//
//
//// ************************************** Create 3D hexahedral and prism meshes ***************************************
//
//
//TPZGeoMesh *CreateThreeDHexPriGMesh(long nnodesthr, REAL Lx, REAL Ly, REAL Lz)
//{
//    TPZGeoMesh * gmesh_ThreeDHexPri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 3; // geometry dimension
//    
//    std::string name("geomesh ThreeDHexPri"); // geometry name
//    gmesh_ThreeDHexPri->SetName(name);
//    gmesh_ThreeDHexPri->SetDimension(geometry_dim);
//    
//    
//    gmesh_ThreeDHexPri->NodeVec().Resize(nnodesthr); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodesthr);
//    
//    TPZVec<long> Hexahedron_topology(8);
//   
//    TPZVec<long> Prism_topology(6);
//    
//    TPZVec <long> Quadrilateral_topology(4);
//    TPZVec<long> Triangle_topology(3);
//    TPZVec <long> Linear_topology(2);
//    TPZVec<long> point_topology(1);
//    
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//    
//    long elementid = 0;
//    int physical_id = 1;
//    
//    // Index of boundary element
//    const int bc_front = -1; // define id for a material (border in front)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_back = -3; // define id for a material (border back)
//    const int bc_left = -4; // define id for a material (border left)
//    const int bc_bottom = -5; // define id for a material (border bottom)
//    const int bc_top = -6; // define id for a material (border top)
//    
//    // 0th element
//    
//    {
//        
//        // 0th node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[0].SetNodeId(0);
//        gmesh_ThreeDHexPri->NodeVec()[0].SetCoord(coord);
//        Hexahedron_topology[0] = 0;
//        
//        // 1st node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[1].SetNodeId(1);
//        gmesh_ThreeDHexPri->NodeVec()[1].SetCoord(coord);
//        Hexahedron_topology[1] = 1;
//        
//        // 2nd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[2].SetNodeId(2);
//        gmesh_ThreeDHexPri->NodeVec()[2].SetCoord(coord);
//        Hexahedron_topology[2] = 2;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[3].SetNodeId(3);
//        gmesh_ThreeDHexPri->NodeVec()[3].SetCoord(coord);
//        Hexahedron_topology[3] = 3;
//        
//        // 4th node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[4].SetNodeId(4);
//        gmesh_ThreeDHexPri->NodeVec()[4].SetCoord(coord);
//        Hexahedron_topology[4] = 4;
//        
//        // 5th node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[5].SetNodeId(5);
//        gmesh_ThreeDHexPri->NodeVec()[5].SetCoord(coord);
//        Hexahedron_topology[5] = 5;
//        
//        // 6th node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[6].SetNodeId(6);
//        gmesh_ThreeDHexPri->NodeVec()[6].SetCoord(coord);
//        Hexahedron_topology[6] = 6;
//        
//        // 7th node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[7].SetNodeId(7);
//        gmesh_ThreeDHexPri->NodeVec()[7].SetCoord(coord);
//        Hexahedron_topology[7] = 7;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//
//    }
//    
//    // 1st element
//    
//    {
//        
//        // 1st node
//        Hexahedron_topology[0] = 4;
//        
//        // 2nd node
//        Hexahedron_topology[1] = 5;
//        
//        // 3rd node
//        Hexahedron_topology[2] = 6;
//        
//        // 4th node
//        Hexahedron_topology[3] = 7;
//        
//        // 5th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[8].SetNodeId(8);
//        gmesh_ThreeDHexPri->NodeVec()[8].SetCoord(coord);
//        Hexahedron_topology[4] = 8;
//        
//        // 6th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[9].SetNodeId(9);
//        gmesh_ThreeDHexPri->NodeVec()[9].SetCoord(coord);
//        Hexahedron_topology[5] = 9;
//        
//        // 7th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[10].SetNodeId(10);
//        gmesh_ThreeDHexPri->NodeVec()[10].SetCoord(coord);
//        Hexahedron_topology[6] = 10;
//        
//        // 8th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[11].SetNodeId(11);
//        gmesh_ThreeDHexPri->NodeVec()[11].SetCoord(coord);
//        Hexahedron_topology[7] = 11;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//
//    }
//    
//    
//    // 2nd element
//    
//    {
//        
//        // 1st node
//        Hexahedron_topology[0] = 3;
//        
//        // 2nd node
//        Hexahedron_topology[1] = 2;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[12].SetNodeId(12);
//        gmesh_ThreeDHexPri->NodeVec()[12].SetCoord(coord);
//        Hexahedron_topology[2] = 12;
//        
//        // 4th node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[13].SetNodeId(13);
//        gmesh_ThreeDHexPri->NodeVec()[13].SetCoord(coord);
//        Hexahedron_topology[3] = 13;
//        
//        // 5th node
//        Hexahedron_topology[4] = 11;
//        
//        // 6th node
//        Hexahedron_topology[5] = 10;
//        
//        // 7th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[14].SetNodeId(14);
//        gmesh_ThreeDHexPri->NodeVec()[14].SetCoord(coord);
//        Hexahedron_topology[6] = 14;
//        
//        // 8th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[15].SetNodeId(15);
//        gmesh_ThreeDHexPri->NodeVec()[15].SetCoord(coord);
//        Hexahedron_topology[7] = 15;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//        
//    }
//    
//    // 3rd element
//    
//    {
//        
//        // 1st node
//        Hexahedron_topology[0] = 1;
//        
//        // 2nd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[16].SetNodeId(16);
//        gmesh_ThreeDHexPri->NodeVec()[16].SetCoord(coord);
//        Hexahedron_topology[1] = 16;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[17].SetNodeId(17);
//        gmesh_ThreeDHexPri->NodeVec()[17].SetCoord(coord);
//        Hexahedron_topology[2] = 17;
//        
//        // 4th node
//        Hexahedron_topology[3] = 2;
//        
//        // 5th node
//        Hexahedron_topology[4] = 9;
//        
//        // 6th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[18].SetNodeId(18);
//        gmesh_ThreeDHexPri->NodeVec()[18].SetCoord(coord);
//        Hexahedron_topology[5] = 18;
//        
//        // 7th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[19].SetNodeId(19);
//        gmesh_ThreeDHexPri->NodeVec()[19].SetCoord(coord);
//        Hexahedron_topology[6] = 19;
//        
//        // 8th node
//        Hexahedron_topology[7] = 10;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//        
//    }
//    
//        // 4th element
//    {
//        
//        // 1st node
//        Prism_topology[0] = 17;
//        
//        // 2nd node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[20].SetNodeId(20);
//        gmesh_ThreeDHexPri->NodeVec()[20].SetCoord(coord);
//        Prism_topology[1] = 20;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[21].SetNodeId(21);
//        gmesh_ThreeDHexPri->NodeVec()[21].SetCoord(coord);
//        Prism_topology[2] = 21;
//        
//        // 4th node
//        Prism_topology[3] = 2;
//        
//        // 5th node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[22].SetNodeId(22);
//        gmesh_ThreeDHexPri->NodeVec()[22].SetCoord(coord);
//        Prism_topology[4] = 22;
//        
//        // 6th node
//        Prism_topology[5] = 12;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (elementid, Prism_topology, physical_id, *gmesh_ThreeDHexPri);
//    }
//    
//    // 5th element
//    {
//        
//        // 1st node
//        Prism_topology[0] = 17;
//        
//        // 2nd node
//        coord[0] = Lx/2; // x coordinate
//        coord[1] = Ly/2; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[23].SetNodeId(23);
//        gmesh_ThreeDHexPri->NodeVec()[23].SetCoord(coord);
//        Prism_topology[1] = 23;
//        
//        // 3rd node
//        Prism_topology[2] = 20;
//        
//        // 4th node
//        Prism_topology[3] = 2;
//        
//        // 5th node
//        Prism_topology[4] = 6;
//        
//        // 6th node
//        Prism_topology[5] = 22;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (elementid, Prism_topology, physical_id, *gmesh_ThreeDHexPri);
//    }
//    
//    
//    // 6th element
//    {
//        
//        // 1st node
//        Prism_topology[0] = 23;
//        
//        // 2nd node
//        Prism_topology[1] = 19;
//        
//        // 3rd node
//        Prism_topology[2] = 20;
//        
//        // 4th node
//        Prism_topology[3] = 6;
//        
//        // 5th node
//        Prism_topology[4] = 10;
//        
//        // 6th node
//        Prism_topology[5] = 22;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (elementid, Prism_topology, physical_id, *gmesh_ThreeDHexPri);
//    }
//    
//    // 7th element
//    {
//        
//        // 1st node
//        Prism_topology[0] = 20;
//        
//        // 2nd node
//        Prism_topology[1] = 19;
//        
//        // 3rd node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPri->NodeVec()[24].SetNodeId(24);
//        gmesh_ThreeDHexPri->NodeVec()[24].SetCoord(coord);
//        Prism_topology[2] = 24;
//        
//        // 4th node
//        Prism_topology[3] = 22;
//        
//        // 5th node
//        Prism_topology[4] = 10;
//        
//        // 6th node
//        Prism_topology[5] = 14;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (elementid, Prism_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//
//    }
//    
//    // 8th element point
//    
//    {
//        
//        point_topology[0] = 18;
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (elementid, point_topology, physical_id, *gmesh_ThreeDHexPri);
//        elementid++;
//        
//    }
//    
//    // 9th element line
//
//    {
//        
//        Linear_topology[0] = 18;
//        
//        Linear_topology[1] = 16;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_ThreeDHexPri);       elementid++;
//    }
//    
//
//    // ********************* boundary ****************
//
//      // in front
//    {
//        {
//            
//            Quadrilateral_topology[0] = 0;
//            
//            Quadrilateral_topology[1] = 1;
//
//            Quadrilateral_topology[2] = 2;
//
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPri); // create boundary element; in front
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 1;
//            
//            Quadrilateral_topology[1] = 16;
//            
//            Quadrilateral_topology[2] = 17;
//            
//            Quadrilateral_topology[3] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPri); // create boundary element; in front
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 2;
//            
//            Quadrilateral_topology[1] = 17;
//            
//            Quadrilateral_topology[2] = 21;
//            
//            Quadrilateral_topology[3] = 12;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPri); // create boundary element; in front
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 3;
//            
//            Quadrilateral_topology[1] = 2;
//            
//            Quadrilateral_topology[2] = 12;
//            
//            Quadrilateral_topology[3] = 13;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPri); // create boundary element; in front
//            elementid++;
//            
//        }
//    }
//    
//    
//    // right
//    {
//        {
//            
//            Quadrilateral_topology[0] = 8;
//            
//            Quadrilateral_topology[1] = 9;
//            
//            Quadrilateral_topology[2] = 5;
//            
//            Quadrilateral_topology[3] = 4;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPri); // create boundary element; right
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 9;
//            
//            Quadrilateral_topology[1] = 18;
//            
//            Quadrilateral_topology[2] = 16;
//            
//            Quadrilateral_topology[3] = 1;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPri); // create boundary element; right
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 4;
//            
//            Quadrilateral_topology[1] = 5;
//            
//            Quadrilateral_topology[2] = 1;
//            
//            Quadrilateral_topology[3] = 0;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPri); // create boundary element; right
//            elementid++;
//            
//        }
//
//    }
//    // back
//    {
//        {
//            
//            Quadrilateral_topology[0] = 8;
//            
//            Quadrilateral_topology[1] = 9;
//            
//            Quadrilateral_topology[2] = 10;
//            
//            Quadrilateral_topology[3] = 11;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPri); // create boundary element; back
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 9;
//            
//            Quadrilateral_topology[1] = 18;
//            
//            Quadrilateral_topology[2] = 19;
//            
//            Quadrilateral_topology[3] = 10;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPri); // create boundary element; back
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 10;
//            
//            Quadrilateral_topology[1] = 19;
//            
//            Quadrilateral_topology[2] = 24;
//            
//            Quadrilateral_topology[3] = 14;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPri); // create boundary element; back
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 11;
//            
//            Quadrilateral_topology[1] = 10;
//            
//            Quadrilateral_topology[2] = 14;
//            
//            Quadrilateral_topology[3] = 15;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPri); // create boundary element; back
//            elementid++;
//            
//        }
//    }
//
//    // left
//    {
//        {
//            
//            Quadrilateral_topology[0] = 15;
//            
//            Quadrilateral_topology[1] = 14;
//            
//            Quadrilateral_topology[2] = 12;
//            
//            Quadrilateral_topology[3] = 13;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPri); // create boundary element; left
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 14;
//            
//            Quadrilateral_topology[1] = 24;
//            
//            Quadrilateral_topology[2] = 20;
//            
//            Quadrilateral_topology[3] = 22;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPri); // create boundary element; left
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 22;
//            
//            Quadrilateral_topology[1] = 20;
//            
//            Quadrilateral_topology[2] = 21;
//            
//            Quadrilateral_topology[3] = 12;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPri); // create boundary element; left
//            elementid++;
//            
//        }
//        
//    }
//    // bottom
//    {
//        {
//            
//            Quadrilateral_topology[0] = 4;
//            
//            Quadrilateral_topology[1] = 8;
//            
//            Quadrilateral_topology[2] = 11;
//            
//            Quadrilateral_topology[3] = 7;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_bottom,*gmesh_ThreeDHexPri); // create boundary element; bottom
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 11;
//            
//            Quadrilateral_topology[1] = 15;
//            
//            Quadrilateral_topology[2] = 13;
//            
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_bottom,*gmesh_ThreeDHexPri); // create boundary element; bottom
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 0;
//            
//            Quadrilateral_topology[1] = 4;
//            
//            Quadrilateral_topology[2] = 7;
//            
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_bottom,*gmesh_ThreeDHexPri); // create boundary element; bottom
//            elementid++;
//            
//        }
//        
//    }
//    
//    // top
//    {
//        {
//            
//            Quadrilateral_topology[0] = 16;
//            
//            Quadrilateral_topology[1] = 18;
//            
//            Quadrilateral_topology[2] = 19;
//            
//            Quadrilateral_topology[3] = 17;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_top,*gmesh_ThreeDHexPri); // create boundary element; top
//            elementid++;
//            
//        }
//        {
//            
//            Triangle_topology[0] = 23;
//            
//            Triangle_topology[1] = 19;
//
//            Triangle_topology[2] = 20;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,bc_top,*gmesh_ThreeDHexPri); // create boundary element; top
//            elementid++;
//        }
//        {
//            
//            Triangle_topology[0] = 19;
//            
//            Triangle_topology[1] = 24;
//            
//            Triangle_topology[2] = 20;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,bc_top,*gmesh_ThreeDHexPri); // create boundary element; top
//            elementid++;
//        }
//        {
//            
//            Triangle_topology[0] = 17;
//            
//            Triangle_topology[1] = 23;
//            
//            Triangle_topology[2] = 20;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,bc_top,*gmesh_ThreeDHexPri); // create boundary element; top
//            elementid++;
//        }
//        {
//            
//            Triangle_topology[0] = 17;
//            
//            Triangle_topology[1] = 20;
//            
//            Triangle_topology[2] = 21;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,bc_top,*gmesh_ThreeDHexPri); // create boundary element; top
//            elementid++;
//        }
//
//        
//    }
//
//    
//    // Build the mesh
//    gmesh_ThreeDHexPri->BuildConnectivity();
//    
//    
//    std::ofstream outgmeshThreeDHexPri("geomesh_ThreeDHexPri.txt");
//    gmesh_ThreeDHexPri->Print(outgmeshThreeDHexPri);
//    
//    std::ofstream vtkgmeshThreeDHexPri("geomesh_ThreeDHexPri.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDHexPri, vtkgmeshThreeDHexPri);
//        return gmesh_ThreeDHexPri;
//
//    }
//
//// ----------------------------------------------------------------------------------------------------
//
//
//// ************************************** Create 3D hexahedral and prism meshes ***************************************
//
//
//TPZGeoMesh *CreateThreeDHexPrytetGMesh(long nnodesthrhpt, REAL Lx, REAL Ly, REAL Lz)
//{
//    TPZGeoMesh * gmesh_ThreeDHexPrytet = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 3; // geometry dimension
//    
//    std::string name("geomesh ThreeDHexPrytet"); // geometry name
//    gmesh_ThreeDHexPrytet->SetName(name);
//    gmesh_ThreeDHexPrytet->SetDimension(geometry_dim);
//    
//    
//    gmesh_ThreeDHexPrytet->NodeVec().Resize(nnodesthrhpt); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodesthrhpt);
//    
//    TPZVec<long> Hexahedron_topology(8);
//    TPZVec <long> Quadrilateral_topology(4);
//    
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//    
//    long elementid = 0;
//    int physical_id = 1;
//    
//    // Index of boundary element
//    const int bc_front = -1; // define id for a material (border in front)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_back = -3; // define id for a material (border back)
//    const int bc_left = -4; // define id for a material (border left)
//    const int bc_bottom = -5; // define id for a material (border bottom)
//    const int bc_top = -6; // define id for a material (border top)
//    
//    // 0th element
//    
//    {
//        
//        // 0th node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[0].SetNodeId(0);
//        gmesh_ThreeDHexPrytet->NodeVec()[0].SetCoord(coord);
//        Hexahedron_topology[0] = 0;
//        
//        // 1st node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[1].SetNodeId(1);
//        gmesh_ThreeDHexPrytet->NodeVec()[1].SetCoord(coord);
//        Hexahedron_topology[1] = 1;
//        
//        // 2nd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[2].SetNodeId(2);
//        gmesh_ThreeDHexPrytet->NodeVec()[2].SetCoord(coord);
//        Hexahedron_topology[2] = 2;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[3].SetNodeId(3);
//        gmesh_ThreeDHexPrytet->NodeVec()[3].SetCoord(coord);
//        Hexahedron_topology[3] = 3;
//        
//        // 4th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[4].SetNodeId(4);
//        gmesh_ThreeDHexPrytet->NodeVec()[4].SetCoord(coord);
//        Hexahedron_topology[4] = 4;
//        
//        // 5th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[5].SetNodeId(5);
//        gmesh_ThreeDHexPrytet->NodeVec()[5].SetCoord(coord);
//        Hexahedron_topology[5] = 5;
//        
//        // 6th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz/2; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[6].SetNodeId(6);
//        gmesh_ThreeDHexPrytet->NodeVec()[6].SetCoord(coord);
//        Hexahedron_topology[6] = 6;
//        
//        // 7th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[7].SetNodeId(7);
//        gmesh_ThreeDHexPrytet->NodeVec()[7].SetCoord(coord);
//        Hexahedron_topology[7] = 7;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
//        elementid++;
//        
//    }
//    
//    // 1st element
//    
//    {
//        
//        // 1st node
//        Hexahedron_topology[0] = 1;
//        
//        // 2nd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[8].SetNodeId(8);
//        gmesh_ThreeDHexPrytet->NodeVec()[8].SetCoord(coord);
//        Hexahedron_topology[1] = 8;
//        
//        // 3rd node
//        coord[0] = Lx; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[9].SetNodeId(9);
//        gmesh_ThreeDHexPrytet->NodeVec()[9].SetCoord(coord);
//        Hexahedron_topology[2] = 9;
//        
//        // 4th node
//        Hexahedron_topology[3] = 2;
//        
//        // 5th node
//        Hexahedron_topology[4] = 5;
//        
//        // 6th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = Ly; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[10].SetNodeId(10);
//        gmesh_ThreeDHexPrytet->NodeVec()[10].SetCoord(coord);
//        Hexahedron_topology[5] = 10;
//        
//        // 7th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = Lz; // Z coordinate
//        gmesh_ThreeDHexPrytet->NodeVec()[11].SetNodeId(11);
//        gmesh_ThreeDHexPrytet->NodeVec()[11].SetCoord(coord);
//        Hexahedron_topology[6] = 11;
//        
//        // 8th node
//        Hexahedron_topology[7] = 6;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexPrytet);
//        elementid++;
//        
//    }
//    
//    // ********************* boundary ****************
//    
//    // in front
//    {
//        {
//            
//            Quadrilateral_topology[0] = 0;
//            
//            Quadrilateral_topology[1] = 1;
//            
//            Quadrilateral_topology[2] = 2;
//            
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPrytet); // create boundary element; in front
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 1;
//            
//            Quadrilateral_topology[1] = 8;
//            
//            Quadrilateral_topology[2] = 9;
//            
//            Quadrilateral_topology[3] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_front,*gmesh_ThreeDHexPrytet); // create boundary element; in front
//            elementid++;
//            
//        }
//        
//    }
//    
//
//    // right
//    {
//        {
//            
//            Quadrilateral_topology[0] = 4;
//            
//            Quadrilateral_topology[1] = 5;
//            
//            Quadrilateral_topology[2] = 1;
//            
//            Quadrilateral_topology[3] = 0;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPrytet); // create boundary element; right
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 5;
//            
//            Quadrilateral_topology[1] = 10;
//            
//            Quadrilateral_topology[2] = 8;
//            
//            Quadrilateral_topology[3] = 1;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_right,*gmesh_ThreeDHexPrytet); // create boundary element; right
//            elementid++;
//            
//        }
//        
//    }
//    
//    // back
//    {
//        {
//            
//            Quadrilateral_topology[0] = 4;
//            
//            Quadrilateral_topology[1] = 5;
//            
//            Quadrilateral_topology[2] = 6;
//            
//            Quadrilateral_topology[3] = 7;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPrytet); // create boundary element; back
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 5;
//            
//            Quadrilateral_topology[1] = 10;
//            
//            Quadrilateral_topology[2] = 11;
//            
//            Quadrilateral_topology[3] = 6;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_back,*gmesh_ThreeDHexPrytet); // create boundary element; back
//            elementid++;
//            
//        }
//        
//    }
//    
//    // left
//    {
//        {
//            
//            Quadrilateral_topology[0] = 7;
//            
//            Quadrilateral_topology[1] = 6;
//            
//            Quadrilateral_topology[2] = 2;
//            
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPrytet); // create boundary element; left
//            elementid++;
//            
//        }
//        {
//            
//            Quadrilateral_topology[0] = 6;
//            
//            Quadrilateral_topology[1] = 11;
//            
//            Quadrilateral_topology[2] = 9;
//            
//            Quadrilateral_topology[3] = 2;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_left,*gmesh_ThreeDHexPrytet); // create boundary element; left
//            elementid++;
//            
//        }
//        
//    }
//    
//    
//    
//    // bottom
//
//    {
//            
//            Quadrilateral_topology[0] = 0;
//            
//            Quadrilateral_topology[1] = 4;
//            
//            Quadrilateral_topology[2] = 7;
//            
//            Quadrilateral_topology[3] = 3;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_bottom,*gmesh_ThreeDHexPrytet); // create boundary element; bottom
//            elementid++;
//            
//   }
//    
//    
//    // top
//   {
//            
//            Quadrilateral_topology[0] = 8;
//            
//            Quadrilateral_topology[1] = 10;
//            
//            Quadrilateral_topology[2] = 11;
//            
//            Quadrilateral_topology[3] = 9;
//            
////            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,bc_top,*gmesh_ThreeDHexPrytet); // create boundary element; top
//            elementid++;
//            
//   }
//
//    
//    // Build the mesh
//    gmesh_ThreeDHexPrytet->BuildConnectivity();
//    std::ofstream outgmeshThreeDHexPrytet("geomesh_ThreeDHexPrytet.txt");
//    gmesh_ThreeDHexPrytet->Print(outgmeshThreeDHexPrytet);
//    
//    std::ofstream vtkgmeshThreeDHexPrytet("geomesh_ThreeDHexPrytet.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDHexPrytet, vtkgmeshThreeDHexPrytet);
//    
//    return gmesh_ThreeDHexPrytet;
//    
//}
//
//
//
//// *************************** (nonlinear elements) ****************************************************
//
//
//void NonOneDElements(){
//    
//    TPZGeoMesh geometry_1DNonL; // Create the objet that will describe the nonlinear geometry (Line).
//    int n_nodes = 3; // number of nodes
//    int geometry_dim = 1; // geometry dimension
//    std::string name("geometry 1DNonL"); // geometry name
//    
//    // setting the object
//    geometry_1DNonL.SetName(name);
//    geometry_1DNonL.SetDimension(geometry_dim);
//    geometry_1DNonL.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonLinear_topology(3);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        // 1st node of line element is located at x ={-10,PI,0.5}
//        
//        coord[0] = -10.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_1DNonL.NodeVec()[0].SetNodeId(node_id);
//        geometry_1DNonL.NodeVec()[0].SetCoord(coord);
//        NonLinear_topology[0] = node_id;
//        
//        // 2nd node of line element is located at x ={10,PI,0.5}
//        node_id++;
//        
//        coord[0] = 10.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_1DNonL.NodeVec()[1].SetNodeId(node_id);
//        geometry_1DNonL.NodeVec()[1].SetCoord(coord);
//        NonLinear_topology[1] = node_id;
//        
//        // 3rd node of line element is located at x ={0.0,PI,0.5}
//        node_id++;
//
//        coord[0] = 0.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_1DNonL.NodeVec()[2].SetNodeId(node_id);
//        geometry_1DNonL.NodeVec()[2].SetCoord(coord);
//        NonLinear_topology[2] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine> (element_id, NonLinear_topology, physical_id, geometry_1DNonL);
//    }
//    
//    geometry_1DNonL.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_1DNonL.txt");
//    geometry_1DNonL.Print(file);
//    
//}
//
//
//// ------------------------------------------------------------------------------------------
//
//void NonTwoDTriElements() {
//    
//    TPZGeoMesh geometry_Non2DTri; // Create the objet that will describe the nonlinear geometry (2D Triangle).
//    int n_nodes = 6; // number of nodes
//    int geometry_dim = 2; // geometry dimension
//    std::string name("geometry Non2DTri"); // geometry name
//    
//    // setting the object
//    geometry_Non2DTri.SetName(name);
//    geometry_Non2DTri.SetDimension(geometry_dim);
//    geometry_Non2DTri.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonTriangle_topology(6);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of triangle element is located at x ={PI,0.0,0.0}
//        
//        coord[0] = M_PI; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DTri.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[0].SetCoord(coord);
//        NonTriangle_topology[0] = node_id;
//        
//        // 2nd node of triangle element is located at x ={0.0,PI,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DTri.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[1].SetCoord(coord);
//        NonTriangle_topology[1] = node_id;
//        
//        // 3rd node of triangle element is located at x ={0.0,0.0,PI}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non2DTri.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[2].SetCoord(coord);
//        NonTriangle_topology[2] = node_id;
//        
//        
//        // 4th node of triangle element is located at x ={PI/2,PI/2,0.0}
//        node_id++;
//
//        coord[0] = M_PI/2; // x coordinate
//        coord[1] = M_PI/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DTri.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[3].SetCoord(coord);
//        NonTriangle_topology[3] = node_id;
//        
//        
//        // 5th node of triangle element is located at x ={0.0,PI/2,PI/2}
//        node_id++;
//        
//        coord[0] = M_PI/2; // x coordinate
//        coord[1] = M_PI/2; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DTri.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[4].SetCoord(coord);
//        NonTriangle_topology[4] = node_id;
//        
//        // 6th node of triangle element is located at x ={PI/2,0.0,PI/2}
//        node_id++;
//        
//        coord[0] = M_PI/2; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI/2; // Z coordinate
//        geometry_Non2DTri.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non2DTri.NodeVec()[5].SetCoord(coord);
//        NonTriangle_topology[5] = node_id;
//        
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticTrig> (element_id, NonTriangle_topology, physical_id, geometry_Non2DTri);
//    }
//    
//    geometry_Non2DTri.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non2DTri.txt");
//    geometry_Non2DTri.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void NonTwoDQuadElements() {
//    
//    TPZGeoMesh geometry_Non2DQuad; // Create the objet that will describe the nonlinear geometry (2D Quadrilateral).
//    int n_nodes = 8; // number of nodes
//    int geometry_dim = 2; // geometry dimension
//    std::string name("geometry Non2DQuad"); // geometry name
//    
//    // setting the object
//    geometry_Non2DQuad.SetName(name);
//    geometry_Non2DQuad.SetDimension(geometry_dim);
//    geometry_Non2DQuad.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonQuadrilateral_topology(8);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of quadrilateral element is located at x ={-1.0,PI,-1.0}
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[0].SetCoord(coord);
//        NonQuadrilateral_topology[0] = node_id;
//        
//        // 2nd node of quadrilateral element is located at x ={1.0,PI,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[1].SetCoord(coord);
//        NonQuadrilateral_topology[1] = node_id;
//        
//        // 3rd node of quadrilateral element is located at x ={1.0,PI,1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[2].SetCoord(coord);
//        NonQuadrilateral_topology[2] = node_id;
//        
//        // 4th node of quadrilateral element is located at x ={-1.0,PI,1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[3].SetCoord(coord);
//        NonQuadrilateral_topology[3] = node_id;
//        
//        // 5th node of quadrilateral element is located at x ={0.0,PI,-1.0}
//        node_id++;
//
//        coord[0] = 0.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[4].SetCoord(coord);
//        NonQuadrilateral_topology[4] = node_id;
//        
//        // 6th node of quadrilateral element is located at x ={1.0,PI,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[5].SetCoord(coord);
//        NonQuadrilateral_topology[5] = node_id;
//        
//        // 7th node of quadrilateral element is located at x ={0.0,PI,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[6].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[6].SetCoord(coord);
//        NonQuadrilateral_topology[6] = node_id;
//        
//        // 8th node of quadrilateral element is located at x ={-1.0,PI,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = M_PI; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non2DQuad.NodeVec()[7].SetNodeId(node_id);
//        geometry_Non2DQuad.NodeVec()[7].SetCoord(coord);
//        NonQuadrilateral_topology[7] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (element_id, NonQuadrilateral_topology, physical_id, geometry_Non2DQuad);
//    }
//    
//    geometry_Non2DQuad.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non2DQuad.txt");
//    geometry_Non2DQuad.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void NonThreeDTetraElements() {
//    
//    TPZGeoMesh geometry_Non3DTetra; // Create the objet that will describe the nonlinear geometry (3D Tetrahedron).
//    int n_nodes = 10; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry Non3DTetra"); // geometry name
//    
//    // setting the object
//    geometry_Non3DTetra.SetName(name);
//    geometry_Non3DTetra.SetDimension(geometry_dim);
//    geometry_Non3DTetra.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonTetrahedron_topology(10);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of tetrahedron element is located at x ={0.0,0.0,0.0}
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[0].SetCoord(coord);
//        NonTetrahedron_topology[0] = node_id;
//        
//        // 2nd node of tetrahedron element is located at x ={1/2,0.0,0.0}
//        node_id++;
//
//        coord[0] = 0.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[1].SetCoord(coord);
//        NonTetrahedron_topology[1] = node_id;
//        
//        // 3rd node of tetrahedron element is located at x ={1.0,0.0,0.0}
//        node_id++;
//
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[2].SetCoord(coord);
//        NonTetrahedron_topology[2] = node_id;
//        
//        // 4th node of tetrahedron element is located at x ={0.0,1/2,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[3].SetCoord(coord);
//        NonTetrahedron_topology[3] = node_id;
//        
//        // 5th node of tetrahedron element is located at x ={1/2,1/2,0.0}
//        node_id++;
//        
//        coord[0] = 0.5; // x coordinate
//        coord[1] = 0.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[4].SetCoord(coord);
//        NonTetrahedron_topology[4] = node_id;
//        
//        // 6th node of tetrahedron element is located at x ={0.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[5].SetCoord(coord);
//        NonTetrahedron_topology[5] = node_id;
//        
//        // 7th node of tetrahedron element is located at x ={0.0,0.0,1/2}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[6].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[6].SetCoord(coord);
//        NonTetrahedron_topology[6] = node_id;
//        
//        // 8th node of tetrahedron element is located at x ={1/2,0.0,1/2}
//        node_id++;
//        
//        coord[0] = 0.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[7].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[7].SetCoord(coord);
//        NonTetrahedron_topology[7] = node_id;
//        
//        // 9th node of tetrahedron element is located at x ={0.0,1/2,1/2}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.5; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[8].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[8].SetCoord(coord);
//        NonTetrahedron_topology[8] = node_id;
//        
//        
//        // 10th node of tetrahedron element is located at x ={0.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DTetra.NodeVec()[9].SetNodeId(node_id);
//        geometry_Non3DTetra.NodeVec()[9].SetCoord(coord);
//        NonTetrahedron_topology[9] = node_id;
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticTetra> (element_id, NonTetrahedron_topology, physical_id, geometry_Non3DTetra);
//    }
//    
//    geometry_Non3DTetra.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non3DTetra.txt");
//    geometry_Non3DTetra.Print(file);
//    
//    std::ofstream filevtk("geometry_Non3DTetra.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(&geometry_Non3DTetra, filevtk);
//    
//    
//}
//
//
//// ------------------------------------------------------------------------------------------
//// This one needs to fix
//
//void NonThreeDPyraElements() {
//    
//    TPZGeoMesh geometry_Non3DPyra; // Create the objet that will describe the nonlinear geometry (3D Pyramid).
//    int n_nodes = 13; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry Non3DPyra"); // geometry name
//    
//    // setting the object
//    geometry_Non3DPyra.SetName(name);
//    geometry_Non3DPyra.SetDimension(geometry_dim);
//    geometry_Non3DPyra.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonPyramid_topology(13);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 0th node of pyramid element is located at x ={1.0,-1.0,0.0}
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[0].SetCoord(coord);
//        NonPyramid_topology[0] = node_id;
//        
//        // 1st node of pyramid element is located at x ={1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[1].SetCoord(coord);
//        NonPyramid_topology[1] = node_id;
//        
//        // 2nd node of pyramid element is located at x ={-1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[2].SetCoord(coord);
//        NonPyramid_topology[2] = node_id;
//        
//        // 3th node of pyramid element is located at x ={-1.0,-1.0,0.0}
//        node_id++;
//
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[3].SetCoord(coord);
//        NonPyramid_topology[3] = node_id;
//        
//        // 4th node of pyramid element is located at x ={0.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[4].SetCoord(coord);
//        NonPyramid_topology[4] = node_id;
//        
//        // 5th node of pyramid element is located at x ={1.0,0.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[5].SetCoord(coord);
//        NonPyramid_topology[5] = node_id;
//        
//        // 6th node of pyramid element is located at x ={0.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[6].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[6].SetCoord(coord);
//        NonPyramid_topology[6] = node_id;
//        
//        // 7th node of pyramid element is located at x ={-1.0,0.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[7].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[7].SetCoord(coord);
//        NonPyramid_topology[7] = node_id;
//        
//        
//        // 8th node of pyramid element is located at x ={0.0,-1.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[8].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[8].SetCoord(coord);
//        NonPyramid_topology[8] = node_id;
//        
//        // 9th node of pyramid element is located at x ={0.5,-0.5,0.5}
//        node_id++;
//        
//        coord[0] = 0.5; // x coordinate
//        coord[1] = -0.5; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[9].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[9].SetCoord(coord);
//        NonPyramid_topology[9] = node_id;
//        
//        // 10th node of pyramid element is located at x ={0.5,0.5,0.5}
//        node_id++;
//        
//        coord[0] = 0.5; // x coordinate
//        coord[1] = 0.5; // Y coordinate
//        coord[2] = -0.5; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[10].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[10].SetCoord(coord);
//        NonPyramid_topology[10] = node_id;
//        
//        // 11th node of pyramid element is located at x ={-0.5,0.5,0.5}
//        node_id++;
//        
//        coord[0] = -0.5; // x coordinate
//        coord[1] = 0.5; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[11].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[11].SetCoord(coord);
//        NonPyramid_topology[11] = node_id;
//        
//        
//        // 12th node of pyramid element is located at x ={-0.5,-0.5,0.5}
//        node_id++;
//        
//        coord[0] = -0.5; // x coordinate
//        coord[1] = -0.5; // Y coordinate
//        coord[2] = 0.5; // Z coordinate
//        geometry_Non3DPyra.NodeVec()[12].SetNodeId(node_id);
//        geometry_Non3DPyra.NodeVec()[12].SetCoord(coord);
//        NonPyramid_topology[12] = node_id;
//        
//        
//        
////        new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (element_id, NonPyramid_topology, physical_id, geometry_Non3DPyra);
//    }
//    
//    geometry_Non3DPyra.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non3DPyra.txt");
//    geometry_Non3DPyra.Print(file);
//    
//}
//
//// ------------------------------------------------------------------------------------------
//
//
//void NonThreeDPrisElements() {
//    
//    TPZGeoMesh geometry_Non3DPris; // Create the objet that will describe the nonlinear geometry (3D Prism).
//    int n_nodes = 15; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry Non3DPris"); // geometry name
//    
//    // setting the object
//    geometry_Non3DPris.SetName(name);
//    geometry_Non3DPris.SetDimension(geometry_dim);
//    geometry_Non3DPris.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonPrism_topology(15);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 0th node of prism element is located at x ={0.0,1.0,1.0}
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPris.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[0].SetCoord(coord);
//        NonPrism_topology[0] = node_id;
//        
//        // 1st node of prism element is located at x ={0.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPris.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[1].SetCoord(coord);
//        NonPrism_topology[1] = node_id;
//        
//        // 2nd node of prism element is located at x ={1.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DPris.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[2].SetCoord(coord);
//        NonPrism_topology[2] = node_id;
//        
//        // 3rd node of prism element is located at x ={0.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[3].SetCoord(coord);
//        NonPrism_topology[3] = node_id;
//        
//        // 4th node of prism element is located at x ={0.0,0.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[4].SetCoord(coord);
//        NonPrism_topology[4] = node_id;
//        
//        // 5th node of prism element is located at x ={1.0,0.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[5].SetCoord(coord);
//        NonPrism_topology[5] = node_id;
//        
//        // 6th node of prism element is located at x ={0.0,0.5,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[6].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[6].SetCoord(coord);
//        NonPrism_topology[6] = node_id;
//        
//        // 7th node of prism element is located at x ={0.5,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[7].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[7].SetCoord(coord);
//        NonPrism_topology[7] = node_id;
//        
//        // 8th node of prism element is located at x ={0.5,0.5,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[8].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[8].SetCoord(coord);
//        NonPrism_topology[8] = node_id;
//        
//        // 9th node of prism element is located at x ={0.0,0.5,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[9].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[9].SetCoord(coord);
//        NonPrism_topology[9] = node_id;
//        
//        
//        // 10th node of prism element is located at x ={0.5,0.0,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[10].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[10].SetCoord(coord);
//        NonPrism_topology[10] = node_id;
//        
//        
//        // 11th node of prism element is located at x ={0.5,0.5,0.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[11].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[11].SetCoord(coord);
//        NonPrism_topology[11] = node_id;
//        
//        // 12th node of prism element is located at x ={0.0,1.0,0.5}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[12].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[12].SetCoord(coord);
//        NonPrism_topology[12] = node_id;
//        
//        // 13th node of prism element is located at x ={0.0,0.0,0.5}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[13].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[13].SetCoord(coord);
//        NonPrism_topology[13] = node_id;
//        
//        // 14th node of prism element is located at x ={1.0,0.0,0.5}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = M_PI; // Z coordinate
//        geometry_Non3DPris.NodeVec()[14].SetNodeId(node_id);
//        geometry_Non3DPris.NodeVec()[14].SetCoord(coord);
//        NonPrism_topology[14] = node_id;
//        
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticPrism> (element_id, NonPrism_topology, physical_id, geometry_Non3DPris);
//    }
//    
//    geometry_Non3DPris.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non3DPris.txt");
//    geometry_Non3DPris.Print(file);
//    
//}
//
//
//// ------------------------------------------------------------------------------------------
//
//
//void NonThreeDHexaElements() {
//    
//    TPZGeoMesh geometry_Non3DHexa; // Create the objet that will describe the nonlinear geometry (3D Hexahedron).
//    int n_nodes = 20; // number of nodes
//    int geometry_dim = 3; // geometry dimension
//    std::string name("geometry Non3DHexa"); // geometry name
//    
//    // setting the object
//    geometry_Non3DHexa.SetName(name);
//    geometry_Non3DHexa.SetDimension(geometry_dim);
//    geometry_Non3DHexa.NodeVec().Resize(n_nodes);
//    
//    
//    int node_id = 0;
//    int element_id = 0;
//    int physical_id = 1;
//    TPZVec<long> NonHexahedron_topology(20);
//    TPZVec<REAL> coord(3,0.0);
//    
//    {
//        
//        // 1st node of hexahedron element is located at x ={-1.0,-1.0,-1.0}
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[0].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[0].SetCoord(coord);
//        NonHexahedron_topology[0] = node_id;
//        
//        // 2nd node of hexahedron element is located at x ={1.0,-1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[1].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[1].SetCoord(coord);
//        NonHexahedron_topology[1] = node_id;
//        
//        // 3rd node of hexahedron element is located at x ={1.0,1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[2].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[2].SetCoord(coord);
//        NonHexahedron_topology[2] = node_id;
//        
//        // 4th node of hexahedron element is located at x ={-1.0,1.0,-1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[3].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[3].SetCoord(coord);
//        NonHexahedron_topology[3] = node_id;
//        
//        // 5th node of hexahedron element is located at x ={-1.0,-1.0,1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[4].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[4].SetCoord(coord);
//        NonHexahedron_topology[4] = node_id;
//        
//        // 6th node of hexahedron element is located at x ={1.0,-1.0,1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[5].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[5].SetCoord(coord);
//        NonHexahedron_topology[5] = node_id;
//        
//        // 7th node of hexahedron element is located at x ={1.0,1.0,1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[6].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[6].SetCoord(coord);
//        NonHexahedron_topology[6] = node_id;
//        
//        // 8th node of hexahedron element is located at x ={-1.0,1.0,1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[7].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[7].SetCoord(coord);
//        NonHexahedron_topology[7] = node_id;
//        
//        // 9th node of hexahedron element is located at x ={0.0,-1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[8].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[8].SetCoord(coord);
//        NonHexahedron_topology[8] = node_id;
//        
//        
//        // 10th node of hexahedron element is located at x ={1.0,0.0,-1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[9].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[9].SetCoord(coord);
//        NonHexahedron_topology[9] = node_id;
//        
//        
//        // 11th node of hexahedron element is located at x ={0.0,1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[10].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[10].SetCoord(coord);
//        NonHexahedron_topology[10] = node_id;
//        
//        
//        // 12th node of hexahedron element is located at x ={-1.0,0.0,-1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[11].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[11].SetCoord(coord);
//        NonHexahedron_topology[11] = node_id;
//        
//        
//        // 13th node of hexahedron element is located at x ={0.0,1.0,-1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[12].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[12].SetCoord(coord);
//        NonHexahedron_topology[12] = node_id;
//        
//        
//        // 14th node of hexahedron element is located at x ={1.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[13].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[13].SetCoord(coord);
//        NonHexahedron_topology[13] = node_id;
//        
//        
//        // 15th node of hexahedron element is located at x ={0.0,1.0,1.0}
//        node_id++;
//        
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[14].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[14].SetCoord(coord);
//        NonHexahedron_topology[14] = node_id;
//        
//        
//        // 16th node of hexahedron element is located at x ={-1.0,0.0,1.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[15].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[15].SetCoord(coord);
//        NonHexahedron_topology[15] = node_id;
//        
//        
//        // 17th node of hexahedron element is located at x ={-1.0,-1.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[16].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[16].SetCoord(coord);
//        NonHexahedron_topology[16] = node_id;
//        
//        
//        // 18th node of hexahedron element is located at x ={1.0,-1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[17].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[17].SetCoord(coord);
//        NonHexahedron_topology[17] = node_id;
//        
//        // 19th node of hexahedron element is located at x ={1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[18].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[18].SetCoord(coord);
//        NonHexahedron_topology[18] = node_id;
//        
//        
//        // 20th node of hexahedron element is located at x ={-1.0,1.0,0.0}
//        node_id++;
//        
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        geometry_Non3DHexa.NodeVec()[19].SetNodeId(node_id);
//        geometry_Non3DHexa.NodeVec()[19].SetCoord(coord);
//        NonHexahedron_topology[19] = node_id;
//        
//        
//
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticCube> (element_id, NonHexahedron_topology, physical_id, geometry_Non3DHexa);
//    }
//    
//    geometry_Non3DHexa.BuildConnectivity();
//    
//    
//    std::ofstream file("geometry_Non3DHexa.txt");
//    geometry_Non3DHexa.Print(file);
//    
//}
//
//
//
//// ************************************** Create 3D nonlinear hexahedral  ***************************************
//
//
//TPZGeoMesh *CreateNonThreeDHexaGMesh(long nnodesthrehex, REAL Lx, REAL Ly, REAL Lz)
//{
//    TPZGeoMesh * gmesh_NonThreeDHexa = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 3; // geometry dimension
//    
//    std::string name("geomesh NonThreeDHexa"); // geometry name
//    gmesh_NonThreeDHexa->SetName(name);
//    gmesh_NonThreeDHexa->SetDimension(geometry_dim);
//    
//    
//    gmesh_NonThreeDHexa->NodeVec().Resize(nnodesthrehex); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodesthrehex);
//    
//    TPZVec<long> NonHexahedron_topology(20);
//    TPZVec<long> NonQuadrilateral_topology(8);
//
//    
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//    
//    long elementid = 0;
//    int physical_id = 1;
//    
//    // Index of boundary element
//    const int bc_front = -1; // define id for a material (border in front)
//    const int bc_right = -2; // define id for a material (border right)
//    const int bc_back = -3; // define id for a material (border back)
//    const int bc_left = -4; // define id for a material (border left)
//    const int bc_bottom = -5; // define id for a material (border bottom)
//    const int bc_top = -6; // define id for a material (border top)
//    
//    
//    // 0th element
//    
//    {
//        
//        // 0th node
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[0].SetNodeId(0);
//        gmesh_NonThreeDHexa->NodeVec()[0].SetCoord(coord);
//        NonHexahedron_topology[0] = 0;
//        
//        // 1st node
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[1].SetNodeId(1);
//        gmesh_NonThreeDHexa->NodeVec()[1].SetCoord(coord);
//        NonHexahedron_topology[1] = 1;
//        
//        // 2nd node
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[2].SetNodeId(2);
//        gmesh_NonThreeDHexa->NodeVec()[2].SetCoord(coord);
//        NonHexahedron_topology[2] = 2;
//        
//        // 3rd node
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = -1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[3].SetNodeId(3);
//        gmesh_NonThreeDHexa->NodeVec()[3].SetCoord(coord);
//        NonHexahedron_topology[3] = 3;
//        
//        // 4th node
//        coord[0] = -1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[4].SetNodeId(4);
//        gmesh_NonThreeDHexa->NodeVec()[4].SetCoord(coord);
//        NonHexahedron_topology[4] = 4;
//        
//        // 5th node
//        coord[0] = 1.0; // x coordinate
//        coord[1] = -1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[5].SetNodeId(5);
//        gmesh_NonThreeDHexa->NodeVec()[5].SetCoord(coord);
//        NonHexahedron_topology[5] = 5;
//        
//        // 6th node
//        coord[0] = 1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[6].SetNodeId(6);
//        gmesh_NonThreeDHexa->NodeVec()[6].SetCoord(coord);
//        NonHexahedron_topology[6] = 6;
//        
//        // 7th node
//        coord[0] = -1.0; // x coordinate
//        coord[1] = 1.0; // Y coordinate
//        coord[2] = 1.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[7].SetNodeId(7);
//        gmesh_NonThreeDHexa->NodeVec()[7].SetCoord(coord);
//        NonHexahedron_topology[7] = 7;
//        
//        // 8th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = -1.5; // Y coordinate
//        coord[2] = -1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[8].SetNodeId(8);
//        gmesh_NonThreeDHexa->NodeVec()[8].SetCoord(coord);
//        NonHexahedron_topology[8] = 8;
//        
//        
//        // 9th node
//        coord[0] = 1.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = -1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[9].SetNodeId(9);
//        gmesh_NonThreeDHexa->NodeVec()[9].SetCoord(coord);
//        NonHexahedron_topology[9] = 9;
//        
//        
//        // 10th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.5; // Y coordinate
//        coord[2] = -1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[10].SetNodeId(10);
//        gmesh_NonThreeDHexa->NodeVec()[10].SetCoord(coord);
//        NonHexahedron_topology[10] = 10;
//        
//        
//        // 11th node
//        coord[0] = -1.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = -1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[11].SetNodeId(11);
//        gmesh_NonThreeDHexa->NodeVec()[11].SetCoord(coord);
//        NonHexahedron_topology[11] = 11;
//        
//        
//        // 12th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.5; // Y coordinate
//        coord[2] = -1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[12].SetNodeId(12);
//        gmesh_NonThreeDHexa->NodeVec()[12].SetCoord(coord);
//        NonHexahedron_topology[12] = 12;
//        
//        
//        // 13th node
//        coord[0] = 1.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[13].SetNodeId(13);
//        gmesh_NonThreeDHexa->NodeVec()[13].SetCoord(coord);
//        NonHexahedron_topology[13] = 13;
//        
//        
//        // 14th node
//        coord[0] = 0.0; // x coordinate
//        coord[1] = 1.5; // Y coordinate
//        coord[2] = 1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[14].SetNodeId(14);
//        gmesh_NonThreeDHexa->NodeVec()[14].SetCoord(coord);
//        NonHexahedron_topology[14] = 14;
//        
//        
//        // 15th node
//        coord[0] = -1.5; // x coordinate
//        coord[1] = 0.0; // Y coordinate
//        coord[2] = 1.5; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[15].SetNodeId(15);
//        gmesh_NonThreeDHexa->NodeVec()[15].SetCoord(coord);
//        NonHexahedron_topology[15] = 15;
//        
//        
//        // 16th node
//        coord[0] = -1.5; // x coordinate
//        coord[1] = -1.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[16].SetNodeId(16);
//        gmesh_NonThreeDHexa->NodeVec()[16].SetCoord(coord);
//        NonHexahedron_topology[16] = 16;
//        
//        
//        // 17th node
//        coord[0] = 1.5; // x coordinate
//        coord[1] = -1.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[17].SetNodeId(17);
//        gmesh_NonThreeDHexa->NodeVec()[17].SetCoord(coord);
//        NonHexahedron_topology[17] = 17;
//        
//        // 18th node
//        coord[0] = 1.5; // x coordinate
//        coord[1] = 1.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[18].SetNodeId(18);
//        gmesh_NonThreeDHexa->NodeVec()[18].SetCoord(coord);
//        NonHexahedron_topology[18] = 18;
//        
//        
//        // 19th node
//        coord[0] = -1.5; // x coordinate
//        coord[1] = 1.5; // Y coordinate
//        coord[2] = 0.0; // Z coordinate
//        gmesh_NonThreeDHexa->NodeVec()[19].SetNodeId(19);
//        gmesh_NonThreeDHexa->NodeVec()[19].SetCoord(coord);
//        NonHexahedron_topology[19] = 19;
//        
//        
//        
//        
////        new TPZGeoElRefPattern< pzgeom::TPZQuadraticCube> (elementid, NonHexahedron_topology, physical_id, *gmesh_NonThreeDHexa);
//    }
//    
//    
//    // ********************* boundary ****************
//
//    
//    // in front
//
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 0;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 1;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 5;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 4;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 8;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 17;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 12;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 16;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_front, *gmesh_NonThreeDHexa);
//    }
//
//    // right
//    
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 1;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 2;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 6;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 5;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 9;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 18;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 13;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 17;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_right, *gmesh_NonThreeDHexa);
//    }
//    
//    // back
//    
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 3;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 2;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 6;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 7;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 10;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 18;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 14;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 19;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_back, *gmesh_NonThreeDHexa);
//    }
//    
//    // left
//    
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 0;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 3;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 7;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 4;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 11;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 19;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 15;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 16;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_left, *gmesh_NonThreeDHexa);
//    }
//    
//    
//    // bottom
//    
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 0;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 1;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 2;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 3;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 8;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 9;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 10;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 11;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_bottom, *gmesh_NonThreeDHexa);
//    }
//    
//    // top
//    
//    {
//        // 1st node
//        NonQuadrilateral_topology[0] = 4;
//        
//        // 2nd node
//        NonQuadrilateral_topology[1] = 5;
//        
//        // 3rd node
//        NonQuadrilateral_topology[2] = 6;
//        
//        // 4th node
//        NonQuadrilateral_topology[3] = 7;
//        
//        // 5th node
//        NonQuadrilateral_topology[4] = 12;
//        
//        // 6th node
//        NonQuadrilateral_topology[5] = 13;
//        
//        // 7th node
//        NonQuadrilateral_topology[6] = 14;
//        
//        // 8th node
//        NonQuadrilateral_topology[7] = 15;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad> (elementid, NonQuadrilateral_topology, bc_top, *gmesh_NonThreeDHexa);
//    }
//    
//
//    // Build the mesh
//    gmesh_NonThreeDHexa->BuildConnectivity();
//    
//    std::ofstream outgmeshNonThreeDHexa("geomesh_NonThreeDHexa.txt");
//    gmesh_NonThreeDHexa->Print(outgmeshNonThreeDHexa);
//    
//    std::ofstream vtkgmeshNonThreeDHexa("geomesh_NonThreeDHexa.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_NonThreeDHexa, vtkgmeshNonThreeDHexa);
//    
//    return gmesh_NonThreeDHexa;
//    
//}
//
//
//
//// ************************************** Create 3D arc  ***************************************
//
//TPZGeoMesh *CreateThreeDarcGMesh(long nnodesarc, REAL Rad)
//
//
//{
//    TPZGeoMesh * gmesh_ThreeDarc = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 2; // geometry dimension
//    REAL ModelRadius = Rad;
//    
//    std::string name("geomesh ThreeDarc"); // geometry name
//    gmesh_ThreeDarc->SetName(name);
//    gmesh_ThreeDarc->SetDimension(geometry_dim);
//    
//    
//    gmesh_ThreeDarc->NodeVec().Resize(nnodesarc); // Resize of the geometry mesh
//    TPZVec<TPZGeoNode> Node(nnodesarc);
//
//    
//    TPZVec<long> Arc_topology(3,0.0);
//    TPZVec<REAL> coord(3,0.0);
//    
//    // Index of element
//    
//    long elementid = 0;
//    
//    int arc1 = -1;
//    int arc2 = -2;
//    int arc3 = -3;
//    int arc4 = -4;
//
//    // Setting node coordantes for Arc
//    int id = 0;
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,0.0 );//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,0.0);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,0.0);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,0.0 );//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,0.0);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,0.0 );//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    gmesh_ThreeDarc->NodeVec()[id].SetNodeId(id);
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y
//    gmesh_ThreeDarc->NodeVec()[id].SetCoord(2,0.0);//coord Z
//    id++;
//    
//    // Create Geometrical Arc #1
//    // Definition of Arc coordenates
//    Arc_topology[0] = 1;
//    Arc_topology[1] = 2;
//    Arc_topology[2] = 5;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D> (elementid, Arc_topology, arc1, *gmesh_ThreeDarc);
//    elementid++;
//    
//    // Create Geometrical Arc #2
//    Arc_topology[0] = 2;
//    Arc_topology[1] = 3;
//    Arc_topology[2] = 6;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D> (elementid, Arc_topology, arc2, *gmesh_ThreeDarc);
//    elementid++;
//    
//    // Create Geometrical Arc #3
//    Arc_topology[0] = 3;
//    Arc_topology[1] = 4;
//    Arc_topology[2] = 7;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D> (elementid, Arc_topology, arc3, *gmesh_ThreeDarc);
//    elementid++;
//    
//    // Create Geometrical Arc #4
//    Arc_topology[0] = 4;
//    Arc_topology[1] = 1;
//    Arc_topology[2] = 8;
//    new TPZGeoElRefPattern< pzgeom::TPZArc3D> (elementid, Arc_topology, arc4, *gmesh_ThreeDarc);
//    elementid++;
//    
//    // Build the mesh
//    gmesh_ThreeDarc->BuildConnectivity();
//    
//    std::ofstream outgmeshThreeDarc("geomesh_ThreeDarc.txt");
//    gmesh_ThreeDarc->Print(outgmeshThreeDarc);
//    
//    std::ofstream vtkgmeshThreeDarc("geomesh_ThreeDarc.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDarc, vtkgmeshThreeDarc);
//    return gmesh_ThreeDarc;
//    
//}
//
//// ************************************** Create 3D elipse  ***************************************
//
//void ThreeEllipsMeshGenerate()
//{
//  
//    TPZGeoMesh * gmesh_ellipsthreed = new TPZGeoMesh;
//    long geometry_dim = 2; // geometry dimension
//    long num_nodes = 4; // number of nodes
//
//    std::string name("geomesh ellipsthreed"); // geometry name
//    
//    gmesh_ellipsthreed->SetName(name);
//    gmesh_ellipsthreed->SetDimension(geometry_dim);
//    gmesh_ellipsthreed->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
//
//
//    
//    TPZManVector<long,4> EllipsSur_topology(4);
//    TPZManVector<long,2> EllipsLin_topology(2);
//    
//    TPZManVector<REAL,3> coord(3,0.);
//
//    // Index of element
//    int physical_id = 1;
//    long elementid = 0;
//    
//    const REAL semiX = 2.0;
//    const REAL semiY = 1.0;
//    
//    
//    TPZManVector<REAL,3> ellipOrigin(3,0.);
//    TPZManVector<REAL,3> semiAxeX(3,0.);
//    TPZManVector<REAL,3> semiAxeY(3,0.);
//    
//    semiAxeX[0] = semiX;
//    semiAxeY[1] = semiY;
//
//    
//    // 0th element
//    {
//    
//    REAL innerRectangleLx = 2.*semiX/sqrt(2.);
//    REAL innerRectangleLy = 2.*semiY/sqrt(2.);
//
//    // 0th node
//    coord[0] = -innerRectangleLx/2.; // x coordinate
//    coord[1] = -innerRectangleLy/2.; // Y coordinate
//    coord[2] = 0.0; // Z coordinate
//    gmesh_ellipsthreed->NodeVec()[0].SetNodeId(0);
//    gmesh_ellipsthreed->NodeVec()[0].SetCoord(coord);
//    EllipsSur_topology[0] = 0;
//
//    // 1st node
//    coord[0] = +innerRectangleLx/2.; // x coordinate
//    coord[1] = -innerRectangleLy/2.; // Y coordinate
//    coord[2] = 0.0; // Z coordinate
//    gmesh_ellipsthreed->NodeVec()[1].SetNodeId(1);
//    gmesh_ellipsthreed->NodeVec()[1].SetCoord(coord);
//    EllipsSur_topology[1] = 1;
//    
//    // 2nd node
//    coord[0] = +innerRectangleLx/2.; // x coordinate
//    coord[1] = +innerRectangleLy/2.; // Y coordinate
//    coord[2] = 0.0; // Z coordinate
//    gmesh_ellipsthreed->NodeVec()[2].SetNodeId(2);
//    gmesh_ellipsthreed->NodeVec()[2].SetCoord(coord);
//    EllipsSur_topology[2] = 2;
//
//    // 3rd node
//    coord[0] = -innerRectangleLx/2.; // x coordinate
//    coord[1] = +innerRectangleLy/2.; // Y coordinate
//    coord[2] = 0.0; // Z coordinate
//    gmesh_ellipsthreed->NodeVec()[3].SetNodeId(3);
//    gmesh_ellipsthreed->NodeVec()[3].SetCoord(coord);
//    EllipsSur_topology[3] = 3;
//
//    new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(elementid, EllipsSur_topology, physical_id, *gmesh_ellipsthreed);
//
//    }
//    
//    {
//    // 1st element
//    EllipsLin_topology[0] = 0;
//    EllipsLin_topology[1] = 1;
//    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge4 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (elementid, EllipsLin_topology, physical_id, *gmesh_ellipsthreed);
//    ellipEdge4->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
//    }
//    
//    {
//    // 2nd element
//    EllipsLin_topology[0] = 1;
//    EllipsLin_topology[1] = 2;
//    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge5 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (elementid, EllipsLin_topology, physical_id, *gmesh_ellipsthreed);
//    ellipEdge5->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
//    }
//    
//    {
//    // 3rd element
//    EllipsLin_topology[0] = 2;
//    EllipsLin_topology[1] = 3;
//    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge6 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (elementid, EllipsLin_topology, physical_id, *gmesh_ellipsthreed);
//    ellipEdge6->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
//    }
//    
//    {
//    // 4th element
//    EllipsLin_topology[0] = 3;
//    EllipsLin_topology[1] = 0;
//    TPZGeoElRefPattern<pzgeom::TPZEllipse3D> * ellipEdge7 = new TPZGeoElRefPattern<pzgeom::TPZEllipse3D> (elementid, EllipsLin_topology, physical_id, *gmesh_ellipsthreed);
//    ellipEdge7->Geom().SetAxes(ellipOrigin,semiAxeX,semiAxeY);
//    }
//    
//    gmesh_ellipsthreed->BuildConnectivity();
//    
//    
//    std::ofstream outgmeshEllipsthree("geometry_ellipsthreed.txt");
//    gmesh_ellipsthreed->Print(outgmeshEllipsthree);
//    
//    
//    std::ofstream vtkgmeshEllipsthree("geomesh_ellipsthreed.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ellipsthreed, vtkgmeshEllipsthree);
//    
//}
//
//
//// ************************************* (Finding nodes and elements) ********************************************
//
//
//
//TPZGeoMesh *CreateOneDFindNoElGMesh(long num_elem, REAL size_elem)
//{
//    TPZGeoMesh * gmesh_OneDFindNoEl = new TPZGeoMesh; // Initilized of TPZGeoMesh class
//    
//    long geometry_dim = 1; // geometry dimension
//    std::string name("geomesh OneDFindNoEl"); // geometry name
//    
//    gmesh_OneDFindNoEl->SetName(name);
//    gmesh_OneDFindNoEl->SetDimension(geometry_dim);
//    
//    long num_nodes = num_elem + 1; // Number of the nodes
//    gmesh_OneDFindNoEl->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
//    
//    
//    const int physical_id = 1; // Define id for material
//    const int bc0 = -1; // Define id for left boundary condition
//    const int bc1 = -2; // Define id for right boundary condition
//    
//    for (long i = 0 ; i < num_nodes; i++)
//    {
//        const REAL valElem = i * size_elem;
//        TPZVec <REAL> coord(3,0.);
//        coord[0] = valElem;
//        coord[1] = 0.0;
//        coord[2] = 0.0;
//        
//        gmesh_OneDFindNoEl->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
//        gmesh_OneDFindNoEl->NodeVec()[i].SetNodeId(i); // The id identification
//    }
//    
//    // Creating linear element and  zero-dimensional boundary element
//    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
//    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
//    long element_id = 0;
//    
//    
//    // Elements
//    
//    for (long iel = 0; iel < num_elem; iel++)
//    {
//        const long inod_l = iel;
//        const long inod_r = iel + 1;
//        Linear_topology[0] = inod_l;
//        Linear_topology[1] = inod_r;
//        
//        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDFindNoEl);
//    }
//    
//    element_id++;
//    
//    
//    // Left boundary condition
//    point_topology[0] = 0;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDFindNoEl);
//    element_id++;
//    
//    
//    // Right boundary condition
//    point_topology[0] = num_nodes-1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDFindNoEl);
//    
//    
//    gmesh_OneDFindNoEl->BuildConnectivity(); // Construct mesh neighbor connectivity
//    
//    std::ofstream outgmeshOneDFindNoEl("geomesh_OneDFindNoEl.txt");
//    gmesh_OneDFindNoEl->Print(outgmeshOneDFindNoEl);
//    
//    std::ofstream vtkgmeshOneDFindNoEl("geomesh_OneDFindNoEl.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDFindNoEl, vtkgmeshOneDFindNoEl);
//    
//    return gmesh_OneDFindNoEl;
//}


// ----------------------------------------------------------------------------------------
