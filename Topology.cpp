//
//  Topology.cpp
//  GeometryDescription
//
//  Created by Manouchehr Sanei on 8/22/19.
//
//

#include "Topology.h"

#include <iostream>

#include "pzvec.h"
#include "tpzintpoints.h"
#include "pzquad.h"

#include "pztrnsform.h"
#include "pzshapequad.h"



/// Type of topology
#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"



/// check






void CheckTopology(){
    
    TPZVec<REAL> blob;
    blob.resize(3);
    
    blob[0] = 10.;
    std::cout << blob[0] << std::endl;
    
    /// Definition of number of corners, sides, dimension and faces.
   
    pztopology::TPZQuadrilateral quad;
    std::cout << "Corner Nodes" << std::endl;
    std::cout << quad.NCornerNodes << std::endl;
    std::cout << "Sides" << std::endl;
    std::cout << quad.NSides << std::endl;
    std::cout << "Dimension" << std::endl;
    std::cout << quad.Dimension << std::endl;
    std::cout << "Faces" << std::endl;
    std::cout << quad.NFaces << std::endl;
    
    
    pztopology::TPZTriangle triangl;
    
    /// Definition of the dimension of each side and associated corner nodes.
    
    int side = 2;
    int dimside = quad.SideDimension(side);
    std::cout << "The dimension of the side" << " " <<  side << " " << "is" << " " << dimside << std::endl;
    int nsidenodes = quad.NSideNodes(side);
    std::cout << "The number of associated corner nodes of the side " << " " <<  side << " " << " is" << " " << nsidenodes << std::endl;
    
    /// Definition of the parametric transformation between sides (note that this function returns a TPZTransform).
    
    TPZVec<REAL> vectorin(2,0.33);
    
    int sidefrom = 6;
    int sideto = 3;
    TPZTransform<> tr = triangl.SideToSideTransform(sidefrom, sideto);
    
    TPZFMatrix<REAL> mult = tr.Mult();
    TPZFMatrix<REAL> sum  = tr.Sum();
    mult.Print(std::cout);
    sum.Print(std::cout);
    
    TPZVec<double> vectorout(1);
    tr.Apply(vectorin, vectorout);
    std::cout<< vectorout[0]<<std::endl;
    
    
    
    /// Creation of integration rules associated to each side (note that the integration rule is a pointer to the TPZIntPoints type).
    
    int sideint = 6;
    
    int order = 3;
    TPZIntPoints * integr = triangl.CreateSideIntegrationRule(sideint, order);
    int npoint = integr->NPoints();
    
    TPZManVector<REAL,3> x;
    TPZManVector<REAL,3> par_space;
    REAL w;

//    for (int ip = 0; ip < npoint; ip++)
//    {
//        integr->Point(ip, par_space, w);
//        
//        std::cout << "weight = " << w << std::endl;
//        std::cout << "Parametric space = " << par_space[ip] << std::endl;
//    }
//    
    
    /// Definition of a transformation index associated with a side.
    int sid = 8;
    
    TPZVec<int64_t> id(4,0);
    id[0]=4;
    id[1]=3;
    id[2]=2;
    id[3]=1;
    
    int val = quad.GetTransformId(8, id);
    
    
    
    
    /// Definition of the permutation index associated with the element.
    
    TPZVec<int> permgather; // output variable
//    quad.GetSideHDivPermutation(transf_id, permgather);
    
    ///  Relationship between sides: which sides are included in the closure of a side and orientation (note that the output variable is a stack object).
    
    TPZStack<int> high;
    quad.HigherDimensionSides(side, high);
    high.Print();
    
    TPZStack<int> low;
    quad.LowerDimensionSides(side, low);
    low.Print();
    
    /// Projection of a point to a side (this function is implemented on TPZShape).
    
    
    int rib = 2;
    TPZVec<REAL> in;
    REAL out;
//    pzshape::TPZShapeQuad::ProjectPoint2dQuadToRib(rib, in, out);
    
    
    
    std::cout << std::endl;
    std::cout << "Finish" << std::endl;

};
