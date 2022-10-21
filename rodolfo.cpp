//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------FLUID-------------------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <cassert>

#include <metis.h>


#include "Element.hpp"
#include "Boundary.hpp"

template<int DIM>
class Fluid{

public:
    typedef Element<DIM> Elements;
    typedef typename Elements::Nodes  Node;
    typedef Boundary<DIM> Boundaries;
    std::vector<Node *>       nodes_;
    std::vector<Node *>       nodesp_;
    std::vector<Elements *>   elements_;
    std::vector<Boundaries *> boundary_;
    typedef typename std::vector<Node *>::iterator     NodeIt;

private:
    std::string inputFile;
    int numElem;
    int numNodes;
    int numNodesP;
    int numTimeSteps;
    int printFreq;    
    int numBoundaries;
    int numBoundElems;
    double velocityInf[3];
    double pressInf;
    double rhoInf;
    double tempInf;
    double viscInf;
    double ktermInf;
    double dTime;
    double fieldForces[3];
    idx_t* part_elem;
    idx_t* part_nodes;

public:
    void dataReading(std::string inputFile);
    void buildLinearMesh();   
    idx_t* domainDecompositionMETIS(); 
    void printVelocity(int step);
    void printPressure(int step);
    int solveSteadyProblem(int iterNumber, double tolerance, int problem_type);
    int solveTransientProblem(int iterNumber, \
                              double tolerance, int problem_type);

    int getNumberElements(){return numElem;}
    int getNumberNodes(){return numNodes;}
    int getNumberNodesP(){return numNodesP;}
    int getNumberTimeSteps(){return numTimeSteps;}
    int getPrintingFrequence(){return printFreq;}
    double* getUndistVelocity(){return velocityInf;}
    double getUndistPressure(){return pressInf;}
    double getUndistDensity(){return rhoInf;}
    double getUndistViscosity(){return viscInf;}
    double getUndistTemperature(){return tempInf;}
    double getUndistTCondutivity(){return ktermInf;}
    double getTimeStep(){return dTime;}
    double* getFieldForces(){return fieldForces;}
    
    std::vector<Node *> getNodesVelocity(){return nodes_;}
    std::vector<Elements *> getElements(){return elements_;}
    std::vector<Node *> getNodesPressure(){return nodesp_;}

    NodeIt nodesBeginIt() {return nodes_.begin();}
    NodeIt nodesEndIt() {return nodes_.end();}

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-------------------------------BUILD LINEAR MESH------------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::buildLinearMesh() {
    
    //Compute the number of nodes in the linear mesh
    int connect_aux[numElem][3];
    for (int jel = 0; jel < numElem; jel++){
        
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();

        connect_aux[jel][0] = connec(0);
        connect_aux[jel][1] = connec(1);
        connect_aux[jel][2] = connec(2);
    };

    numNodesP=3;

    for (int jel = 1; jel < numElem; jel++){
        
        int iflag1 = 0;
        int iflag2 = 0;
        int iflag3 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][0] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][0] == connect_aux[jel2][1]) ||  \
                (connect_aux[jel][0] == connect_aux[jel2][2])) {
                iflag1 += 1;
            };
        };
        (iflag1 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][1] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][1] == connect_aux[jel2][1]) ||  \
                (connect_aux[jel][1] == connect_aux[jel2][2])) {
                iflag2 += 1;
            };
        };
        (iflag2 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][2] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][2] == connect_aux[jel2][1]) ||  \
                (connect_aux[jel][2] == connect_aux[jel2][2])) {
                iflag3 += 1;
            };
        };
        (iflag3 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

    };

    //Creates the linear mesh element connectivity
    typename Elements::ConnectivityP connect;

    connect(0) = 0;
    connect(1) = 1;
    connect(2) = 2;

    elements_[0]->setConnectivityP(connect);
    
    int iflag1 = 2;

    int iflag2 = 0;

    for (int jel = 1; jel < numElem; jel++){
        
        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 3; j++){
                if (connect_aux[jel][0] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(0) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(0) = iflag1;
        };
        
        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 3; j++){
                if (connect_aux[jel][1] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(1) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(1) = iflag1;
        };

        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 3; j++){
                if (connect_aux[jel][2] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(2) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(2) = iflag1;
        };

        elements_[jel]->setConnectivityP(connect);

    };

    //Sets the nodal coordinates in the linear mesh
    iflag1 = -1;
    int index = 0;
    for (int jel = 0; jel<numElem; jel++){
        typename Node::VecLocD x;
        typename Elements::Connectivity connectv;
        typename Elements::ConnectivityP connectp;
        
        connectv = elements_[jel] -> getConnectivity();
        connectp = elements_[jel] -> getConnectivityP();

        for (int i=0; i<3; i++){
            if(connectp(i) > iflag1){
                x = nodes_[connectv(i)] -> getCoordinates();
                
                Node *node = new Node(x, index++);
                nodesp_.push_back(node);
                iflag1 += 1;
            };
        };
    };

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::cout << " NodesV = " << numNodes \
              << " NodesP = " << numNodesP << std::endl;
    
    //Send linear mesh information to the elements
    for (int jel = 0; jel<numElem; jel++){
        elements_[jel] -> setNodesP(nodesp_);
    };    

};

template<>
void Fluid<3>::buildLinearMesh() {
    
    //Compute the number of nodes in the linear mesh
    int connect_aux[numElem][4];
    for (int jel = 0; jel < numElem; jel++){
        
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();

        connect_aux[jel][0] = connec(0);
        connect_aux[jel][1] = connec(1);
        connect_aux[jel][2] = connec(2);
        connect_aux[jel][3] = connec(3);
    };

    numNodesP=4;

    for (int jel = 1; jel < numElem; jel++){
        
        int iflag1 = 0;
        int iflag2 = 0;
        int iflag3 = 0;
        int iflag4 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][0] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][0] == connect_aux[jel2][1]) ||  \
                (connect_aux[jel][0] == connect_aux[jel2][2]) ||  \
                (connect_aux[jel][0] == connect_aux[jel2][3])) {
                iflag1 += 1;
            };
        };
        (iflag1 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][1] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][1] == connect_aux[jel2][1]) ||  \
                (connect_aux[jel][1] == connect_aux[jel2][2]) ||  \
                (connect_aux[jel][1] == connect_aux[jel2][3])) {
                iflag2 += 1;
            };
        };
        (iflag2 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][2] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][2] == connect_aux[jel2][1]) ||  
                (connect_aux[jel][2] == connect_aux[jel2][2]) ||  \
                (connect_aux[jel][2] == connect_aux[jel2][3])) {
                iflag3 += 1;
            };
        };
        (iflag3 == 0 ? numNodesP += 1 : numNodesP = numNodesP);

        for (int jel2 = 0; jel2 < jel; jel2++){
            if ((connect_aux[jel][3] == connect_aux[jel2][0]) ||  \
                (connect_aux[jel][3] == connect_aux[jel2][1]) ||  
                (connect_aux[jel][3] == connect_aux[jel2][2]) ||  \
                (connect_aux[jel][3] == connect_aux[jel2][3])) {
                iflag4 += 1;
            };
        };
        (iflag4 == 0 ? numNodesP += 1 : numNodesP = numNodesP);
    };

    //Creates the element connectivity
    typename Elements::ConnectivityP connect;

    connect(0) = 0;
    connect(1) = 1;
    connect(2) = 2;
    connect(3) = 3;

    elements_[0]->setConnectivityP(connect);
    
    int iflag1 = 3;

    int iflag2 = 0;

    for (int jel = 1; jel < numElem; jel++){
        
        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 4; j++){
                if (connect_aux[jel][0] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(0) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(0) = iflag1;
        };
        
        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 4; j++){
                if (connect_aux[jel][1] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(1) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(1) = iflag1;
        };

        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 4; j++){
                if (connect_aux[jel][2] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(2) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(2) = iflag1;
        };

        iflag2 = 0;

        for (int jel2 = 0; jel2 < jel; jel2++){
            for (int j = 0; j < 4; j++){
                if (connect_aux[jel][3] == connect_aux[jel2][j]) { 
                    typename Elements::ConnectivityP connectp;
                    connectp = elements_[jel2] -> getConnectivityP();
                    connect(3) = connectp(j);
                    iflag2 += 1;
                };
            };
        };
        if (iflag2 == 0){
            iflag1 += 1;
            connect(3) = iflag1;
        };

        elements_[jel]->setConnectivityP(connect);

    };

    //Sets the nodal coordinates in the linear mesh
    iflag1 = -1;
    int index = 0;
    for (int jel = 0; jel<numElem; jel++){
        typename Node::VecLocD x;
        typename Elements::Connectivity connectv;
        typename Elements::ConnectivityP connectp;
        
        connectv = elements_[jel] -> getConnectivity();
        connectp = elements_[jel] -> getConnectivityP();

        for (int i=0; i<4; i++){
            if(connectp(i) > iflag1){
                x = nodes_[connectv(i)] -> getCoordinates();
                
                Node *node = new Node(x, index++);
                nodesp_.push_back(node);
                iflag1 += 1;
            };
        };              
    };

    std::cout << " NodesV = " << numNodes \
              << " NodesP = " << numNodesP << std::endl;

    //Send linear mesh information to the elements
    for (int jel = 0; jel<numElem; jel++){
        elements_[jel] -> setNodesP(nodesp_);
    };    
};

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<>
idx_t* Fluid<2>::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[(4*dd-2)*numEl];
    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];


    for (idx_t i = 0; i < numEl+1; i++){
        elem_start[i]=(4*dd-2)*i;
    };
    for (idx_t jel = 0; jel < numEl; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();        
        
        for (idx_t i=0; i<(4*dd-2); i++){
        elem_connec[(4*dd-2)*jel+i] = connec(i);
        };
    };

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                              NULL, NULL, &one, &ssize, NULL, NULL,    \
                              &objval, part_elem, part_nodes);

    mirrorData << std::endl << "Divisao dos elementos entre os processos" \
               << std::endl;
    for(int i = 0; i < numElem; i++){
        mirrorData << "part elem = " << part_elem[i] << " " \
                   << i << " " << std::endl;
    };

    mirrorData << std::endl << "Divisao dos nos entre os processos" \
               << std::endl;
    for(int i = 0; i < numNodes; i++){
        mirrorData << "part nodes = " << part_nodes[i] << " " \
                   << i << " " << std::endl;
    };
    
    return part_elem;

};

template<>
idx_t* Fluid<3>::domainDecompositionMETIS() {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numElem;
    idx_t numNd = numNodes;
    idx_t dd = 3;
    idx_t ssize = size;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[(4*dd-2)*numEl];
    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];

    for (idx_t i = 0; i < numEl+1; i++){
        elem_start[i]=(4*dd-2)*i;
    };
    for (idx_t jel = 0; jel < numEl; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();        
        
        for (idx_t i=0; i<(4*dd-2); i++){
        elem_connec[(4*dd-2)*jel+i] = connec(i);
        };
    };

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                              NULL, NULL, &one, &ssize, NULL, NULL,    \
                              &objval, part_elem, part_nodes);

    mirrorData << std::endl << "Divisao dos elementos entre os processos" \
               << std::endl;
    for(int i = 0; i < numElem; i++){
        mirrorData << "part elem = " << part_elem[i] << " " \
                   << i << " " << std::endl;
    };

    mirrorData << std::endl << "Divisao dos nos entre os processos" \
               << std::endl;
    for(int i = 0; i < numNodes; i++){
        mirrorData << "part nodes = " << part_nodes[i] << " " \
                   << i << " " << std::endl;
    };

    return part_elem;
    
};

//------------------------------------------------------------------------------
//----------------------------PRINT VELOCITY RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::printVelocity(int step) {

    //    std::cout << "Printing Velocity Results" << std::endl;
    
    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVel"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std :: endl
           << "<VTKFile type=\"UnstructuredGrid\">" << std :: endl
           << "  <UnstructuredGrid>" << std :: endl
           << "  <Piece NumberOfPoints=\"" << numNodes
           << "\"  NumberOfCells=\"" << numElem
           << "\">" << std :: endl;

    // write nodal coordinates
    output_v << "    <Points>" << std :: endl
           << "      <DataArray type=\"Float64\" "
           << "NumberOfComponents=\"3\" format=\"ascii\">" << std :: endl;

    for (int i=0; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();
        output_v << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_v << "      </DataArray>" << std :: endl
           << "    </Points>" << std :: endl;

    // write element connectivity
    output_v << "    <Cells>" << std :: endl
           << "      <DataArray type=\"Int32\" "
           << "Name=\"connectivity\" format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        typename Elements::Connectivity connec;
        connec=elements_[i]->getConnectivity();
        
        output_v << connec(0) << " " << connec(1) << " " << connec(2) << " "  \
               << connec(3) << " " << connec(4) << " " << connec(5) <<  \
            std::endl;
    };
    output_v << "      </DataArray>" << std :: endl;
  
    // write offsets in the data array
    output_v << "      <DataArray type=\"Int32\""
           << " Name=\"offsets\" format=\"ascii\">" << std :: endl;

    int aux = 0;
    for (int i=0; i<numElem; i++){
        output_v << aux + 6 << std::endl;
        aux += 6;
    };
    output_v << "      </DataArray>" << std :: endl;
  
    // write element types
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
           << "format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        output_v << 22 << std::endl;
    };

    output_v << "      </DataArray>" << std :: endl
           << "    </Cells>" << std :: endl;


    //Write point results
    output_v << "    <PointData>" << std :: endl;
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 << "Name=\"Velocity\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numNodes; i++){
        output_v << nodes_[i] -> getVelocity(0) << " "          \
                 << nodes_[i] -> getVelocity(1) << " " << 0.0 << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

   //  output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 // << "Name=\"Stress\" format=\"ascii\">" << std :: endl;
   //  for (int i=0; i<Nnos; i++){
   //      output_v << 0.0 << " " << 0.0 << " " << 0.0 << endl;
   //  };
    
   //  output_v << "      </DataArray> " << std::endl;

    output_v << "    </PointData>" << std :: endl; 

    //Write cell results
    output_v << "    <CellData>" << std :: endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	 << "Name=\"Process\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numElem; i++){
        output_v << part_elem[i] << std::endl;
    };
    
    output_v << "      </DataArray> " << std::endl;
    output_v << "    </CellData>" << std :: endl; 



    // finalise
    output_v << "  </Piece>" << std :: endl
           << "  </UnstructuredGrid>" << std :: endl
           << "</VTKFile>" << std :: endl;

};

template<>
void Fluid<3>::printVelocity(int step) {

    //    std::cout << "Printing Velocity Results" << std::endl;

    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaVel"+result+".vtu";
    
    std::fstream output_v(s.c_str(), std::ios_base::out);

    output_v << "<?xml version=\"1.0\"?>" << std :: endl
           << "<VTKFile type=\"UnstructuredGrid\">" << std :: endl
           << "  <UnstructuredGrid>" << std :: endl
           << "  <Piece NumberOfPoints=\"" << numNodes
           << "\"  NumberOfCells=\"" << numElem
           << "\">" << std :: endl;

    // write nodal coordinates
    output_v << "    <Points>" << std :: endl
           << "      <DataArray type=\"Float64\" "
           << "NumberOfComponents=\"3\" format=\"ascii\">" << std :: endl;

    for (int i=0; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();
        output_v << x(0) << " " << x(1) << " " << x(2) << std::endl;        
    };
    output_v << "      </DataArray>" << std :: endl
           << "    </Points>" << std :: endl;

    // write element connectivity
    output_v << "    <Cells>" << std :: endl
           << "      <DataArray type=\"Int32\" "
           << "Name=\"connectivity\" format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        typename Elements::Connectivity connec;
        connec=elements_[i]->getConnectivity();
        
        output_v << connec(0) << " " << connec(1) << " " << connec(2) << " "  \
                 << connec(3) << " " << connec(4) << " " << connec(7) << " "  \
                 << connec(5) << " " << connec(6) << " " << connec(9) << " "  \
                 << connec(8) << std::endl;
    };
    output_v << "      </DataArray>" << std :: endl;
  
    // write offsets in the data array
    output_v << "      <DataArray type=\"Int32\""
           << " Name=\"offsets\" format=\"ascii\">" << std :: endl;

    int aux = 0;
    for (int i=0; i<numElem; i++){
        output_v << aux + 10 << std::endl;
        aux += 10;
    };
    output_v << "      </DataArray>" << std :: endl;
  
    // write element types
    output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
           << "format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        output_v << 24 << std::endl;
    };

    output_v << "      </DataArray>" << std :: endl
           << "    </Cells>" << std :: endl;


    //Write point results
    output_v << "    <PointData>" << std :: endl;
    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 << "Name=\"Velocity\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numNodes; i++){
        output_v << nodes_[i] -> getVelocity(0) << " "          \
                 << nodes_[i] -> getVelocity(1) << " "          \
                 << nodes_[i] -> getVelocity(2) << std::endl;
    };
    output_v << "      </DataArray> " << std::endl;

   //  output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 // << "Name=\"Stress\" format=\"ascii\">" << std :: endl;
   //  for (int i=0; i<Nnos; i++){
   //      output_v << 0.0 << " " << 0.0 << " " << 0.0 << endl;
   //  };
    
   //  output_v << "      </DataArray> " << std::endl;

    output_v << "    </PointData>" << std :: endl; 

    //Write cell results
    output_v << "    <CellData>" << std :: endl;

    output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	 << "Name=\"Process\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numElem; i++){
        output_v << part_elem[i] << std::endl;
    };
    
    output_v << "      </DataArray> " << std::endl;
    output_v << "    </CellData>" << std :: endl; 



    // finalise
    output_v << "  </Piece>" << std :: endl
           << "  </UnstructuredGrid>" << std :: endl
           << "</VTKFile>" << std :: endl;

};

//------------------------------------------------------------------------------
//----------------------------PRINT PRESSURE RESULTS----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::printPressure(int step) {

    //    std::cout << "Printing Pressure Results" << std::endl;

    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaPres"+result+".vtu";
    
    std::fstream output_p(s.c_str(), std::ios_base::out);

    int aux;

    output_p << "<?xml version=\"1.0\"?>" << std :: endl
           << "<VTKFile type=\"UnstructuredGrid\">" << std :: endl
           << "  <UnstructuredGrid>" << std :: endl
           << "  <Piece NumberOfPoints=\"" << numNodesP
           << "\"  NumberOfCells=\"" << numElem
           << "\">" << std :: endl;

    // write nodal coordinates
    output_p << "    <Points>" << std :: endl
           << "      <DataArray type=\"Float64\" "
           << "NumberOfComponents=\"3\" format=\"ascii\">" << std :: endl;

    for (int i=0; i<numNodesP; i++){
        typename Node::VecLocD x;
        x=nodesp_[i]->getCoordinates();
        output_p << x(0) << " " << x(1) << " " << 0.0 << std::endl;        
    };
    output_p << "      </DataArray>" << std :: endl
           << "    </Points>" << std :: endl;

    // write element connectivity
    output_p << "    <Cells>" << std :: endl
           << "      <DataArray type=\"Int32\" "
           << "Name=\"connectivity\" format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        typename Elements::ConnectivityP connec;
        connec=elements_[i]->getConnectivityP();
        
        output_p << connec(0) << " " << connec(1) << " " \
                        << connec(2) << std::endl;
    };
    output_p << "      </DataArray>" << std :: endl;
  
    // write offsets in the data array
    output_p << "      <DataArray type=\"Int32\""
           << " Name=\"offsets\" format=\"ascii\">" << std :: endl;

    aux = 0;
    for (int i=0; i<numElem; i++){
        output_p << aux + 3 << std::endl;
        aux += 3;
    };
    output_p << "      </DataArray>" << std :: endl;
  
    // write element types
    output_p << "      <DataArray type=\"UInt8\" Name=\"types\" "
           << "format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        output_p << 5 << std::endl;
    };

    output_p << "      </DataArray>" << std :: endl
           << "    </Cells>" << std :: endl;


    //Write point results
    output_p << "    <PointData>" << std :: endl;
    output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	 << "Name=\"Pressure\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numNodesP; i++){
        output_p << nodesp_[i] -> getPressure() << std::endl;
    };
    output_p << "      </DataArray> " << std::endl;

   //  output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 // << "Name=\"Stress\" format=\"ascii\">" << std :: endl;
   //  for (int i=0; i<Nnos; i++){
   //      output_p << 0.0 << " " << 0.0 << " " << 0.0 << endl;
   //  };
    
   //  output_p << "      </DataArray> " << std::endl;

    output_p << "    </PointData>" << std :: endl; 

    //Write cell results
    output_p << "    <CellData>" << std :: endl;

    output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Process\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numElem; i++){
        output_p << part_elem[i] << std::endl;
    };
    output_p << "      </DataArray> " << std::endl;
    output_p << "    </CellData>" << std :: endl; 



    // finalise
    output_p << "  </Piece>" << std :: endl
           << "  </UnstructuredGrid>" << std :: endl
           << "</VTKFile>" << std :: endl;
};

template<>
void Fluid<3>::printPressure(int step) {

    //    std::cout << "Printing Pressure Results" << std::endl;

    std::string result;
    std::ostringstream convert;

    convert << step+100000;
    result = convert.str();
    std::string s = "saidaPres"+result+".vtu";
    
    std::fstream output_p(s.c_str(), std::ios_base::out);

    int aux;

    output_p << "<?xml version=\"1.0\"?>" << std :: endl
           << "<VTKFile type=\"UnstructuredGrid\">" << std :: endl
           << "  <UnstructuredGrid>" << std :: endl
           << "  <Piece NumberOfPoints=\"" << numNodesP
           << "\"  NumberOfCells=\"" << numElem
           << "\">" << std :: endl;

    // write nodal coordinates
    output_p << "    <Points>" << std :: endl
           << "      <DataArray type=\"Float64\" "
           << "NumberOfComponents=\"3\" format=\"ascii\">" << std :: endl;

    for (int i=0; i<numNodesP; i++){
        typename Node::VecLocD x;
        x=nodesp_[i]->getCoordinates();
        output_p << x(0) << " " << x(1) << " " << x(2) << std::endl;        
    };
    output_p << "      </DataArray>" << std :: endl
           << "    </Points>" << std :: endl;

    // write element connectivity
    output_p << "    <Cells>" << std :: endl
           << "      <DataArray type=\"Int32\" "
           << "Name=\"connectivity\" format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        typename Elements::ConnectivityP connec;
        connec=elements_[i]->getConnectivityP();
        
        output_p << connec(0) << " " << connec(1) << " "  \
                 << connec(2) << " " << connec(3) << std::endl;
    };
    output_p << "      </DataArray>" << std :: endl;
  
    // write offsets in the data array
    output_p << "      <DataArray type=\"Int32\""
           << " Name=\"offsets\" format=\"ascii\">" << std :: endl;

    aux = 0;
    for (int i=0; i<numElem; i++){
        output_p << aux + 4 << std::endl;
        aux += 4;
    };
    output_p << "      </DataArray>" << std :: endl;
  
    // write element types
    output_p << "      <DataArray type=\"UInt8\" Name=\"types\" "
           << "format=\"ascii\">" << std :: endl;
    
    for (int i=0; i<numElem; i++){
        output_p << 10 << std::endl;
    };

    output_p << "      </DataArray>" << std :: endl
           << "    </Cells>" << std :: endl;

    //Write point results
    output_p << "    <PointData>" << std :: endl;
    output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
             << "Name=\"Pressure\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numNodesP; i++){
        output_p << nodesp_[i] -> getPressure() << std::endl;
    };
    output_p << "      </DataArray> " << std::endl;

   //  output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	 // << "Name=\"Stress\" format=\"ascii\">" << std :: endl;
   //  for (int i=0; i<Nnos; i++){
   //      output_p << 0.0 << " " << 0.0 << " " << 0.0 << endl;
   //  };
    
   //  output_p << "      </DataArray> " << std::endl;

    output_p << "    </PointData>" << std :: endl; 

    //Write cell results
    output_p << "    <CellData>" << std :: endl;

    output_p << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	 << "Name=\"Process\" format=\"ascii\">" << std :: endl;
    for (int i=0; i<numElem; i++){
        output_p << part_elem[i] << std::endl;
    };
    output_p << "      </DataArray> " << std::endl;
    output_p << "    </CellData>" << std :: endl; 

    // finalise
    output_p << "  </Piece>" << std :: endl
           << "  </UnstructuredGrid>" << std :: endl
           << "</VTKFile>" << std :: endl;

};

//------------------------------------------------------------------------------
//----------------------------READS FLUID INPUT FILE----------------------------
//------------------------------------------------------------------------------
template<>
void Fluid<2>::dataReading(std::string inputFile) {

    std::cout << "Reading data.." << std::endl;

    std::string mirror;
    std::string line;   
    mirror = "mirror.txt";
    
    //Defines input and output files    
    std::ifstream inputData(inputFile.c_str());
    std::ofstream mirrorData(mirror.c_str());

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numElem >> numNodes >> numTimeSteps >> printFreq;
    mirrorData << "Number of Elements     = " << numElem << std::endl;
    mirrorData << "Number of Nodes        = " << numNodes << std::endl;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density, temperature, viscosity and thermal condutivity
    inputData >> rhoInf >> tempInf >> viscInf >> ktermInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Temperature= " << tempInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;
    mirrorData << "Undisturbed Term. Cond.= " << ktermInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime;

    mirrorData << "Time Step              = " << dTime << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    int dimension=2;

    int index = 0;

    //Read nodal coordinates
    nodes_.reserve(numNodes);

    for (int inode=0; inode<numNodes; inode++){

        ublas::bounded_vector<double,2> x;

        for (int i=0; i<dimension; i++){
                inputData >> x(i);
        };
        Node *node = new Node(x, index++);
        nodes_.push_back(node);
    };

    mirrorData << "Nodal Coordinates" << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();       
        for (int j=0; j<dimension; j++){
            mirrorData << x(j) << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
    };

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read element connectivity
    elements_.reserve(numElem);
    index = 0;
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connect;

        for (int i=0; i<4*dimension-2; i++){
            inputData >> connect(i);
            connect(i)=connect(i)-1;
        };
        Elements *el = new Elements(index++,connect,nodes_);
        elements_.push_back(el);
    };
    
        
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();       
        for (int i=0; i<4*dimension-2; i++){
            mirrorData << connec(i) << " ";
        };
        mirrorData << std::endl;
    };

    //Sets Viscosity, density, time step, time integration 
    //scheme and field forces values
    for (int i=0; i<numElem; i++){
        elements_[i] -> setViscosity(viscInf);
        elements_[i] -> setDensity(rhoInf);
        elements_[i] -> setTimeStep(dTime);
        elements_[i] -> setTimeIntegrationScheme(1.);
        elements_[i] -> setFieldForce(fieldForces);
    };
    
    //Creates the linear mesh
    buildLinearMesh();
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);
    
    //Read number of boundary groups
    inputData >> numBoundaries >> numBoundElems;
    mirrorData << "Number of Boundary groups = " << numBoundaries << std::endl;
    mirrorData << "Number of Boundary Elements = " << numBoundElems <<std::endl;

    boundary_.reserve(numBoundElems);

    //Read boundary information: connectivity and prescribed value
    for (int ibound=0; ibound<numBoundaries; ibound++){

        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
        bool constrain[3];
        double value[3];
        int numGroup;
        Boundaries::BoundConnect connectB;
        index = 0;
        
        inputData >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] \
                  >> value[0] >> value[1] >> value[2];
        
        mirrorData << "Constrains = " << constrain[0] << " " << constrain[1] \
                   << " " << constrain[2] << std::endl;
        mirrorData << "Values = " << value[0] << " " << value[1]  \
                   << " " << value[2] << std::endl;
        
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);

        mirrorData << "Connectivity = " << std::endl;
        
        for (int i=0; i<numGroup; i++){
            inputData >> connectB(0) >> connectB(1) >> connectB(2);
            
            connectB(0) -= 1;
            connectB(1) -= 1;
            connectB(2) -= 1;

            Boundaries *bound = new Boundaries(connectB, index++,  \
                                               constrain, value);
            boundary_.push_back(bound);

            mirrorData << connectB(0) << " " << connectB(1) \
                       << " " << connectB(2) << std::endl;
        };
       
    };

    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        Boundaries::BoundConnect connectB;
        connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB(0);
        int no2 = connectB(1);
        int no3 = connectB(2);
        
        
        
        
        (boundary_[ibound] -> getConstrain(0)){
            nodes_[no1] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no2] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no3] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
        };
        if (boundary_[ibound] -> getConstrain(1)){
            nodes_[no1] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no2] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no3] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
        };               
    };
    
    //Sets nodal constrains
    for (int i=0; i<numNodes; i++){
        mirrorData<< "Constrains " << i \
                  << " " << nodes_[i]->getConstrains(0) \
                  << " " << nodes_[i]->getConstrainValue(0) \
                  << " " << nodes_[i]->getConstrains(1) \
                  << " " << nodes_[i]->getConstrainValue(1) << std::endl;

    }; 
};

template<>
void Fluid<3>::dataReading(std::string inputFile) {

    std::cout << "Reading data.." << std::endl;

    std::string mirror;
    std::string line;
    mirror = "mirror.txt";
    
    //Defines input and output files    
    std::ifstream inputData(inputFile.c_str());
    std::ofstream mirrorData(mirror.c_str());

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read number of nodes, elements, time steps and printing frequence
    inputData >> numElem >> numNodes >> numTimeSteps >> printFreq;
    mirrorData << "Number of Elements     = " << numElem << std::endl;
    mirrorData << "Number of Nodes        = " << numNodes << std::endl;
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;

    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read undisturbed density, temperature, viscosity and thermal condutivity
    inputData >> rhoInf >> tempInf >> viscInf >> ktermInf;

    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Temperature= " << tempInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;
    mirrorData << "Undisturbed Term. Cond.= " << ktermInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read time step lenght
    inputData >> dTime;

    mirrorData << "Time Step              = " << dTime << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];

    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
               << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    int dimension = 3;

    int index = 0;

    //Read nodal coordinates
    nodes_.reserve(numNodes);

    for (int inode=0; inode<numNodes; inode++){
        typename Node::VecLocD x;
        for (int i=0; i<dimension; i++){
                inputData >> x(i);
        };
        Node *node = new Node(x, index++);
        nodes_.push_back(node);
    };
   
    mirrorData << "Nodal Coordinates" << std::endl;
    for (int i = 0 ; i<numNodes; i++){
        typename Node::VecLocD x;
        x=nodes_[i]->getCoordinates();
        for (int j=0; j<dimension; j++){
            mirrorData << x(j) << " ";
        };
        mirrorData << std::endl;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
    };

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read element connectivity
    elements_.reserve(numElem);
    index = 0;
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connect;
        
        // inputData >> connect(0) >> connect(1) >> connect(2) >> connect(3) >> connect(4) >> connect(7) >> connect(5) >> connect(6) >> connect(8) >> connect(9);
        inputData >> connect(3) >> connect(0) >> connect(1) >> connect(2) \
                  >> connect(6) >> connect(4) >> connect(9) >> connect(8) \
                  >> connect(7) >> connect(5);
 
       for (int i=0; i<4*dimension-2; i++){
            connect(i) -= 1;
        };

        Elements *el = new Elements(index++,connect,nodes_);
        elements_.push_back(el);
    };
        
    mirrorData << std::endl << "Element Connectivity" << std::endl;        
    for (int jel=0; jel<numElem; jel++){
        typename Elements::Connectivity connec;
        connec=elements_[jel]->getConnectivity();
        for (int i=0; i<4*dimension-2; i++){
            mirrorData << connec(i) << " ";
        };
        mirrorData << std::endl;
    };

    //Sets Viscosity, density, time step, time integration 
    //scheme and field forces values
    for (int i=0; i<numElem; i++){
        elements_[i] -> setViscosity(viscInf);
        elements_[i] -> setDensity(rhoInf);
        elements_[i] -> setTimeStep(dTime);
        elements_[i] -> setTimeIntegrationScheme(1.);
        elements_[i] -> setFieldForce(fieldForces);
    };

    //Creates the linear mesh
    buildLinearMesh();

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);
    
    //Read number of boundary groups
    inputData >> numBoundaries >> numBoundElems;
    mirrorData << "Number of Boundary groups = " << numBoundaries << std::endl;
    mirrorData << "Number of Boundary Elements = " << numBoundElems <<std::endl;

    boundary_.reserve(numBoundElems);

    for (int ibound=0; ibound<numBoundaries; ibound++){

        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
        bool constrain[3];
        double value[3];
        int numGroup;
        Boundaries::BoundConnect connectB;
        index = 0;
        
        inputData >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] \
                  >> value[0] >> value[1] >> value[2];
        
        mirrorData << "Constrains = " << constrain[0] << " " << constrain[1] \
                   << " " << constrain[2] << std::endl;
        mirrorData << "Values = " << value[0] << " " << value[1]  \
                   << " " << value[2] << std::endl;
        
        getline(inputData,line);getline(inputData,line);getline(inputData,line);
        getline(inputData,line);getline(inputData,line);

        mirrorData << "Connectivity = " << std::endl;
        
        for (int i=0; i<numGroup; i++){
            inputData >> connectB(0) >> connectB(1) >> connectB(2) \
                      >> connectB(3) >> connectB(4) >> connectB(5);

            connectB(0) -= 1;
            connectB(1) -= 1;
            connectB(2) -= 1;
            connectB(3) -= 1;
            connectB(4) -= 1;
            connectB(5) -= 1;

            Boundaries *bound = new Boundaries(connectB, index++,  \
                                               constrain, value);
            boundary_.push_back(bound);

            mirrorData << connectB(0) << " " << connectB(1) << " " \
                       << connectB(2) << " " << connectB(3) << " " \
                       << connectB(4) << " " << connectB(5) << std::endl;
        };
    };

    for (int ibound = 0; ibound < numBoundElems; ibound++){
        
        Boundaries::BoundConnect connectB;
        connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB(0);
        int no2 = connectB(1);
        int no3 = connectB(2);
        int no4 = connectB(3);
        int no5 = connectB(4);
        int no6 = connectB(5);
        
        if (boundary_[ibound] -> getConstrain(0)){
            nodes_[no1] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no2] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no3] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no4] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no5] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
            nodes_[no6] -> setConstrains(0,boundary_[ibound]        \
                                         -> getConstrainValue(0));
        };
        if (boundary_[ibound] -> getConstrain(1)){
            nodes_[no1] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no2] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no3] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no4] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no5] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
            nodes_[no6] -> setConstrains(1,boundary_[ibound]        \
                                         -> getConstrainValue(1));
        };       
        if (boundary_[ibound] -> getConstrain(2)){
            nodes_[no1] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
            nodes_[no2] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
            nodes_[no3] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
            nodes_[no4] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
            nodes_[no5] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
            nodes_[no6] -> setConstrains(2,boundary_[ibound]        \
                                         -> getConstrainValue(2));
        };               
    };
    
    //Sets nodal constrains
    for (int i=0; i<numNodes; i++){
        mirrorData<< "Constrains " << i << " " << nodes_[i]->getConstrains(0) \
                  << " " << nodes_[i]->getConstrainValue(0)             \
                  << " " << nodes_[i]->getConstrains(1)                 \
                  << " " << nodes_[i]->getConstrainValue(1)             \
                  << " " << nodes_[i]->getConstrains(2)                 \
                  << " " << nodes_[i]->getConstrainValue(2) << std::endl;
    };

};

//------------------------------------------------------------------------------
//--------------------------SOLVE STEADY FLUID PROBLEM--------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveSteadyProblem(int iterNumber, double tolerance,\
                                 int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione, iterations;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    idx_t* part_elem = domainDecompositionMETIS();

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, \
                            2*numNodes+numNodesP, 2*numNodes+numNodesP, \
                            numNodesP,NULL,numNodesP*3,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodesP);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                if (problem_type == 1) elements_[jel] -> getSteadyStokes();
                
                if (problem_type == 2) {
                    if (inewton == 0) {
                        elements_[jel] -> getSteadyNavierStokes();
                    } else {
                        elements_[jel] -> getSteadyNavierStokes();
                    };
                };
                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                typename Elements::ConnectivityP connecP;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                connecP = elements_[jel] -> getConnectivityP();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j+1),ADD_VALUES);
                        };
                    };
                    
                    //Matrix Q and Qt
                    for (int j=0; j<3; j++){
                        if (fabs(Ajac(2*i  ,12+j)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,12+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(12+j,2*i  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,12+j)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,12+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(12+j,2*i+1),ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                    };
                    
                    if (fabs(Rhs(2*i+1)) >= 1.e-8){
                        int dof_i = 2*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);
                    };
                };
                for (int i=0; i<3; i++){
                    if (fabs(Rhs(12+i)) >= 1.e-8){
                        int dof_i = 2*numNodes + connecP(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-7,1.e-10,PETSC_DEFAULT,
                                10000);CHKERRQ(ierr);
        
        //        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
      
        ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        
        ierr = PCSetType(pc, PCJACOBI);CHKERRQ(ierr);
        
        //ierr = KSPSetPCSide(ksp, PC_RIGHT);
        //ierr = KSPSetType(ksp,KSPTFQMR); CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        ierr = KSPGetTotalIterations(ksp, &iterations);

        std::cout << "GMRES Iterations = " << iterations << std::endl;
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double u_ [2];
        double p_;
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            Ii = 2*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 2*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            nodes_[i] -> incrementVelocity(0,u_[0]);
            nodes_[i] -> incrementVelocity(1,u_[1]);
        };
        for (int i = 0; i<numNodesP; i++){
            Ii = 2*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodesp_[i] -> incrementPressure(p_);
        };
        
        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "Du Norm = " << val << std::endl;
            std::cout<<"Iteration = " << inewton << std::endl;
        };

        if(val <= tolerance){
            break;            
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
             
    };

    if (rank == 0) {
        printVelocity(1);
        printPressure(1);
    };

    return 0;
};

template<>
int Fluid<3>::solveSteadyProblem(int iterNumber, \
                                 double tolerance, int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    int               rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    idx_t* part_elem = domainDecompositionMETIS();

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    for (int inewton = 0; inewton < iterNumber; inewton++){
        //Create PETSc matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, \
                            3*numNodes+numNodesP, 3*numNodes+numNodesP, \
                            numNodesP*8,NULL,numNodesP*8,NULL,&A);CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,3*numNodes+numNodesP);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //    std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;      
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                if (problem_type == 1) elements_[jel] -> getSteadyStokes();
                
                if (problem_type == 2) {
                    if (inewton == 0) {
                        elements_[jel] -> getSteadyNavierStokes();
                    } else {
                        elements_[jel] -> getSteadyNavierStokes();
                    };
                };


                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                typename Elements::ConnectivityP connecP;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                connecP = elements_[jel] -> getConnectivityP();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<10; i++){
                    for (int j=0; j<10; j++){
                        if (fabs(Ajac(3*i  ,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i  ,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i  ,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j+2),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j+2),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j+2),ADD_VALUES);
                        };
                    };
                    
                    //Matrix Q and Qt
                    for (int j=0; j<4; j++){
                        if (fabs(Ajac(3*i  ,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i+2),ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector
                    if (fabs(Rhs(3*i  )) >= 1.e-8){
                        int dof_i = 3*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i  ),ADD_VALUES);
                    };
                    if (fabs(Rhs(3*i+1)) >= 1.e-8){
                        int dof_i = 3*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i+1),ADD_VALUES);
                    };
                    if (fabs(Rhs(3*i+2)) >= 1.e-8){
                        int dof_i = 3*connec(i)+2;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i+2),ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
                                30000);CHKERRQ(ierr);
        
        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        
        ierr = KSPGetPC(ksp,&pc);
        
        ierr = PCSetType(pc, PCNONE);
        
        //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
               
        //Updates nodal values
        double u_ [3];
        double p_;
        Ione = 1;
        for (int i = 0; i<numNodes; i++){
            Ii = 3*i;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[0] = val;
            Ii = 3*i+1;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[1] = val;
            Ii = 3*i+2;
            ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
            u_[2] = val;
            nodes_[i] -> incrementVelocity(0,u_[0]);
            nodes_[i] -> incrementVelocity(1,u_[1]);
            nodes_[i] -> incrementVelocity(2,u_[2]);
        };
        for (int i = 0; i<numNodesP; i++){
            Ii = 3*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodesp_[i] -> incrementPressure(p_);
        };

        //Computes the solution vector norm
        ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout << "Du Norm = " << val << std::endl;
            std::cout<<"Iteration = " << inewton << std::endl;
        };

        if(val <= tolerance){
            break;            
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
             
    };

    if (rank == 0) {
        printVelocity(1);
        printPressure(1);
    };

    return 0;
};

//------------------------------------------------------------------------------
//-------------------------SOLVE TRANSIENT FLUID PROBLEM------------------------
//------------------------------------------------------------------------------
template<>
int Fluid<2>::solveTransientProblem(int iterNumber, double tolerance,\
                                 int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
   
    int rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    idx_t* part_elem = domainDecompositionMETIS();

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };
        
    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) std::cout << "----------- TIME STEP = "          \
                                 << iTimeStep << " ------------" << std::endl;

        for (int i = 0; i < numNodes; i++){
            double accel[2], u[2], uprev[2];
            
            //Compute acceleration
            u[0] = nodes_[i] -> getVelocity(0);
            u[1] = nodes_[i] -> getVelocity(1);

            uprev[0] = nodes_[i] -> getPreviousVelocity(0);
            uprev[1] = nodes_[i] -> getPreviousVelocity(1);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            
            nodes_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodes_[i] -> setPreviousVelocity(u);

            //Computes predictor
            u[0] += dTime * accel[0];
            u[1] += dTime * accel[1];

            nodes_[i] -> setVelocity(u);
            
            //Updates acceleration
            nodes_[i] -> setPreviousAcceleration(accel);

        };
        
        
        
        

    for (int inewton = 0; inewton < iterNumber; inewton++){

        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, \
                            2*numNodes+numNodesP, 2*numNodes+numNodesP, \
                            numNodesP,NULL,numNodesP*3,NULL,&A); CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,2*numNodes+numNodesP);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                if (problem_type == 1) elements_[jel] -> getTransientStokes();
                
                if (problem_type == 2) elements_[jel] -> \
                                           getTransientNavierStokes();
                
                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                typename Elements::ConnectivityP connecP;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                connecP = elements_[jel] -> getConnectivityP();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<6; i++){
                    for (int j=0; j<6; j++){
                        if (fabs(Ajac(2*i  ,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j  )) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i  ,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,2*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,2*j+1)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,2*j+1),ADD_VALUES);
                        };
                    };
                    
                    //Matrix Q and Qt
                    for (int j=0; j<3; j++){
                        if (fabs(Ajac(2*i  ,12+j)) >= 1.e-8){
                            int dof_i = 2*connec(i);
                            int dof_j = 2*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i  ,12+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(12+j,2*i  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(2*i+1,12+j)) >= 1.e-8){
                            int dof_i = 2*connec(i)+1;
                            int dof_j = 2*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(2*i+1,12+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(12+j,2*i+1),ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector
                    if (fabs(Rhs(2*i  )) >= 1.e-8){
                        int dof_i = 2*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i  ),ADD_VALUES);
                    };
                    
                    if (fabs(Rhs(2*i+1)) >= 1.e-8){
                        int dof_i = 2*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(2*i+1),ADD_VALUES);
                    };
                };
                for (int i=0; i<3; i++){
                    if (fabs(Rhs(12+i)) >= 1.e-8){
                        int dof_i = 2*numNodes + connecP(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(12+i),ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
                                PETSC_DEFAULT);CHKERRQ(ierr);
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        
        ierr = KSPGetPC(ksp,&pc);
        
        ierr = PCSetType(pc, PCJACOBI);

        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);

        //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
                
        //Updates nodal values
        double p_;
        double duNorm = 0.;
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            if (nodes_[i] -> getConstrains(0) == false){
                Ii = 2*i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                duNorm += val*val;
                nodes_[i] -> incrementVelocity(0,val);
            }; 

            if (nodes_[i] -> getConstrains(1) == false){
                Ii = 2*i+1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                duNorm += val*val;
                nodes_[i] -> incrementVelocity(1,val);
            };
        };


        

        for (int i = 0; i<numNodesP; i++){
            Ii = 2*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodesp_[i] -> incrementPressure(p_);
        };
        
        //Computes the solution vector norm
        //ierr = VecNorm(u,NORM_2,&val);CHKERRQ(ierr);
        
        if(rank == 0){
            std::cout<<"Iteration = " << inewton << \
                "      Du Norm = " << sqrt(duNorm) << " " << tolerance << std::endl;
        };

        if(sqrt(duNorm) <= tolerance){
            break;            
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
             
    };

    if (rank == 0) {
        printVelocity(iTimeStep);
        printPressure(iTimeStep);
    };

    };

    return 0;
};

template<>
int Fluid<3>::solveTransientProblem(int iterNumber, \
                                 double tolerance, int problem_type) {

    Mat               A;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ione;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    int               rank;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    idx_t* part_elem = domainDecompositionMETIS();

    //Check if the problem type can be computed
    if ((problem_type > 2) || (problem_type < 1)){
        std::cout << "WRONG PROBLEM TYPE." << std::endl;
        return 0;
    };

    for (int iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
        
        if (rank == 0) std::cout << "----------- TIME STEP = "          \
                                 << iTimeStep << " ------------" << std::endl;

        for (int i = 0; i < numNodes; i++){
            double accel[3], u[3], uprev[3];
            
            //Compute acceleration
            u[0] = nodes_[i] -> getVelocity(0);
            u[1] = nodes_[i] -> getVelocity(1);
            u[2] = nodes_[i] -> getVelocity(2);

            uprev[0] = nodes_[i] -> getPreviousVelocity(0);
            uprev[1] = nodes_[i] -> getPreviousVelocity(1);
            uprev[2] = nodes_[i] -> getPreviousVelocity(2);
            
            accel[0] = (u[0] - uprev[0]) / dTime;
            accel[1] = (u[1] - uprev[1]) / dTime;
            accel[2] = (u[2] - uprev[2]) / dTime;
            
            nodes_[i] -> setAcceleration(accel);
            
            //Updates velocity
            nodes_[i] -> setPreviousVelocity(u);

            //Computes predictor
            u[0] += dTime * accel[0];
            u[1] += dTime * accel[1];
            u[2] += dTime * accel[2];

            nodes_[i] -> setVelocity(u);
            
            //Updates acceleration
            nodes_[i] -> setPreviousAcceleration(accel);

        };


    for (int inewton = 0; inewton < iterNumber; inewton++){
        //Create PETSc matrix
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, \
                            3*numNodes+numNodesP, 3*numNodes+numNodesP, \
                            numNodesP*8,NULL,numNodesP*8,NULL,&A);CHKERRQ(ierr);
        
        ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);
        
        //Create PETSc vectors
        ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
        ierr = VecSetSizes(b,PETSC_DECIDE,3*numNodes+numNodesP);CHKERRQ(ierr);
        ierr = VecSetFromOptions(b);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
        ierr = VecDuplicate(b,&All);CHKERRQ(ierr);
        
        //    std::cout << "Istart = " << Istart << " Iend = " << Iend << std::endl;      
        
        for (int jel = 0; jel < numElem; jel++){   
            
            if (part_elem[jel] == rank) {
                
                //Compute Element matrix
                if (problem_type == 1) elements_[jel] -> getSteadyStokes();
                
                if (problem_type == 2) {
                    if (inewton == 0) {
                        elements_[jel] -> getSteadyStokes();
                    } else {
                        elements_[jel] -> getSteadyNavierStokes();
                    };
                };


                typename Elements::LocalMatrix Ajac;
                typename Elements::LocalVector Rhs;
                typename Elements::Connectivity connec;
                typename Elements::ConnectivityP connecP;
                
                //Gets element connectivity, jacobian and rhs 
                connec = elements_[jel] -> getConnectivity();
                connecP = elements_[jel] -> getConnectivityP();
                Ajac = elements_[jel] -> getJacNRMatrix();
                Rhs = elements_[jel] -> getRhsVector();
                
                //Disperse local contributions into the global matrix
                //Matrix K and C
                for (int i=0; i<10; i++){
                    for (int j=0; j<10; j++){
                        if (fabs(Ajac(3*i  ,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j  )) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i  ,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j+1)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j)+1;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i  ,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,3*j+2),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,3*j+2),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,3*j+2)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*connec(j)+2;
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,3*j+2),ADD_VALUES);
                        };
                    };
                    
                    //Matrix Q and Qt
                    for (int j=0; j<4; j++){
                        if (fabs(Ajac(3*i  ,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i);
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i  ,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i  ),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+1,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i)+1;
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+1,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i+1),ADD_VALUES);
                        };
                        if (fabs(Ajac(3*i+2,30+j)) >= 1.e-8){
                            int dof_i = 3*connec(i)+2;
                            int dof_j = 3*numNodes + connecP(j);
                            ierr = MatSetValues(A,1,&dof_i,1,&dof_j,    \
                                                &Ajac(3*i+2,30+j),ADD_VALUES);
                            ierr = MatSetValues(A,1,&dof_j,1,&dof_i,    \
                                                &Ajac(30+j,3*i+2),ADD_VALUES);
                        };
                    };
                    
                    //Rhs vector
                    if (fabs(Rhs(3*i  )) >= 1.e-8){
                        int dof_i = 3*connec(i);
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i  ),ADD_VALUES);
                    };
                    if (fabs(Rhs(3*i+1)) >= 1.e-8){
                        int dof_i = 3*connec(i)+1;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i+1),ADD_VALUES);
                    };
                    if (fabs(Rhs(3*i+2)) >= 1.e-8){
                        int dof_i = 3*connec(i)+2;
                        ierr = VecSetValues(b,1,&dof_i,&Rhs(3*i+2),ADD_VALUES);
                    };
                };
            };
        };
        
        //Assemble matrices and vectors
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
        
        //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        
        //Create KSP context to solve the linear system
        ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        
        ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
        
        ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,
                                30000);CHKERRQ(ierr);
        
        ierr = KSPGMRESSetRestart(ksp, 10); CHKERRQ(ierr);
        
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        
        ierr = KSPGetPC(ksp,&pc);
        
        ierr = PCSetType(pc, PCNONE);
        
        //    ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
        
        //ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKERRQ(ierr);
        
        //Gathers the solution vector to the master process
        ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
        
        ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
        
        ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
               
        //Updates nodal values
        double p_;
        double duNorm = 0.;
        Ione = 1;

        for (int i = 0; i < numNodes; ++i){
            if (nodes_[i] -> getConstrains(0) == false){
                Ii = 3*i;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                duNorm += val*val;
                nodes_[i] -> incrementVelocity(0,val);
            }; 

            if (nodes_[i] -> getConstrains(1) == false){
                Ii = 3*i+1;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                duNorm += val*val;
                nodes_[i] -> incrementVelocity(1,val);
            };

            if (nodes_[i] -> getConstrains(2) == false){
                Ii = 3*i+2;
                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
                duNorm += val*val;
                nodes_[i] -> incrementVelocity(2,val);
            };
        };

        for (int i = 0; i<numNodesP; i++){
            Ii = 3*numNodes+i;
            ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
            p_ = val;
            nodesp_[i] -> incrementPressure(p_);
        };

        if(rank == 0){
            std::cout<<"Iteration = " << inewton << \
                "      Du Norm = " << sqrt(duNorm) <<std::endl;
        };

        if(sqrt(duNorm) <= tolerance){
            break;            
        };

        ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
        ierr = VecDestroy(&b); CHKERRQ(ierr);
        ierr = VecDestroy(&u); CHKERRQ(ierr);
        ierr = VecDestroy(&All); CHKERRQ(ierr);
        ierr = MatDestroy(&A); CHKERRQ(ierr);
             
    };

    if (rank == 0) {
        printVelocity(iTimeStep);
        printPressure(iTimeStep);
    };

    };
    return 0;
};



#endif
