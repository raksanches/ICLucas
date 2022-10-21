#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>

class node

{
    public:
    double x, y, z;

};

int main(){

    int numElem, numNodes, numTimesteps,dimension;
    std::vector<double> cx;
    std::vector<double> cy;
    std::vector<double> cz;

    std::vector<int> conec; 
    std::string input_file, output_file,line;

    output_file = "mirror_coords.txt";
    input_file = "coords.txt";

    std::ifstream inputData(input_file.c_str());
    std::ofstream outputData(output_file.c_str());

    outputData << "Coordenadas  \n\n"; 

    //coordenadas

    for (int i=0; i<4; i++){

        node node;
        inputData >> node.x >> node.y >> node.z;
        //outputData << node.x << " " << node.y << " "<< node.z << " \n\n";
        cx.push_back(node.x);
        cy.push_back(node.y);
        cz.push_back(node.z);
    }

    for (int i=0; i<4; i++){

        outputData << cx[i] << " " << cy[i] << " " << cz[i] << "\n\n";

    }

    getline(inputData,line);getline(inputData,line);getline(inputData,line);getline(inputData,line);

    node centro;

    inputData >> centro.x >> centro.y >> centro.z;
    outputData << centro.x << " " << centro.y << " " << centro.z;

    return 0;
}