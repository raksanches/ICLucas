#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>

class connection

{
    public:
    double a, b;

};

class node

{
    public:
    double x, y, z, rot;

};


int main(){

    double pi = 4.*atan(1);
    int numElem, numNodes, numTimesteps;
    double x[3],v[3],k[3],c[3],cForce[3],vForce[3],w[3],beta,gama,timestep,m;
    std::string input_file, output_file,line, output_file_2;
    std::vector<double> cx;
    std::vector<double> cy;
    std::vector<double> cz;
    std::vector<double> conec1;
    std::vector<double> conec2; 

    output_file_2 = "displacement.dat";
    output_file = "mirror.txt";
    input_file = "input_test.txt";

    std::ifstream inputData(input_file.c_str());
    std::ofstream outputData(output_file.c_str());
    std::ofstream output(output_file_2.c_str());

    getline(inputData,line),getline(inputData,line),getline(inputData,line);

    //aspectos gerais
    inputData >> numElem >> numNodes >> numTimesteps;

    outputData << "Número de elementos: " << numElem << std::endl;
    outputData << "Número de nós: " << numNodes << std::endl;
    outputData << "Número de timesteps: " << numTimesteps << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);;

    //condicoes iniciais    
    inputData >> v[0] >> v[1] >> v[2];

    outputData << "Velocidade em X: " << v[0] << std::endl; 
    outputData << "Velocidade em Y: " << v[1] << std::endl; 
    outputData << "Velocidade rotacional: " << v[2] << std::endl; 

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);

    //inputData >> a[0] >> a[1] >> a[2];

    //outputData << "Aceleracao inicial em X: " << a[0] << std::endl; 
    //outputData << "Aceleracao inicial em Y: " << a[1] << std::endl; 
    //outputData << "Aceleracao radial inicial: " << a[2] << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);

    inputData >> x[0] >> x[1] >> x[2];
    
    outputData << "Deslocamento inicial X: " << x[0] << std::endl; 
    outputData << "Descolamento inicial Y: " << x[1] << std::endl; 
    outputData << "Rotacao inicial: " << x[2] << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //caracteristicas da secao
    inputData >> m >> k[0] >> k[1] >> k[2] >> c[0] >> c[1] >> c[2];

    outputData << "Massa: " << m << std::endl;
    outputData << "K X: " << k[0] << std::endl; 
    outputData << "K Y: " << k[1] << std::endl; 
    outputData << "Kteta: " << k[2] << std::endl;
    outputData << "Cx " << c[0] << std::endl;
    outputData << "Cy " << c[1] << std::endl;
    outputData << "Cteta " << c[2] << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //metodo numerico
    inputData >> timestep >> gama >> beta;

    outputData << "Timestep: " << timestep << std::endl; 
    outputData << "Gama: " << gama << std::endl; 
    outputData << "Beta: " << beta << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);   
    getline(inputData,line); 

    //forcas aplicadas
    inputData >> cForce[0] >> vForce[0] >> w[0] >> cForce[1] >> vForce[1] >> w[1] >> cForce[2] >> vForce[2] >> w[2];

    outputData << "Fx(A+Bsinwt) A:" << cForce[0]<< " B: " << vForce[0] << " w: " << w[0] << std::endl; 
    outputData << "Fy(A+Bsinw) A:" << cForce[1]<< " B: " << vForce[1] << " w: " << w[1] << std::endl;
    outputData << "My(A+Bsinwt) A:" << cForce[2]<< " B:" << vForce[2] << " w:" << w[2] << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //coordenadas
    outputData << "Coordenadas \n";
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

    //conectividade
    outputData << "Conectividade \n";
    for (int i=0; i<4; i++){

        connection connection;
        inputData >> connection.a >> connection.b;
        //outputData << conection.a << " " << conection.b;
        conec1.push_back(connection.a);
        conec2.push_back(connection.b);

    }

    for (int i=0; i<4; i++){

        outputData << conec1[i] << " " << conec2[i] << "\n\n";

    }
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //centro da secao 
    node cscenter_initial;
    node cscenter;

    inputData >> cscenter_initial.x >> cscenter_initial.y >> cscenter_initial.z >> cscenter_initial.rot;
    outputData << "Centro \n";
    outputData << cscenter_initial.x << " " << cscenter_initial.y << " " << cscenter_initial.z << " \n\n";

    cscenter.x = cscenter_initial.x;
    cscenter.y = cscenter_initial.y;
    cscenter.z = cscenter_initial.z;

    //calculo do alfa
    std::vector<double> alfa;
    std::vector<double> alfaInit;
    std::vector<double> sectionNodesRadius;

     for (int i=0; i<numNodes; i++){

        double alfa_node;
        int quad;
        double delta_cx = cx[i] - cscenter_initial.x;
        double delta_cy = cy[i] - cscenter_initial.y;
        double l_node = sqrt(delta_cx*delta_cx+delta_cy*delta_cy);
        sectionNodesRadius.push_back(l_node);

        if (delta_cx > 0.){
            if (delta_cy> 0.){
                alfa_node = asin(delta_cy/l_node);
                quad = 1;
            }else{
                alfa_node =  2. * pi + asin(delta_cy/l_node);
                quad = 4;
            }
        }else{
            if (delta_cy>0){
                alfa_node = pi - asin(delta_cy/l_node);
                quad = 2;
            }else{
                alfa_node = pi - asin(delta_cy/l_node);
                quad = 3;
            }
        }
        
        alfa.push_back(alfa_node);
        alfaInit.push_back(alfa_node);
    }
    //ACELERACAO INICIAL
    double t = 0.;

    std::vector<double> a;

    for (int i=0; i<3; i++){

        double aceleration = (1/m)*(-k[i]*x[i]-c[i]*v[i]+(cForce[i]+vForce[i]*sin(0)));
        a.push_back(aceleration);
        outputData << "a(" << i << ") = " << aceleration <<  " \n";
    }


    //MOVIMENTACAO DO CENTRO DE MASSA:
    double xPrev[3], vPrev[3], aPrev[3];
    for (int i=0; i<3; i++){
	
			xPrev[i] = x[i];
			vPrev[i] = v[i];
			aPrev[i] = a[i];	
		}

    outputData << "\nTime: " << t << " delta_x: " << x[0] << " delta_y: " << x[1] << " teta: " << x[2] << std::endl;

    for (int i=0; i<numNodes; i++){

            outputData << "Node: " << i << " x: " << cx[i] << " y: " << cy[i] << " \n";

    }

    for (int iTime=1; iTime<=numTimesteps; iTime++){

		t=t+timestep;

        double force[3];
        for (int i=0;i<3;i++){

                force[i]=cForce[i]+vForce[i]*sin(w[i]*t);
        }

        for (int i=0; i<3;i++){


                x[i] = (force[i] + ((1./beta)*(m/(timestep*timestep))+(gama/beta)*(c[i]/timestep))*xPrev[i]
                +((1./beta)*(m/timestep)+((gama/beta)-1.)*c[i])*vPrev[i]+((1./(2.*beta)-1.)*m+(gama/(2.*beta)-1.)*c[i]*timestep)*aPrev[i])
                /((1./beta)*(m/(timestep*timestep))+k[i]+(gama/beta)*(c[i]/timestep));

                v[i]=((gama/beta)/timestep)*x[i]-((gama/beta)/timestep)*xPrev[i]+(1-(gama/beta))*vPrev[i]
                +(1-(gama/(2.*beta)))*timestep*aPrev[i];

                a[i] = (force[i] - k[i]*x[i]-c[i]*v[i])/m;
 
                xPrev[i] = x[i];
                vPrev[i] = v[i];
                aPrev[i] = a[i];	
                    
        }
    
        outputData << "\n Time: " << t << " delta_x: " << x[0] << " delta_y: " << x[1] << " delta_teta: " << x[2] << " \n\n";
        outputData << " velocidade:" << v[0] << "aceleracao: " << a[0] << " \n\n";

        //NOVAS COORDENDAS!!

        //centro.rot = centro.rot + x[2];
        //outputData << "Desloc x: " << x[0] << " y: " << x[1] << " teta: " << x[2] << "\n";

        //nodes
        for (int i=0; i<numNodes; i++){
            
            double l_node = sectionNodesRadius[i]; 
            //novo alfa
            alfa[i]=alfaInit[i]+x[2];
            //efeito da rotacao
            cx[i]= cscenter.x+l_node*cos(alfa[i]);
            cy[i]= cscenter.y+l_node*sin(alfa[i]);
            //translacao nodes
            cx[i] += x[0];
            cy[i] += x[1];           
            //centro

            outputData << "Node: " << i << " x: " << cx[i] << " y: " << cy[i] << " \n";

        }
            cscenter.x = cscenter_initial.x + x[0];
            cscenter.y = cscenter_initial.y + x[1];

        std::string result;
        std::ostringstream convert;

        convert << iTime + 10000000;

        result = convert.str();
        std::string s = "crossSection"+result+".vtu";
    
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
            output_v << cx[i] << " " << cy[i] << " " << 0.0 << std::endl;        
        };
        output_v << "      </DataArray>" << std :: endl
               << "    </Points>" << std :: endl;

        // write element connectivity
        output_v << "    <Cells>" << std :: endl
           << "      <DataArray type=\"Int32\" "
           << "Name=\"connectivity\" format=\"ascii\">" << std :: endl;
    
        for (int i=0; i<numElem; i++){
        
        output_v << conec1[i] << " " << conec2[i] << std::endl;
        };
        output_v << "      </DataArray>" << std :: endl;
  
        // write offsets in the data array
        output_v << "      <DataArray type=\"Int32\""
           << " Name=\"offsets\" format=\"ascii\">" << std :: endl;

        int aux = 0;
        for (int i=0; i<numElem; i++){
            output_v << aux + 2 << std::endl;
            aux += 2;
        };
        output_v << "      </DataArray>" << std :: endl;
  
        // write element types
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
               << "format=\"ascii\">" << std :: endl;
    
        for (int  i=0; i<numElem; i++){
            output_v << 4 << std::endl;
        };

        output_v << "      </DataArray>" << std :: endl
               << "    </Cells>" << std :: endl;


//     //Write point results
//     output_v << "    <PointData>" << std :: endl;
//     output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
// 	 << "Name=\"Velocity\" format=\"ascii\">" << std :: endl;
//     for (int i=0; i<numNodes; i++){
//         output_v << nodes_[i] -> getVelocity(0) << " "          \
//                  << nodes_[i] -> getVelocity(1) << " " << 0.0 << std::endl;
//     };
//     output_v << "      </DataArray> " << std::endl;

//    //  output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
// 	 // << "Name=\"Stress\" format=\"ascii\">" << std :: endl;
//    //  for (int i=0; i<Nnos; i++){
//    //      output_v << 0.0 << " " << 0.0 << " " << 0.0 << endl;
//    //  };
    
//    //  output_v << "      </DataArray> " << std::endl;

//     output_v << "    </PointData>" << std :: endl; 

//     //Write cell results
//     output_v << "    <CellData>" << std :: endl;

//     output_v << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
// 	 << "Name=\"Process\" format=\"ascii\">" << std :: endl;
//     for (int i=0; i<numElem; i++){
//         output_v << part_elem[i] << std::endl;
//     };
    
//     output_v << "      </DataArray> " << std::endl;
//     output_v << "    </CellData>" << std :: endl; 



    // finalise
    output_v << "  </Piece>" << std :: endl
           << "  </UnstructuredGrid>" << std :: endl
           << "</VTKFile>" << std :: endl;



        output << t << " "<< x[0] << std::endl; 
    } //t (Time step)

    return 0;

}
