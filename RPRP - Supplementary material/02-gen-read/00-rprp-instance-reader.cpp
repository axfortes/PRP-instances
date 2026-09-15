#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <string>
#include <string.h>
#include <deque>
#include <random>
#include "boost/dynamic_bitset.hpp"
#include <chrono>
#include <ctime>

//bibliotecas cplex
#include <ilcplex/ilocplex.h>
#include <ilcplex/cplexx.h>
#include <ilcplex/cplex.h>

#define myrand ((float)(random())/(float)(RAND_MAX) )

using namespace std;
using namespace chrono;

ILOSTLBEGIN
//coloring the text
char DGRAY[]={0x1b,'[','1',';','3','0','m',0};

char red[]={0x1b,'[','0',';','3','1','m',0};
char RED[]={0x1b,'[','1',';','3','1','m',0};
//char uRED[]={0x1b,'[','4',';','3','1','m',0};

char green[]={0x1b,'[','0',';','3', '2','m',0};
char GREEN[]={0x1b,'[','1',';','3', '2','m',0};

char yellow[]={0x1b,'[','0',';','3', '3', 'm',0};
char YELLOW[]={0x1b,'[','1',';','3', '3', 'm',0};

char blue[]={0x1b,'[','0',';','3','4','m',0};
char BLUE[]={0x1b,'[','1',';','3','4','m',0};

//char Upurple[]={0x1b,'[','4',';','3','5','m',0};
char purple[]={0x1b,'[','0',';','3','5','m',0};
char PURPLE[]={0x1b,'[','1',';','3','5','m',0};

char cyan[]={0x1b,'[','0',';','3','6','m',0};
char CYAN[]={0x1b,'[','1',';','3','6','m',0};

char lgray[]={0x1b,'[','0',';','3','7','m',0};
char LGRAY[]={0x1b,'[','1',';','3','7','m',0};

char normal[]={0x1b,'[','0',';','3','9','m',0};
char NORMAL[]={0x1b,'[','1',';','3','9','m',0};


//instance data
typedef struct{
    int n;      //numero de clientes
    int P;      //numeto de produtos
    int T;      //numero de periodos;

    //informacoes sobre a planta produtiva
    std::vector<int> C;     //capacidade de produção de cada fabrica de cada item
    std::vector<double> l;  //custo de setup/ativacao de cada item em cada fabrica
    std::vector<double> u;  //custo unitário de produção de cada item em cada fabrica;

    //informacoes dos clientes
    std::vector<std::vector<std::vector<int> > >d;   //demanda periodica de cada item por cada cliente
    std::vector<std::vector<double> > a;   //tempo de viagem de i para j
    std::vector<int> s;                //tempo de servico de i
    std::vector<std::vector<double> > B;//custo de atrasar o produto p para o cliente i
    std::vector<std::vector<int> > D; //demanda total

    //informacoes da frota
    int V;                              //numero de veiculos
    std::vector<int> Q;                 //capacidade maxima de transporte de determinado veiculo
    std::vector<int> e;                 //custo de setup/ativação de cada veiculo
    int H;                              //horizonte de tempo para entrega, de 12 horas (6h - 18h)
    int Qmax;

    //informacoes gerais
    std::vector<std::vector<double> > c; //custo de transporte de i para j
    std::vector<double> x;               //coordenada x do no i
    std::vector<double> y;               //coordenada y do no i
    std::vector<std::vector<double> > h; //custo de estocagem do produto p no no i
    std::vector<std::vector<int> > U;    //limite superior de estocagem no item de cada planta
    std::vector<std::vector<int> > I0;   //inventario inicial em cada no, deve satisfazer a demanda do periodo

    string name;
    string pathIn;
    string pathOut;


}data;

void help();
void read(data &d, string sPath, string sName);

void help(){
  std::cout << std::endl << std::endl << "exec [data file] \n " << std::endl;
  exit(1);
}

void read(data &d){
    string sPath = d.pathIn + d.name;
    //std::cout << "\n\tPath = " << sPath;
    char *cPath = new char[sPath.length()+1];
    memcpy(cPath, sPath.c_str(), sPath.length() + 1);
    //std::cout << "\n\tPath = " << cPath;
    //fim conversao string para char

    //abertura do arquivo de dados
    ifstream arq(cPath);
    if (!arq.is_open()) help();

    //getchar();

    //iniciando leitura dos parametros
    arq >> d.n;      //number of clients
    arq >> d.T;      //number of periods
    arq >> d.P;      //number of commodities
    arq >> d.V;      //number of vehicles
    arq >> d.H;      //maximum time route size

    //reservando posicoes 0 e n+1 para nos artificiais
    d.x = std::vector<double> (d.n + 2,0.0);  //coordenada x do no i
    d.y = std::vector<double> (d.n + 2,0.0);  //coordenada x do no i
    d.c = std::vector<std::vector<double> > (d.n + 2, std::vector<double> (d.n + 2, 0.0));  //custo (i,j)
    d.a = std::vector<std::vector<double> > (d.n + 2, std::vector<double> (d.n + 2, 0.0));  //tempo (i,j)
    d.s = std::vector<int> (d.n + 1, 10); //tempo de atendimento

    //estoque inicial
    d.I0 = std::vector<std::vector<int> > (d.P + 1, std::vector<int> (d.n + 1,0));

    //demanda periodica e acumulada
    d.d = std::vector<std::vector<std::vector<int> > > (d.T + 1, std::vector<std::vector <int> > (d.P + 1, std::vector<int>(d.n + 1, 0)));
    //d.D = std::vector<std::vector<int> > (d.P + 1, std::vector<int>(d.n + 1, 0)); //original
    d.D = std::vector<std::vector<int> > (d.T + 1, std::vector<int>(d.P + 1, 0));

    //informacoes sobre producao e fabrica
    d.u = std::vector<double> (d.P + 1,0); //custo de producao unitario
    d.C = std::vector<int> (d.P + 1,0);       //capacidade de cada linha de produção
    d.l = std::vector<double> (d.P + 1,0); //custo de preparação de produção (setup)
    d.medCostP = 0.0;

    //informacoes sobre os veiculos
    d.Q = std::vector<int> (d.V + 1, 0);        //carga maxima carregada
    d.e = std::vector<int> (d.V + 1, 0);        //custo de ativação

    //custos e limites de estocagem e custo de atraso
    d.h = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n + 1,0.0));
    d.U = std::vector<std::vector<int> > (d.P + 1, std::vector<int> (d.n + 1, 0));
    d.B = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n + 1,0.0));


    //informacoes sobre os veiculos
    d.Qmax = 0;        
    for (int v = 1; v <= d.V; v++){
        arq >> d.Q[v];
        arq >> d.e[v];
        if (d.Q[v] > d.Qmax) d.Qmax = d.Q[v];        
        d.maxLoad += d.Q[v];
    }    

    //lendo informacoes dos nos clientes e plantas
    for (int i = 0; i <= d.n; i++){
        arq >> d.x[i];
        arq >> d.y[i];
        //std::cout << "\ni: " <<  i << " - x: " << d.x[i] << " - y: " << d.y[i] << " ";
        if (i == 0){
            for (int p = 1; p <= d.P; p++){
                arq >> d.I0[p][i];
                arq >> d.h[p][i];
                arq >> d.U[p][i];
                arq >> d.C[p];
                arq >> d.l[p];
                arq >> d.u[p];
                d.medCostP += d.l[p];
                //std::cout << "\n>>p: " << p << " - d.I0: " << d.I0[p][i] << " - h: " << d.h[p][i] << " - U: " << d.U[p][i] << " - u: " << d.u[p] << " - l: " << d.l[p] << " - C: " << d.C[p];
            }
            d.x[d.n+1] = d.x[i];
            d.y[d.n+1] = d.y[i];
        }else if (i > 0){
            arq >> d.s[i];
            for (int p = 1; p <= d.P; p++){
                arq >> d.I0[p][i];
                arq >> d.h[p][i];
                arq >> d.U[p][i];
                arq >> d.B[p][i];
                //std::cout << "\n>>p: " << p << " - d.I0: " << d.I0[p][i] << " - h: " << d.h[p][i] << " - U: " << d.U[p][i] << " - B: " << d.B[p][i];
            }
        }
    }

    //demandas de cada cliente por cada produto em cada periodo
    for (int p = 1; p <= d.P; p++)
        for (int i = 1; i <= d.n; i++)
            for (int t = 1; t <= d.T; t++) arq >> d.d[t][p][i];


    for (int t = 1; t <= d.T; t++) 
        for (int p = 1; p <= d.P; p++)
            for (int i = 1; i <= d.n; i++) if (t > 1) d.D[t][p] += d.d[t][p][i];


    for (int i = 0; i <= d.n + 1; i++)
        for (int j = 0; j <= d.n + 1; j++)
            d.c[i][j] = round(sqrt((d.x[i] - d.x[j])*(d.x[i] - d.x[j]) + (d.y[i] - d.y[j])*(d.y[i] - d.y[j])));

    d.a = std::vector<std::vector<double> > (d.c);


    arq.close();
}

int main(int argc, char *argv[]){
    data d;
    d.name = argv[1];   //instance name
            
    d.pathIn = "/X/";
    read(d);
    
    return 0;    
}
