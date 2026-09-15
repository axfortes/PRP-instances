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
using chrono::system_clock;
duration<int,ratio<1> > segundo(1); //duraçao de um segundo
duration<int,ratio<60> > minuto(1); //duraçao de um minuto
duration<int,ratio<300> > cincomin(1); //duraçao de um minuto
duration<int,ratio<1200> > vinte(1); //duraçao de uma hora
duration<int,ratio<1800> > meiahora(1); //duraçao de uma hora
duration<int,ratio<3600> > hora(1); //duraçao de uma hora
duration<int,ratio<7200> > duashoras(1); //duraçao de duas hora

//armazenando o sucesso de cada busca
std::vector<int> iers (9,0); //inter-route success
std::vector<int> iars (8,0); //intra-route success
std::vector<int> sini (5,0); //initial VRP solutions
int totalOfPerturbations = 0;
int totalOfImprovements = 0;

//parametros calibrados pelo irace, config 14
int maxIter = 300;  //numero maximo de iterações do laco principal
int maxIterILS = 15;//numero maximo de iteracoes do ILS
int maxPertP = 3; //numero de perturbacoes do tipo Swap(K,K)
int maxWeight = 5;  //peso maximo adotado na perturbacao dos alpha
double alpha = 0.1; //distancia máxima da solucao corrente

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

std::random_device rd;     //Will be used to obtain a seed for the random number engine
std::mt19937_64 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

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

    //informacoes da frota
    int V;                              //numero de veiculos
    std::vector<int> Q;                 //capacidade maxima de transporte de determinado veiculo
    std::vector<int> e;                 //custo de setup/ativação de cada veiculo
    int H;                              //horizonte de tempo para entrega, de 12 horas (6h - 18h)
    int Qmax;
    int Qmin;
    int maxLoad;
    std::vector<int> accumLoad;

    //informacoes gerais
    std::vector<std::vector<double> > c; //custo de transporte de i para j
    std::vector<double> x;               //coordenada x do no i
    std::vector<double> y;               //coordenada y do no i
    std::vector<std::vector<double> > h; //custo de estocagem do produto p no no i
    std::vector<std::vector<int> > U;    //limite superior de estocagem no item de cada planta
    std::vector<std::vector<int> > I0;   //inventario inicial em cada no, deve satisfazer a demanda do periodo    

    std::vector<std::vector<int> > N;           //closest neighbors
    std::vector<boost::dynamic_bitset<> >Nb;    //binary set of closest neighbors
    std::vector<double> Delta;                  //estimated visitation cost    
        
    int ntMax;                                      //numero maximo de vizinhancas

    std::vector<std::vector<int> > M_;              //será usado nos cortes no-Good de integralidade para pkt

    string name;
    string pathIn;
    string pathOut;

    //parametros para calibrar
    int maxTime;
    int maxIter;

    double minGAP;
    double maxGAP;

    duration<double> totTime;
}data;

//model information
typedef struct {
    IloNum ub;          //upper bound
    IloNum lb;          //lower bound
    IloNum gap;
    IloTimer *crono;
    IloCplex cplx;
    IloModel model;
    //binary variables
    IloNumVarArray g; 	//defined in k,t
    IloNumArray g_;
    IloNumVarArray x; 	//defined in (i,j), v, t
    IloNumArray x_;
    IloNumVarArray y; 	//defined in j,k,t
    IloNumArray y_;
    IloNumVarArray z; 	//defined in i,v,t
    IloNumArray z_;
    //integer variables
    IloNumVarArray b; 	//defined in i,k,t
    IloNumArray b_;
    IloNumVarArray f; 	//defined in
    IloNumArray f_;
    IloNumVarArray I; 	//defined in i,k,t
    IloNumArray I_;
    IloNumVarArray o;	//defined in i,k,t
    IloNumArray o_;
    IloNumVarArray p; 	//defined in k,t
    IloNumArray p_;
    IloNumVarArray q;	//defined in i,k,v,t
    IloNumArray q_;
    IloNumVarArray w;	//defined in (i,j), v, t
    IloNumArray w_;

    IloRangeArray constraints;
    IloRangeArray r01Prod;
    IloRangeArray r02InvBalFac;
    IloRangeArray r03InvBalCust;
    IloRangeArray r04InvLim;
    IloRangeArray r05MaxLoadPerVehicle;
    IloRangeArray r06MaxLoadDeliv;
    IloRangeArray r07VMaxNumberOfVisit;

    IloObjective objective;
} cpx;

//routing information class
typedef struct {
    std::vector<std::vector<int> > q;                       //aggregated load for each client per pediod period
    std::vector<std::vector<boost::dynamic_bitset<> > > z;  // if some customer i is visited by vehicle v at period t
    std::vector<int> V;                                     //maximum number of necessary vehicles per period
    std::vector<std::vector<bool> > g;                      //which vehicles are used at each period
    std::vector<int> D;                                     //total load
} routing;

typedef struct {
    double f;
    std::vector<std::vector<int> > p;                   //producao unitaria
    std::vector<std::vector<int> > y;                   //setup de producao
    std::vector<std::vector<std::vector<int> > > o;     //ofertas periodicas
    std::vector<std::vector<std::vector<int> > > I;     //nivel de estoque
    std::vector<std::vector<std::vector<int> > > b;     //nivel de postergaçao
    std::vector<boost::dynamic_bitset<> > z;            //visitacao    
} prodinv;

//individual route information
typedef struct{
    double rCost;                       //total cost
    int rLoad;                          //total load
    std::deque<int> rOrder;             //visitation order
    boost::dynamic_bitset<> rComp;    //composition
    int rFST;                           //first visited
    int rLST;                           //last visited
    int rTTT;                           //total travel time
} route;

//vehicle routing solution/subproblem (vrs)
typedef struct{
    std::vector<route> R;               //routes
    double c;                           //total solution cost
    std::vector<int> cpp;               //clients per partition
    std::vector<std::pair<int,int> >ps; //predecessor, successor
    bool anyRoute;
} vrs;

typedef struct{
    cpx ppp;
    std::vector<vrs> vrp;
    routing r;
    prodinv p;
    std::vector<double> fo;
    double bestTime;    
} solprp;

//clarke-wright saving
typedef struct{
    std::pair<int,int> stop;    //indicates the nodes i (first) and j (second) that compose the saving
    double save;                //indicates the total SAVED gathering i and j
} cws;

/*VLNS structs */
//node class for very large-scale neighborhoods exchanges
typedef struct{
    int pred;               //predecessor
    int succ;               //sucessor
    int indx;               //position that wil, be replaced, inserted of erased
} pos;

// //arc class for very large-scale neighborhoods exchanges
typedef struct{
    int i;      //arc start
    int j;      //arc end
    int w;      //weight
    int t;      //type (2 - allocation), (3 - remotion), (4 - complementary) or (1 - substuitution)
    int n;
} arc;

typedef struct{
    int nArcs;               // contador de arcos
    int nNods;               // contador de nós
    double maxCost;          // arco com o maior custo
    int src,snk;             // origem e destino do caminho
    std::vector<arc> arcs;   // lista de arcos do grafo de melhoria
    std::vector<pos> nPos;   // para cada no i cliente, regista qual a posição ele ocupará na nova partição
    std::vector<std::vector<int> >ftd;   // indice dos arcos que deixam o no i, Fecho Transitivo Direto
} net;

typedef struct{
    int k;     //product with load changed
    int p;     //quantity transferred
    int n;     //movement number id
} mov;

//function heads
//auxiliares
void help();

//data section
void read(data &d, string sPath, string sName);
void createNeighborhood(std::vector<std::vector<double> > c, int n, std::vector<std::vector<int> > &N, int nSize, std::vector<boost::dynamic_bitset<> > &Nb);
void selection_sort(std::vector<double> &cost,  std::vector<int> &Ni, int tam);

//milp production-inventory-designation planning (Production-Planning Model)
void pppModel(data d, cpx &c);
void getInfoPPP (data d, cpx &c);
void pppRecovInfo(data d, cpx c, routing &r, prodinv &pid, std::vector<double> &f);
bool buildAndSolvePPP(data d, cpx &ppp);
void updateDelta(data d, std::vector<vrs> vrp, prodinv p, cpx &ppp);
bool fixSetups(data d, bool fix, prodinv p, cpx &ppp);

//copy information methods
routing cpRout(routing r0);
prodinv cpProd(prodinv p0);
solprp cpSol(solprp s0);
prodinv newProd(data d);
routing newRout(data d);

//checking solutions costs, exporting and printing
void export2Sol(data d, solprp s);
std::vector<double> checkSolCost(data d, prodinv p, std::vector<vrs> vrp);


/*Routing Methods*/
double checkRoutCost(data d, std::deque<int> rota);
int checkRoutTTT(data d, std::deque<int> rota);
int checkLoad(std::deque<int> r, std::vector<int> q);

route newRoute(int n);
route newRoute(double cost, int ttt, int load, std::deque<int> order, boost::dynamic_bitset<> bit);
route cpyRoute(route r0);
void printRoute(data d, std::vector<int> q, route r, int v);
void printRoute(data d, route r, int v);

vrs newVRS(int n, int V);
vrs cpyVRS(vrs s0);
void printVRS(data d, routing r, vrs s, int T);
bool checkVRSCost (data d, routing r, vrs s, int t, string name);
void updatePredSucc(int n, std::deque<int> rOrder, std::vector<std::pair<int,int> > &ps);

//GREEDY initialçl solution
vrs mgr(data d, routing &r, int T, bool &failed);
void greedyRoute(data d, std::vector<int> q, int &v, boost::dynamic_bitset<> &N, vrs &s);
vrs lgr(data d, routing &r, int T, bool &failed);

/*CLARKE-WRIGHT SAVING HEURISTIC methods*/
void createCWSavingList(data d, routing r, std::vector<cws> &savings, int T);
void sortCW(std::vector<cws> &savings);
void clearSavListCW (std::vector<cws> &savings, int i);
vrs pcw(data d, routing &r, int T, bool &failed);
vrs scw(data d, routing &r, int T, bool &failed);

vrs S0(data d, int t, routing &r, bool &failed);

//LOCAL SEARCH ROUTING METHODS
//INTER-ROUTE HEURISTICS
bool soloCustAlloc (data d, vrs &vrsol, routing &r, int t, int v1);
//vlns
arc newArc(int i, int j, int w, int t, int n);
pos newNodePos(int p, int s, int i);
net improvementGraph(data d, routing r, vrs s, int T);
std::vector<int> labelCorrecting(data d, net g);
void recoverCycle(int i, std::vector<int> pred, std::vector<int> minArcList, std::vector<int> &cycle, boost::dynamic_bitset<> &cBit);
bool cycleValidation(data d, std::vector<int> cpp, net g, std::vector<int> &cycle);
void buildSolution(data d, std::vector<int> &cycle, net g, routing &r, vrs &sol, int T);
net clearGraph();
bool vlns(data d, routing &r, int T, vrs &s);

bool shift(data d, vrs &s, routing &r, int t, int K);
bool cross(data d, vrs &s, routing &r, int t);
bool shift20(data d, vrs &s, routing &r, int t);
bool swap21(data d, vrs &s, routing &r, int t);
bool swap22(data d, vrs &s, routing &r, int t);
bool interRouteOpt (data d, routing &r, vrs &s, int t);

//INTRA-ROUTE HEURISTICS
double OnePnt(data d, std::vector<std::pair<int,int> > &ps, route &r);
double TwoPnt(data d, std::vector<std::pair<int,int> > &ps, route &r);
double ThrPnt(data d, std::vector<std::pair<int,int> > &ps, route &r, int v);
double TwoOpt(data d, std::vector<std::pair<int,int> > &ps, route &r);
double OrOpt(data d, std::vector<std::pair<int,int> > &ps, route &r, int Or, int v);
bool intraRouteOpt(data d, std::vector<bool> mudou, vrs &s);

//Vehicle Routing Problem Optmization/Improvement/Local Searches
bool BuildAndOptimizeVRS(data d, std::vector<vrs> &vrsol, routing &r, double &totVRPcost);

//perturbacao do PPP
mov newMove(int i, int k, int o, int n);
std::vector<std::vector<std::vector<mov> > > genTransfProdMoves (data d, prodinv p, std::vector<vrs> vrp, int &n);
bool modifyProdPlan(data d, prodinv &p, std::vector<vrs> vrp, cpx &ppp);
bool correctBoundsProdVars(data d, cpx &ppp);

//Tactical ILS
double tils(data d);

//Erro de leitura de dados
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

    //std::cout << "\nV = " << d.V << "\nP = " << d.P;

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

    //informacoes sobre producao e fabrica
    d.u = std::vector<double> (d.P + 1,0); //custo de producao unitario
    d.C = std::vector<int> (d.P + 1,0);       //capacidade de cada linha de produção
    d.l = std::vector<double> (d.P + 1,0); //custo de preparação de produção (setup)    

    //informacoes sobre os veiculos
    d.Q = std::vector<int> (d.V + 1, 0);        //carga maxima carregada
    d.e = std::vector<int> (d.V + 1, 0);        //custo de ativação

    //custos e limites de estocagem e custo de atraso
    d.h = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n + 1,0.0));
    d.U = std::vector<std::vector<int> > (d.P + 1, std::vector<int> (d.n + 1, 0));
    d.B = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n + 1,0.0));

    //informacoes sobre os veiculos
    d.Qmax = 0;
    d.Qmin = 99999999;
    d.maxLoad = 0;
    d.accumLoad = std::vector<int> (d.V + 1, 0);    
    for (int v = 1; v <= d.V; v++){
        arq >> d.Q[v];
        arq >> d.e[v];
        if (d.Q[v] > d.Qmax) d.Qmax = d.Q[v];
        if (d.Q[v] < d.Qmin) d.Qmin = d.Q[v];
        d.maxLoad += d.Q[v];
        d.accumLoad[v] += (d.accumLoad[v-1] + d.Q[v]);        
        //std::cout << "\nv: " << v << " - Q: " << d.Q[v] << " - e: " << d.e[v] << " -> Accumulated: " << d.accumLoad[v];
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

    for (int i = 0; i <= d.n + 1; i++)
        for (int j = 0; j <= d.n + 1; j++)
            d.c[i][j] = round(sqrt((d.x[i] - d.x[j])*(d.x[i] - d.x[j]) + (d.y[i] - d.y[j])*(d.y[i] - d.y[j])));

    d.a = std::vector<std::vector<double> > (d.c);

    //creating Neighborhoods of customers i
    int sizeNi = 10;

    d.N = std::vector<std::vector<int > > (d.n);
    d.Nb = std::vector<boost::dynamic_bitset<> > (d.n, boost::dynamic_bitset<> (d.n,0));
    createNeighborhood(d.c, d.n, d.N, sizeNi, d.Nb);
    //Neighborhoods created

    d.Delta = std::vector<double> (d.n + 1, 0);
    for (int i = 1; i <= d.n; i++){
        int j = d.N[i-1][0];
        int k = d.N[i-1][1];
        d.Delta[i] = std::min(2*d.c[0][i],0.5*(d.c[i][j] + d.c[k][i]));
        //std::cout << "\nDelta[" << i << "] = " << d.Delta[i];
    }
    //std::cout << std::endl; 

    d.M_ = std::vector<std::vector<int> > (d.T + 1, std::vector<int> (d.P + 1, 0));
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            int dMaxK = 0;
            for (int e = t; e <= d.T; e++)
                for (int i = 1; i <= d.n; i++) dMaxK += d.d[e][k][i];
            int M = std::min(d.C[k],dMaxK);
            d.M_[t][k] = M;
        }
    }

    d.ntMax = 5;

    arq.close();
}
void createNeighborhood(std::vector<std::vector<double> > c, int n, std::vector<std::vector<int > > &N, int nSize, std::vector<boost::dynamic_bitset<> >&Nb){
    for (int i = 1; i <= n; i++){
        std::vector<double> cost;
        std::vector<int> Ni;

        for (int j = 1; j <= n; j++) {
            if (j != i) {
                cost.push_back(c[i][j]);
                Ni.push_back(j);
            }
        }
        selection_sort(cost, Ni, Ni.size());

        for (int k = 0; k < (nSize-1); k++) N[i-1].push_back(Ni[k]);

        N[i-1].push_back(i);
        std::sort (N[i-1].begin(), N[i-1].begin()+nSize);

        for(int j = 0; j < N[i-1].size(); j++){
            boost::dynamic_bitset<> t (n,1);
            t <<= (N[i-1][j]-1);
            Nb[i-1] = (Nb[i-1] | t);
        }
    }
}
void selection_sort(std::vector<double> &cost,  std::vector<int> &Ni, int tam){
    int i, j, imin, iaux;
    double daux;
    for (i = 0; i < (tam-1); i++){  // para todas as posições
        imin = i;
        for (j = (i+1); j < tam; j++)
            if(cost[j] < cost[imin])   //avalio o valor armazenado
                imin = j;

        if (i != imin){
            daux = cost[i];
            cost[i] = cost[imin];
            cost[imin] = daux;

            iaux = Ni[i];
            Ni[i] = Ni[imin];
            Ni[imin] = iaux;
        }
    }
}

//milp production-inventory-designation planning (Production-Planning Model)
void pppModel(data d, cpx &c){
    //std::cout << "\n\tConstruindo modelo";
    std::uniform_int_distribution<> W(1,maxWeight);
    
    //variable vector's size
    int sY = d.P*d.T;                 //p tem o mesmo dominio        
    int sQ = d.n*d.P*d.V*d.T;            //n para representar a carga entregue aos clientes
    int sI = (d.n + 1)*d.P*d.T;          //n+1 para representar a planta e clientes
    int sB = d.n*d.P*d.T;                //n para representar os clientes

    IloEnv env = c.cplx.getEnv();

    // begin variables' declaration
    //binaries   
    c.y = IloNumVarArray(env, sY, 0.0, 1.0, ILOINT);
    c.y_= IloNumArray(env,sY);

    //nonnegative
    c.q = IloNumVarArray(env, sQ, 0.0, IloInfinity, ILOFLOAT);
    c.q_= IloNumArray(env, sQ);
    
    c.p = IloNumVarArray(env, sY, 0.0, IloInfinity, ILOFLOAT);
    c.p_= IloNumArray(env,sY);

    c.I = IloNumVarArray(env, sI, 0.0, IloInfinity, ILOFLOAT);
    c.I_= IloNumArray(env,sI);

    c.b = IloNumVarArray(env, sB, 0.0, IloInfinity, ILOFLOAT);
    c.b_= IloNumArray(env,sB);
    //end variables' declaration

    //turning to 0 de backorder at last period
    for (int k = 1; k <= d.P; k++)
        for (int i = 1; i <= d.n; i++)
            c.b[d.P*d.n*(d.T-1) + d.n*(k-1) + (i-1)].setUB(0);
    
    //inventory limits
    for (int t = 1; t <= d.T; t++)
        for (int k = 1; k <= d.P; k++)
            for (int i = 0; i <= d.n; i++)
                c.I[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)].setUB(d.U[k][i]);
                
    //delivery lots limits    
    for (int t = 1; t <= d.T; t++){
        for (int v = 1; v <= d.V; v++){
            for (int k = 1; k <= d.P; k++){
                for (int i = 1; i <= d.n; i++){
                    int rDem = 0;
                    int M = std::min(d.U[k][i],d.Q[v]);
                    for (int e = t; e <= d.T; e++) rDem += d.d[e][k][i];
                    M = std::min(M,rDem);                    
                    c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)].setUB(M);                                    
                }
            }
        }
    }
            
    //variables' name
    char name[50];

    //begin objective function declaration
    IloExpr ObjFunc(env);
    for (int t = 1; t <= d.T; t++){
        //production, setup, storage and backorder components costs
        for (int k = 1; k <= d.P; k++){
            //production and setup
            ObjFunc += d.u[k] * c.p[d.P*(t-1) + (k-1)];
            sprintf(name,"p[%d][%d]",t,k);
            c.p[d.P*(t-1) + (k-1)].setName(name);

            ObjFunc += d.l[k] * c.y[d.P*(t-1) + (k-1)];
            sprintf(name,"y[%d][%d]",t,k);
            c.y[d.P*(t-1) + (k-1)].setName(name);

            //inventory and backorder
            for (int i = 0; i <= d.n; i++){
                ObjFunc += d.h[k][i] * c.I[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)];
                sprintf(name,"I[%d][%d][%d]",t,k,i);
                c.I[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)].setName(name);

                if (i > 0){
                    ObjFunc += d.B[k][i] * c.b[d.P*d.n*(t-1) + d.n*(k-1) + (i-1)];
                    sprintf(name,"b[%d][%d][%d]",t,k,i);
                    c.b[d.P*d.n*(t-1) + d.n*(k-1) + (i-1)].setName(name);
                    
                    for (int v = 1; v <= d.V; v++){
                        ObjFunc += W(gen)*int(std::ceil(d.Delta[i] + d.e[v]/d.Q[v] + std::min(d.u[k] + d.h[k][0], d.u[k] + d.h[k][i]))) * c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];
                        sprintf(name,"q[%d][%d][%d][%d]",t,v,k,i);
                        c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)].setName(name);
                    }
                }
            }
        }
    }

    c.objective = IloAdd(c.model,IloMinimize(env, ObjFunc));
    //end objective function declaration

    //begin constraints' declaration
    //c.constraints = IloRangeArray(env);
    c.r01Prod = IloRangeArray(env);
    c.r02InvBalFac = IloRangeArray(env);
    c.r03InvBalCust = IloRangeArray(env);    
    c.r05MaxLoadPerVehicle = IloRangeArray(env);
    //c.r06MaxLoadDeliv = IloRangeArray(env);
  
    //PRODUCTION AND INVENTORY
    //maximum production lot size (3.2)
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            IloExpr rProd(env);
            IloRange ctrProd;
            sprintf(name,"rProd(%d,%d)", t,k);

            rProd += c.p[d.P*(t-1) + (k-1)];
            rProd -= d.M_[t][k]*c.y[d.P*(t-1) + (k-1)];

            ctrProd = (rProd <= 0);
            ctrProd.setName(name);
            //c.constraints.add(ctrProd);
            c.r01Prod.add(ctrProd);
            rProd.end();
        }
    }

    //inventory balance at factory (3.3)
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            IloExpr rInvFac(env);
            IloRange ctrInvFac;

            rInvFac += c.I[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (0)];
            if (t == 1)
                rInvFac -= d.I0[k][0];
            else if (t > 1)
                rInvFac -= c.I[d.P*(d.n + 1)*(t-2) + (d.n + 1)*(k-1) + (0)];

            rInvFac -= c.p[d.P*(t-1) + (k-1)];

            for (int i = 1; i <= d.n; i++)
                for (int v = 1; v <= d.V; v++)
                    rInvFac += c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];

            sprintf(name,"rInvBalFac(%d,%d,%d)", t,k,0);
            ctrInvFac = (0 <= rInvFac <= 0);
            ctrInvFac.setName(name);
            //c.constraints.add(ctrInvFac);
            c.r02InvBalFac.add(ctrInvFac);
            rInvFac.end();
        }
    }

    //inventory balance at customers (3.4)
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            for (int i = 1; i <= d.n; i++){
                IloExpr rInvCust(env);
                IloRange ctrInvCust;
                sprintf(name,"rInvBalCust(%d,%d,%d)", t,k,i);
                rInvCust += c.I[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)];
                if (t == 1)
                    rInvCust -= d.I0[k][i];
                else if (t > 1)
                    rInvCust -= c.I[d.P*(d.n + 1)*(t-2) + (d.n + 1)*(k-1) + (i)];                    

                for (int v = 1; v <= d.V; v++)
                    rInvCust -= c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];

                rInvCust += d.d[t][k][i];

                rInvCust -= c.b[d.P*d.n*(t-1) + d.n*(k-1) + (i-1)];
                if (t > 1)  //backorder at period 0 is 0!
                    rInvCust += c.b[d.P*d.n*(t-2) + d.n*(k-1) + (i-1)];

                ctrInvCust = (0 <= rInvCust <= 0);
                ctrInvCust.setName(name);
                //c.constraints.add(ctrInvCust);
                c.r03InvBalCust.add(ctrInvCust);
                rInvCust.end();
            }
        }
    }

    //ROUTING
    //maximum load carried by each vehicle v at each period t (3.08)
    for (int t = 1; t <= d.T; t++){
        for (int v = 1; v <= d.V; v++){
            IloExpr rMaxLoadPerVehicle(env);
            IloRange ctrMaxLoadPerVehicle;
            sprintf(name,"rMaxLoadPerVehicle(%d,%d)",t,v);

            for (int k = 1; k <= d.P; k++)
                for (int i = 1; i <= d.n; i++)
                    rMaxLoadPerVehicle += c.q[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];            

            ctrMaxLoadPerVehicle = (rMaxLoadPerVehicle <= d.Q[v]);
            ctrMaxLoadPerVehicle.setName(name);
            //c.constraints.add(ctrMaxLoadPerVehicle);
            c.r05MaxLoadPerVehicle.add(ctrMaxLoadPerVehicle);
            rMaxLoadPerVehicle.end();
        }
    }
    
    //c.model.add(c.constraints);
    c.model.add(c.r01Prod);
    c.model.add(c.r02InvBalFac);
    c.model.add(c.r03InvBalCust);    
    c.model.add(c.r05MaxLoadPerVehicle);
    //c.model.add(c.r06MaxLoadDeliv);
    
    //c.cplx.exportModel("ModelOriginal.lp");
}
void getInfoPPP (data d, cpx &ppp){
    //std::cout << "\n\n=================================" << RED << "\n\tLSPDS solution";
    try{
        ppp.ub = (double) ppp.cplx.getObjValue();
        ppp.lb = ppp.ub;        
        ppp.cplx.getValues(ppp.p_,ppp.p);
        ppp.cplx.getValues(ppp.I_,ppp.I);
        ppp.cplx.getValues(ppp.b_,ppp.b);
        ppp.cplx.getValues(ppp.q_,ppp.q);
        ppp.cplx.getValues(ppp.y_,ppp.y);        
    }
    catch(IloException& ex){
        std::cerr << RED << "\nError: " << ex << std::endl;
        exit(1);
    }
}
void pppRecovInfo(data d, cpx ppp, routing &r, prodinv &p, std::vector<double> &fo){    
    //no modelo, a somatória do que é entregue não pode ser maior que a somatória das capacidades
    //std::cout << normal <<  "\n\t\tProd/Inv/Rout >> " << normal;        
    r = newRout(d);
    p = newProd(d);
    fo = std::vector<double> (4, 0.0);
    
    //accumulated load of each customer per period
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            p.y[t][k] = int(ppp.y_[d.P*(t-1) + (k-1)]);
            p.p[t][k] = int(ppp.p_[d.P*(t-1) + (k-1)]);

            fo[1] += (d.l[k] * ppp.y_[d.P*(t-1) + (k-1)]);
            fo[1] += (d.u[k] * ppp.p_[d.P*(t-1) + (k-1)]);
                            
            for (int i = 0; i <= d.n; i++){
                p.I[t][k][i] = double(ppp.I_[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)]);
                fo[2] += (d.h[k][i] * ppp.I_[d.P*(d.n + 1)*(t-1) + (d.n + 1)*(k-1) + (i)]);            
                if (i > 0){
                    p.b[t][k][i] = double(ppp.b_[d.P*d.n*(t-1) + d.n*(k-1) + (i-1)]);
                    fo[2] += (d.B[k][i] * ppp.b_[d.P*d.n*(t-1) + d.n*(k-1) + (i-1)]);
                    for (int v = 1; v <= d.V; v++){
                        r.q[t][i] += ppp.q_[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];
                        p.o[t][k][i] += ppp.q_[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];
                    }
                }
            }
        }
        
        //recovering the clients visited
        for (int i = 1; i <= d.n; i++) 
            if (r.q[t][i] > 0) p.z[t].set(i-1);
    }

    fo[0] = fo[1] + fo[2];
    p.f = fo[0];
    //std::cout << "\n\tppp = " << pppCost << " = ";
    //std::cout << fo[1] << " (prod) + " << fo[2] << " (inv)";    
}
bool buildAndSolvePPP(data d, cpx &ppp){
    try{
        IloEnv env;
        IloModel model(env);
        IloCplex cplx(model);
        IloTimer crono(env);

        cplx.setParam(IloCplex::TiLim,600);
        cplx.setParam(IloCplex::EpGap,0.05);
        cplx.setWarning(cplx.getEnv().getNullStream());
        cplx.setOut(cplx.getEnv().getNullStream());

        ppp.cplx = cplx;
        ppp.model = model;
        ppp.crono = &crono;
        pppModel(d,ppp);

        ppp.cplx.solve();
        ppp.crono->stop();
        getInfoPPP(d,ppp);        
        ppp.cplx.setParam(IloCplex::TiLim,120);
        ppp.cplx.setParam(IloCplex::EpGap,0.025);
        return true;        
    }
    catch(IloException& ex){
        std::cout << d.name << " infeasible";
        exit(1);
        return false;
        //std::cerr << RED << "\nError: " << ex << std::endl;
    }	
}

//fix and unfixing setups
bool fixSetups(data d, bool fix, prodinv p, cpx &ppp){
    IloNum um = 1;
    IloNum zr = 0;
    if (fix){
        bool fixed = false;
        //std::cout << "\nFixing setups";
        try {
            for (int t = 1; t <= d.T; t++){
                for (int k = 1; k <= d.P; k++){
                    if (p.y[t][k] == 1)
                        ppp.y[d.P*(t-1) + (k-1)].setBounds(um,um);
                    else if (p.y[t][k] == 0)
                        ppp.y[d.P*(t-1) + (k-1)].setBounds(zr,zr);
                }
            }
            return true;
        }catch(IloException& ex){ return false; }
    }else {
        bool unfixed = false;
        //std::cout << "\nUnfixing setups";
        try {
            for (int t = 1; t <= d.T; t++)
                for (int k = 1; k <= d.P; k++)
                    ppp.y[d.P*(t-1) + (k-1)].setBounds(zr,um);
            return true;
        }catch(IloException& ex){ return false; }
    }
}
void updateDelta(data d, std::vector<vrs> vrp, prodinv p, cpx &ppp){
    std::uniform_int_distribution<> W(1,maxWeight);
    
    IloNumArray delta (ppp.cplx.getEnv(), d.n*d.P*d.V*d.T);
        
    for (int t = 1; t <= d.T; t++){
        for (int v = 1; v <= d.V; v++){
            for (int i = 1; i <= d.n; i++){
                //se o cliente i eh visitado pela rova v
                if (vrp[t].R[v].rComp.test(i-1)){
                    for (int k = 1; k <= d.P; k++){
                        delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil(d.e[v]/d.Q[v]));
                        int h = vrp[t].ps[i].first;
                        int j = vrp[t].ps[i].second;
                        delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil((d.c[h][i] + d.c[i][j])/vrp[t].R[v].rCost));
                        if (t > 1){
                            if (p.I[t][k][i] - p.I[t-1][k][i] > 0) delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil(d.h[k][i]/(p.I[t][k][i] - p.I[t-1][k][i])));
                            if (p.p[t][k] - p.p[t-1][k] > 0) delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil(d.u[k]/(p.p[t][k] - p.p[t-1][k])));
                        }else{
                            if (p.I[t][k][i] > 0) delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil(d.h[k][i]/p.I[t][k][i]));
                            if (p.p[t][k] > 0) delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] += int(std::ceil(d.u[k]/p.p[t][k]));
                        }                        
                        delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] = W(gen) * delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)];
                    }
                }
                else
                    for (int k = 1; k <= d.P; k++) 
                        delta[d.V*d.P*d.n*(t-1) + d.P*d.n*(v-1) + d.n*(k-1) + (i-1)] = W(gen) * int(std::ceil(d.Delta[i] + d.e[v]/d.Q[v] + std::min(d.u[k] + d.h[k][0], d.u[k] + d.h[k][i])));                                        
            }            
        }        
    }
    
    ppp.objective.setLinearCoefs(ppp.q,delta);
    //ppp.cplx.exportModel("ModelwithDeltas.lp");
    //exit(1);
    //return delta; 
}

//copy information methods
routing cpRout(routing r0){
    routing r1;
    r1.q = std::vector<std::vector<int> > (r0.q);
    r1.z = std::vector<std::vector<boost::dynamic_bitset<> > > (r0.z);
    r1.V = std::vector<int> (r0.V);
    r1.g = std::vector<std::vector<bool> > (r0.g);
    r1.D = std::vector<int> (r0.D);
    return r1;
}
prodinv cpProd(prodinv p0){
    prodinv p1;
    p1.f = p0.f;
    p1.p = std::vector<std::vector<int> > (p0.p);
    p1.y = std::vector<std::vector<int> > (p0.y);                   //setup de producao
    p1.o = std::vector<std::vector<std::vector<int> > > (p0.o);     //ofertas periodicas    
    p1.I = std::vector<std::vector<std::vector<int> > > (p0.I);     //nivel de estoque
    p1.b = std::vector<std::vector<std::vector<int> > > (p0.b);     //nivel de postergaçao    
    p1.z = std::vector<boost::dynamic_bitset<> > (p0.z);            //visitacao
    return p1;
}
solprp cpSol(solprp s0){
    solprp s1;
    s1.ppp = s0.ppp;
    s1.vrp = std::vector<vrs> (s0.vrp);
    s1.r = cpRout(s0.r);
    s1.p = cpProd(s0.p);
    s1.fo = std::vector<double> (s0.fo);    
    return s1;
}
prodinv newProd(data d){
    prodinv p;
    p.o = std::vector<std::vector<std::vector<int> > > (d.T + 1, std::vector<std::vector <int> > (d.P + 1, std::vector<int>(d.n + 1, 0)));    
    p.I = std::vector<std::vector<std::vector<int> > > (d.T + 1, std::vector<std::vector <int> > (d.P + 1, std::vector<int>(d.n + 1, 0)));
    p.b = std::vector<std::vector<std::vector<int> > > (d.T + 1, std::vector<std::vector <int> > (d.P + 1, std::vector<int>(d.n + 1, 0)));    
    p.p = std::vector<std::vector<int> > (d.T + 1, std::vector<int> (d.P + 1,0));
    p.y = std::vector<std::vector<int> > (d.T + 1, std::vector<int> (d.P + 1,0));
    p.z = std::vector<boost::dynamic_bitset<> > (d.T + 1, boost::dynamic_bitset<>(d.n, 0));
    return p;    
}
routing newRout(data d){
    routing r;
    r.q = std::vector<std::vector<int> > (d.T + 1, std::vector<int> (d.n + 1,0));
    r.V = std::vector<int> (d.T + 1,0);
    r.g = std::vector<std::vector<bool> > (d.T + 1,std::vector<bool> (d.V + 1, false));
    r.z = std::vector<std::vector<boost::dynamic_bitset<> > > (d.T + 1, std::vector<boost::dynamic_bitset<> > (d.V + 1, boost::dynamic_bitset<>(d.n, 0)));    
    return r;    
}
//end of Copy methods

//checking solutions costs, exporting and printing
void export2Sol(data d, solprp s){
    std::vector<vrs> vrp (s.vrp);
    routing r = cpRout(s.r);
    prodinv p = cpProd(s.p);

    string textRemove = ".dat";
    string nName = d.name;
    std::size_t found = nName.rfind(textRemove);
    if (found!=std::string::npos) nName.replace(found,nName.length(),"");

    ofstream arq;

    string name = nName + ".out";
    string outputfile = d.pathOut + name;
    //std::cout << "\n\toutputfile: " << outputfile;

    //convertendo string para char
    char *cName = new char[outputfile.length()+1];
    memcpy(cName, outputfile.c_str(), outputfile.length() + 1);


    arq.open(cName);

    if (!arq.is_open()) help();

    arq << "\t" << d.name << "\n\trunning = " << d.totTime.count() << "\n\tbest = " << s.bestTime;
    
    arq << "\n\t TOTAL COST: " << s.fo[0] << " = " << s.fo[1] << " + " << s.fo[2] << " + " << s.fo[3] << std::endl;

    //accessing basic variables one by one
    for (int t = 1; t <= d.T; t++){
        arq << "\n\t =====================================================" << "\n\t t = " << t;
        arq << "\n\t Setup and Production ";
        for (int k = 1; k <= d.P; k++)
            if (p.y[t][k] > 0)
                arq << "\n\t\t k = " << k << ": y = " << p.y[t][k] << "(l=" << d.l[k] << ") | p = " << p.p[t][k] << "(u=" << d.u[k] << ")";

        arq << "\n\n\t Inventory, Backorder and Delivery";
        std::vector<bool> printed (d.n + 1, false);
        for (int i = 0; i <= d.n; i++){
            arq << "\n\t\t i = " << i;
            for (int k = 1; k <= d.P; k++){
                arq << "\n\t\t\tk = " << k << ": I = " << p.I[t][k][i] << "(h=" << d.h[k][i] << ") ";
                if (i > 0) arq << "\t| b = " << p.b[t][k][i] << "(B=" << d.B[k][i] << ")\t| o = " << p.o[t][k][i];
            }
        }

        //routing component
        arq << std::endl << "\n\n\t Vehicle Activation, Visitation and Route";
        arq << "\n\t ====================================================\n\t VRP[" << t << "] solution -> Cost = " << vrp[t].c;
        for (int v = 1; v <= d.V; v++){
            if (r.g[t][v]){
                arq << "\n\t\t> v = " << v << " - bs: " << vrp[t].R[v].rComp << "\n\t\t\tTTT = " << vrp[t].R[v].rTTT << "\n\t\t\tCost = " << d.e[v] << " + " << vrp[t].R[v].rCost-d.e[v] << " = " << vrp[t].R[v].rCost << "\n\t\t\tload/Q = " << vrp[t].R[v].rLoad << "/" << d.Q[v] << "\n\t\t\torder: ";
                for (int c = 0; c < vrp[t].R[v].rOrder.size(); c++) arq << vrp[t].R[v].rOrder[c] << " ";
                arq << "\n\t\t\ttot" << " >> ";
                for (int c = 0; c < vrp[t].R[v].rOrder.size(); c++) arq << r.q[t][vrp[t].R[v].rOrder[c]] << " ";
                arq << "\n\t\t\t(p,s)" << " >> ";
                for (int c = 0; c < vrp[t].R[v].rOrder.size(); c++) arq << "(" << vrp[t].ps[vrp[t].R[v].rOrder[c]].first << "," << vrp[t].ps[vrp[t].R[v].rOrder[c]].second << ") ";
            }
        }
        arq << "\n\t =====================================================";
    }

    arq << std::endl;
    arq.close();
}
std::vector<double> checkSolCost(data d, prodinv p, std::vector<vrs> vrp){
    std::vector<double> fo (4, 0.0);
    for (int t = 1; t <= d.T; t++){
        for (int k = 1; k <= d.P; k++){
            fo[1] += (d.l[k] * p.y[t][k] + d.u[k] * p.p[t][k]);
            for (int i = 0; i <= d.n; i++) {
                fo[2] += (d.h[k][i] * p.I[t][k][i]);
                if (i > 0) fo[2] += (d.B[k][i] * p.b[t][k][i]);
            }
        }
        fo[3] += vrp[t].c;
    }
    fo[0] = fo[1] + fo[2] + fo[3];
    return fo;
}

/*Routing Methods*/
double checkRoutCost(data d, std::deque<int> rota){
    double cost = 0.0;
    if (rota.size() > 1){
        for (int ind = 0; ind < rota.size(); ind++){
            if (ind == 0)
                cost += (d.c[0][rota[ind]] + d.c[rota[ind]][rota[ind+1]]);
            else if (ind == (rota.size()-1))
                cost += d.c[rota[ind]][d.n+1];
            else
                cost += d.c[rota[ind]][rota[ind+1]];
        }
    }else if (rota.size() == 1)
        cost += (d.c[0][rota[0]] + d.c[rota[0]][d.n+1]);
    else{
        //std::cout << "\n\trota vazia!";
        return 99999999;
    }
    return cost;
}
int checkRoutTTT(data d, std::deque<int> rota){
    int TTT = 0;
    if (rota.size() > 1){
        for (int ind = 0; ind < rota.size(); ind++){
            if (ind == 0)
                TTT += (d.a[0][rota[ind]] + d.a[rota[ind]][rota[ind + 1]]);
            else if (ind == (rota.size()-1))
                TTT += d.a[rota[ind]][0];
            else
                TTT += (d.a[rota[ind]][rota[ind + 1]]);

            TTT += d.s[rota[ind]];
        }
    }else if (rota.size() == 1)
        TTT += (d.a[0][rota[0]] + d.s[rota[0]] + d.c[rota[0]][0]);
    else{
        //std::cout << "\n\trota vazia!";
        return 99999999;
    }
    return TTT;
}
int checkLoad(std::deque<int> r, std::vector<int> q){
    int demand = 0;
    for (int pos = 0; pos < r.size(); pos++) demand += q[r[pos]];
    return demand;
}

route newRoute(int n){
    route rt;
    rt.rComp = boost::dynamic_bitset<> (n,0);
    rt.rOrder = std::deque<int> (0);
    rt.rCost = 0;
    rt.rLoad = 0;
    rt.rFST = 0;
    rt.rLST = 0;
    rt.rTTT = 0;
    return rt;
}
route newRoute(double cost, int ttt, int load, std::deque<int> order, boost::dynamic_bitset<> bit){
    route R1;
    R1.rCost = cost;
    R1.rTTT = ttt;
    R1.rLoad = load;
    R1.rOrder = std::deque<int> (order);
    R1.rFST = R1.rOrder.front();
    R1.rLST = R1.rOrder.back();
    R1.rComp = boost::dynamic_bitset<> (bit);
    return R1;
}
route cpyRoute(route r0){
    route rt;
    rt.rComp = boost::dynamic_bitset<> (r0.rComp);
    rt.rOrder = std::deque<int> (r0.rOrder);
    rt.rCost = r0.rCost;
    rt.rLoad = r0.rLoad;
    rt.rFST = r0.rFST;
    rt.rLST = r0.rLST;
    rt.rTTT = r0.rTTT;
    return rt;
}
void printRoute(data d, std::vector<int> q, route r, int v){
    std::cout << "\n\t v = " << v << " >> cost = " << r.rCost << " = " << d.e[v] << " + " << r.rCost - d.e[v] << " >> ttt = " << r.rTTT << "/" << d.H << " >> load/Q = " << r.rLoad << "/" << d.Q[v] << "\n\t order = ";
    for (int c = 0; c < r.rOrder.size(); c++) std::cout << r.rOrder[c] << " ";
    std::cout << "\n\t dem = ";
    for (int c = 0; c < r.rOrder.size(); c++) std::cout << q[r.rOrder[c]] << " ";
    std::cout << "\n\t |r| = " << r.rComp.count() << " >> " << r.rComp;
    std::cout << std::endl;
}
void printRoute(data d, route r, int v){
    std::cout << "\n\t v = " << v << " >> cost = " << r.rCost << " = " << d.e[v] << " + " << r.rCost - d.e[v] << " >> ttt = " << r.rTTT << "/" << d.H << " >> load/Q = " << r.rLoad << "/" << d.Q[v] << "\n\t order = ";
    for (int c = 0; c < r.rOrder.size(); c++) std::cout << r.rOrder[c] << " ";
    std::cout << "\t |r| = " << r.rComp.count() << " >> " << r.rComp;
    std::cout << std::endl;
}

vrs newVRS(int n, int V){
    vrs s;
    s.R = std::vector<route> (V + 1, newRoute(n));
    s.c = 0.0;
    s.cpp = std::vector<int> (n + 1, 0);
    s.ps = std::vector<std::pair<int,int> > (n + 2, std::pair<int,int> (0,0));
    s.anyRoute = false;
    return s;
}
vrs cpyVRS(vrs s0){
    vrs s1;
    s1.R = std::vector<route> (s0.R);
    s1.c = s0.c;
    s1.cpp = std::vector<int> (s0.cpp);
    s1.ps = std::vector<std::pair<int,int> > (s0.ps);
    s1.anyRoute = s0.anyRoute;
    return s1;
}
void printVRS(data d, routing r, vrs s, int T){
    std::cout << "\n\t====================================================\n\tVRP[" << T << "] solution -> Cost = " << s.c;
    for (int v = 1; v <= d.V; v++){
        if (r.g[T][v]){
            std::cout << "\n\t\t> v = " << v << " |r| = " << s.R[v].rComp.count() <<  " >> r: " << s.R[v].rComp << "\n\t\t\tCost = " << d.e[v] << " + " << s.R[v].rCost-d.e[v] << " = " << s.R[v].rCost << "\n\t\t\tTTT/H = " << s.R[v].rTTT << "/" << d.H << "\n\t\t\tload/Q = " << s.R[v].rLoad << "/" << d.Q[v] << "\n\t\t\torder: ";
            for (int c = 0; c < s.R[v].rOrder.size(); c++) std::cout << s.R[v].rOrder[c] << " ";
            std::cout << "\n\t\t\tvehicle: ";
            for (int c = 0; c < s.R[v].rOrder.size(); c++) std::cout << s.cpp[s.R[v].rOrder[c]] << " ";
            std::cout << "\n\t\t\tq" << " >> ";
            for (int c = 0; c < s.R[v].rOrder.size(); c++) std::cout << r.q[T][s.R[v].rOrder[c]] << " ";
        }
    }

    std::cout << "\n\t==================================================== ";
}
bool checkVRSCost (data d, routing r, vrs s, int t, string name){
    bool ok = true;
    double costVRS = 0;
    for (int v = 1; v <= d.V; v++){
        if (r.g[t][v]) costVRS += (d.e[v] + checkRoutCost(d,s.R[v].rOrder));
        else{
            if (checkRoutCost(d,s.R[v].rOrder) > 0 && checkRoutCost(d,s.R[v].rOrder) < 99999999){
                //std::cout << RED << "\n\tErro na conferência, custo deveria ser nulo e é = " << checkRoutCost(d,s.R[v].rOrder) << normal;
                //help();
                return false;
            }
        }
    }    
        
    if (costVRS != s.c)  {
        /*
        std::cout << RED << "\n\t\tAlgo errado na " << name << ", custo = " << s.c << " <> custo real = " << costVRS;

        for (int v = 1; v <= d.V; v++){
            if (r.g[t][v]){
                double cost = d.e[v] + checkRoutCost(d,s.R[v].rOrder);
                if (cost != s.R[v].rCost){
                    std::cout << RED << "\n\t\tdeu zica nos custos de roteamento [v=" << v << "]" << " - certo: " << cost << " <> " << " - errado: " << s.R[v].rCost << normal;
                    std::cout << red;
                    printRoute(d,s.R[v],v);
                    std::cout << normal;
                    for (int pos = 0; pos < s.R[v].rOrder.size(); pos++){
                        int i = s.R[v].rOrder[pos];
                        int pi = s.ps[i].first;
                        int si = s.ps[i].second;
                        std::cout << red << "\n\t\t pred = " << pi << " >> i = " << i << " >> succ = " << si;
                    }
                }
            }
        }
        help();
        //*/
        //std::cout << normal;
        return false;
    }
    //std::cout << normal;
    return ok;
}
void updatePredSucc(int n, std::deque<int> rOrder, std::vector<std::pair<int,int> > &ps){
    for (int c = 0; c < rOrder.size(); c++){
        if (rOrder.size() > 1){
            if (c == 0)
                ps[rOrder[c]] = std::pair<int,int> (0, rOrder[c+1]);
            else if (c == rOrder.size()-1)
                ps[rOrder[c]] = std::pair<int,int> (rOrder[c-1],n + 1);
            else ps[rOrder[c]] = std::pair<int,int> (rOrder[c-1],rOrder[c+1]);
        }else
            ps[rOrder[c]] = std::pair<int,int> (0, n + 1);
    }
}


//GREEDY initial solution
//Maximum load/maxima vehicle capacity Greedy Route
vrs mgr(data d, routing &r, int T, bool &failed){
    //std::cout << "\n\t\t\tmgr t = " << T << " >> ";
    vrs s = newVRS(d.n,d.V);

    boost::dynamic_bitset<> N (d.n, 0);        
    for (int i = 1; i <= d.n; i++) if (r.q[T][i] > 0) N.set(i-1);
    
    if (N.none()) return s;
                
    std::vector<double> q;
    std::vector<int> ii;
    
    for (int j = 1; j <= d.n; j++){
        q.push_back(double(r.q[T][j]));
        ii.push_back(j);
    }
    selection_sort(q, ii, ii.size());
    
    //for (int j = 0; j < d.n; j++) std::cout << "\n\ti = " << ii[j] << " >> q = " << q[j];
    
    int loopCnt = 0;
    int maxLoop = 100;
    
    int v = d.V;
    while (N.any() && v > 0){
        int oldV = v;
        
        boost::dynamic_bitset<> R (d.n,0);
        int load = d.Q[v];
        int ttt = 0;
        int h = 0;
        std::deque<int> rOrder;
        
        for (int j = 0; j < d.n; j++){
            int i = ii[j];
            if(N.test(i-1)){
                if (r.q[T][i] <= load && ttt - d.a[h][0] + d.a[h][i] + d.s[i] + d.a[i][0] <= d.H){
                    R.set(i-1);
                    N.reset(i-1);
                    rOrder.push_back(i);
                    load = load - r.q[T][i];
                    ttt = ttt - d.a[h][0] + d.a[h][i] + d.s[i] + d.a[i][0];
                    h = i;
                    s.cpp[i] = v;
                }
            }
        }
        
        if (R.count() > 0){
            //salvando nova rota
            if (ttt != checkRoutTTT(d,rOrder)) std::cout << RED << "\n\tttt1 = " << ttt << " >> TTT = " << checkRoutTTT(d,rOrder);
            s.R[v] = newRoute(d.e[v] + checkRoutCost(d,rOrder), ttt, d.Q[v] - load, rOrder, R);
            s.c += s.R[v].rCost;
            if (s.R[v].rOrder.size() > 0) updatePredSucc(d.n,s.R[v].rOrder,s.ps);
            /*
            std::cout << red;
            printRoute(d,r.q[T],s.R[v],v);
            std::cout << normal;
            //*/
            R.reset();
            if (N.any()) v--;
            load = d.Q[v];
            ttt = 0;
            h = 0;
        }
        
        bool term = false;        
        if (oldV == v){
            loopCnt++;
            //std::cout << "\n\t N = " << N << " >> v = " << v << " >> Q = " << d.Q[v];
            //for(int i = 1; i <= d.n; i++) if (N.test(i-1)) std::cout << "\nq[" << i << "]=" << q[i];
            //std::cout << "\n\t loop = " << loopCnt;            
            if (loopCnt > maxLoop && v < d.V){
                v--;
                loopCnt = 0;
            }else
                term = true;
        }
        
        if (term) break;
    }
    
    if (N.any()){
        failed = true;
        s.anyRoute = false;
        return newVRS(d.n,d.V);
    }else{
        failed = false;
        //atualizando routing r
        r.V[T] = 0;
        for (int v = 1; v <= d.V; v++){
            //std::cout << BLUE << "\n\t\tv = " << v << " >> " << s.R[v].rComp << normal;
            if (s.R[v].rComp.any()){
                r.g[T][v] = true;
                r.z[T][v] = boost::dynamic_bitset<> (s.R[v].rComp);
                r.V[T]++;
                //std::cout << yellow << "\n\tg = " << r.g[T][v];                
            }else{
                r.g[T][v] = false;
                r.z[T][v] = boost::dynamic_bitset<> (d.n, 0);                
            }            
        }
        s.anyRoute = true;        
    }
        
    //if (!checkVRSCost(d,r,s,T,"GREEDY")) failed = true;

    return s;
}
void greedyRoute(data d, std::vector<int> q, int &v, boost::dynamic_bitset<> &N, vrs &s){
    //std::cout << "\n\t\t\t\tgreedy";
    int loopCnt = 0;
    int maxLoop = 100;
    while (N.count() > 0 && v >= 1) {
        int oldV = v;
        
        boost::dynamic_bitset<> R (d.n,0);
        int load = d.Q[v];
        int ttt = 0;
        int h = 0;
        std::deque<int> rOrder;

        for(int i = 1; i <= d.n; i++){
            if(N.test(i-1)){
                if (q[i] <= load && ttt - d.a[h][0] + d.a[h][i] + d.s[i] + d.a[i][0] <= d.H){
                    R.set(i-1);
                    N.reset(i-1);
                    rOrder.push_back(i);
                    load = load - q[i];
                    ttt = ttt - d.a[h][0] + d.a[h][i] + d.s[i] + d.a[i][0];
                    h = i;
                    s.cpp[i] = v;
                }
            }
        }

        if (R.count() > 0){
            //salvando nova rota
            //if (ttt != checkRoutTTT(d,rOrder)) std::cout << RED << "\n\tttt1 = " << ttt << " >> TTT = " << checkRoutTTT(d,rOrder);
            s.R[v] = newRoute(d.e[v] + checkRoutCost(d,rOrder), ttt, d.Q[v] - load, rOrder, R);
            s.c += s.R[v].rCost;
            if (s.R[v].rOrder.size() > 0) updatePredSucc(d.n,s.R[v].rOrder,s.ps);
            /*
            std::cout << red;
            printRoute(d,q,s.R[v],v);
            std::cout << normal;
            //*/
            R.reset();
            if (N.any()) v--;
            load = d.Q[v];
            ttt = 0;
            h = 0;
        }
        
        bool term = false;        
        if (oldV == v){
            loopCnt++;
            //std::cout << "\n\t N = " << N << " >> v = " << v << " >> Q = " << d.Q[v];
            //for(int i = 1; i <= d.n; i++) if (N.test(i-1)) std::cout << "\nq[" << i << "]=" << q[i];
            //std::cout << "\n\t loop = " << loopCnt;            
            if (loopCnt > maxLoop && v > 0){
                v--;
                loopCnt = 0;
            }else
                term = true;
        }
        
        if (term) break;
    }
    //std::cout << "\n\tEND GREEDY";
}
//lexicographic greedy route
vrs lgr(data d, routing &r, int T, bool &failed){
    //std::cout << "\n\t\t\tlgr t = " << T << " >> ";
    vrs s = newVRS(d.n,d.V);
    
    //conjunto de clientes ainda não alocados
    boost::dynamic_bitset<> N (d.n, 0);
    for (int i = 1; i <= d.n; i++) if (r.q[T][i] > 0) N.set(i-1);
    if (N.none()) return s;
    
    int v = d.V;
    greedyRoute(d,r.q[T],v,N,s);
    //std::cout << " >> |N| = " << N.count();
    
    if (N.any()){
        failed = true;
        s.anyRoute = false;        
    }else{
        failed = false;
        //atualizando routing r
        r.V[T] = 0;
        for (int v = 1; v <= d.V; v++){
            //std::cout << BLUE << "\n\t\tv = " << v << " >> " << s.R[v].rComp << normal;
            if (s.R[v].rComp.any()){
                r.g[T][v] = true;
                r.z[T][v] = boost::dynamic_bitset<> (s.R[v].rComp);
                r.V[T]++;
                //std::cout << yellow << "\n\tg = " << r.g[T][v];                
            }else{
                r.g[T][v] = false;
                r.z[T][v] = boost::dynamic_bitset<> (d.n, 0);                
            }            
        }
        s.anyRoute = true;        
    }
    
    //if (!checkVRSCost(d,r,s,T,"GREEDY")) failed = true;

    return s;
}

/*CLARKE-WRIGHT SAVING HEURISTIC methods*/
//cria lista de economias para cada periodo
void createCWSavingList(data d, routing r, std::vector<cws> &savings, int T) {
    for (int i = 1; i <= d.n; i++){
        for (int j = 1; j <= d.n; j++){
            if(i != j && i < j){
                if (r.q[T][i] + r.q[T][j] <= d.Qmax && r.q[T][i] > 0 && r.q[T][j] > 0){
                    cws s;
                    s.stop.first = i;
                    s.stop.second = j;
                    s.save = d.c[0][i] + d.c[j][0] - d.c[i][j];
                    savings.push_back(s); //storing the pairs of routess
                }
            }
        }
    }

    //ordena crescentemente os ganhos
    if (savings.size() > 0) sortCW(savings);
}
//ordena em ordem crescente de ganho
void sortCW(std::vector<cws> &savings) {
    for (int fixo = 0; fixo < savings.size() - 1; fixo++) {
        int menor = fixo;

        for (int i = menor + 1; i < savings.size(); i++)
            if (savings[i].save < savings[menor].save) menor = i;

        if (menor != fixo) {
            cws swp = savings[fixo];
            savings[fixo] = savings[menor];
            savings[menor] = swp;
        }
    }
}
void clearSavListCW (std::vector<cws> &savings, int i){
    for (int s = savings.size()-1; s >= 0; s--){
        cws sav = savings[s];
        if (sav.stop.first == i || sav.stop.second == i) savings.erase(savings.begin() + s);
    }
}
//Parallel Clarke-Wright procedure
vrs pcw(data d, routing &r, int T, bool &failed){
    //std::cout << "\n\t\t\tPCW t = " << T << " >> ";
    
    vrs s = newVRS(d.n, d.V);
    
    //indicates which customers already are at some route
    boost::dynamic_bitset<> N (d.n,0);        
    for (int i = 1; i <= d.n; i++) 
        if (r.q[T][i] == 0) N.set(i-1);
        
    if (N.all()) return s;
    
    std::vector<cws> savings;
    createCWSavingList(d,r,savings,T);
    
    if (savings.size() == 0) {
        failed = true;
        return s;
    }
    
    std::vector<cws> savCopy (savings); //corrigindo containers
    //std::cout << "\n\tSaving list = " << savings.size();
    //getchar();
    std::vector<bool> rCreated (d.V + 1, false);
    
        
    //std::cout << "\n\tN = " << N;
    
    int counter = 0;
    do{
        //std::cout << "\n.";
        if (savCopy.size() > 0){
            cws best = savCopy[savCopy.size()-1];
            savCopy.resize(savCopy.size()-1);
            bool added = false;
            int i = best.stop.first;
            int j = best.stop.second;
            //std::cout << "\t(" << i << "," << j << ")";
            for (int v = d.V; v >= 1 && !added; v--){
                if (!rCreated[v] && !N.test(i-1) && !N.test(j-1)){
                    if(r.q[T][i] + r.q[T][j] <= d.Q[v] && d.a[0][i] + d.s[i] + d.a[i][j] + d.s[j] + d.a[j][0] <= d.H){ 
                        //eliminating all unecessary savings
                        clearSavListCW(savCopy,i);
                        clearSavListCW(savCopy,j);
                        
                        boost::dynamic_bitset<> b (d.n,0);
                        b.set(i-1);
                        b.set(j-1);
                        std::deque<int> rt = {i,j};
                        //creating a new route
                        s.R[v] = newRoute(d.e[v] + d.c[0][i] + d.c[i][j] + d.c[j][0], d.a[0][i] + d.s[i] + d.a[i][j] + d.s[j] + d.a[j][0], r.q[T][i] + r.q[T][j], rt, b);
                        /*   
                        std::cout << blue;
                        printRoute(d,s.R[v],v);
                        std::cout << normal;       
                        //*/
                        rCreated[v] = true;
                        
                        r.z[T][v] = boost::dynamic_bitset<> (b);
                        r.g[T][v] = true;
                        
                        N.set(i-1);
                        N.set(j-1);
                        added = true;
                        counter++;
                    }
                }
            }
        }else{
                /*    
                std::cout << RED << "\n\tThe copy of the saving list is empty! t = " << T << "\n\t\tvehicles activated: ";
                for (int v=1; v <= d.V; v++)
                    std::cout << yellow <<"\n\t\t\tv = " << v << " - necessary? " << r.g[T][v] << " - created? "  << rCreated[v] << " -> assignment = " << r.z[T][v] << " | current = " << s.R[v].rComp;

                std::cout << "\n\tEstá sem rota: ";
                for (int i = 1; i <= d.n; i++)
                    if (!N.test(i-1)) std::cout << i << "(" << r.q[T][i] << ") ";
                std::cout << "\n";
                
                for (int v = 1; v <= d.V; v++)
                    if (!rCreated[v]) std::cout << "v = " << v << " >> Q = " << d.Q[v];
                //*/
                
            if (d.n - N.count() == 1 && counter < d.V){
                for (int i = 1; i <= d.n; i++){
                    if (!N.test(i-1)){
                        bool added = false;
                        for (int v = 1; v <= d.V && !added; v++){
                            if (!rCreated[v] && r.q[T][v] <= d.Q[v]) {
                                std::deque<int> rt = {i};
                                s.R[v] = newRoute(d.e[v] + d.c[0][i] + d.c[i][0], d.a[0][i] + d.s[i] + d.a[i][0], r.q[T][i], rt, boost::dynamic_bitset<> (d.n,pow(2,i)));
                                    
                                rCreated[v] = true;
                                r.z[T][v] = boost::dynamic_bitset<> (s.R[v].rComp);
                                r.g[T][v] = true;
                                N.set(i-1);
                                added = true;
                                counter++;
                                /*
                                std::cout << red;
                                printRoute(d,s.R[v],v);
                                std::cout << normal;
                                //*/
                            }
                        }
                    }
                }
            }
            
            if (N.any()){
                failed = true;
                return newVRS(d.n,d.V);                
            }
            
            //getchar();
            //savCopy = savings;            
        }
        
        if (N.all()) break;
    }while(counter < d.V);
    
    if (N.all()){
        //std::cout << "\n\tTodos alocados";
        for (int v = 1; v <= d.V; v++){
            if (r.g[T][v]){
                s.c += s.R[v].rCost;
                for (int i = 1; i <= d.n; i++)
                    if (s.R[v].rComp.test(i-1)) s.cpp[i] = v;                    
                updatePredSucc(d.n, s.R[v].rOrder, s.ps);
            }
        }
        
        if (!checkVRSCost(d,r,s,T,"PCW")) {
            failed = true;
            return newVRS(d.n,d.V);            
        }else return s;
    }
    
    //std::cout << "\n\tN = " << N;
    
    //exit(1);
    int index = savings.size()-1;
    counter = 0;
    boost::dynamic_bitset<> oldISR = N;
    int lstSize = savings.size();
    
    while(N.any()){
        //std::cout << "\n.";
        cws nxtCW = savings[index];
        int i = nxtCW.stop.first;
        int j = nxtCW.stop.second;
        //std::cout << "\t(" << i << "," << j << ")";
        bool added = false;
        for (int v = d.V; v >= 1 && !added; v--){
            if (r.g[T][v]){
                int exclui = 0;
                boost::dynamic_bitset<> b (s.R[v].rComp);
                std::deque<int> rota (s.R[v].rOrder);
                if (N[j-1] == 0 && i == s.R[v].rOrder.front() && s.R[v].rLoad + r.q[T][j] <= d.Q[v] && s.R[v].rComp[j-1] == 0 && s.R[v].rTTT - d.a[0][i] + d.s[j] + d.a[0][j] + d.a[j][i] <= d.H){
                    exclui = i;
                    rota.push_front(j);
                    b.set(j-1);
                    N.set(j-1);                    
                }
                else if (N[j-1] == 0 && i == s.R[v].rOrder.back() && s.R[v].rLoad + r.q[T][j] <= d.Q[v] && s.R[v].rComp[j-1] == 0 && s.R[v].rTTT - d.a[i][0] + d.s[j] + d.a[i][j] + d.a[j][0] <= d.H){
                    exclui = i;
                    rota.push_back(j);
                    b.set(j-1);
                    N.set(j-1);                    
                }else if (N[i-1] == 0 && j == s.R[v].rOrder.front() && s.R[v].rLoad + r.q[T][i] <= d.Q[v] && s.R[v].rComp[i-1] == 0 && s.R[v].rTTT - d.a[0][j] + d.s[i] + d.a[0][i] + d.a[i][j] <= d.H){
                    exclui = j;
                    rota.push_front(i);
                    b.set(i-1);
                    N.set(i-1);                    
                }else if (N[i-1] == 0 && j == s.R[v].rOrder.back() && s.R[v].rLoad + r.q[T][i] <= d.Q[v] && s.R[v].rComp[i-1] == 0 && s.R[v].rTTT - d.a[j][0] + d.s[i] + d.a[j][i] + d.a[i][0] <= d.H){
                    exclui = j;
                    rota.push_back(i);
                    b.set(i-1);
                    N.set(i-1);                    
                }
                
                if (exclui > 0){
                    clearSavListCW(savings,exclui);
                    
                    //creating a new route
                    s.R[v] = newRoute(d.e[v] + checkRoutCost(d,rota),checkRoutTTT(d,rota), checkLoad(rota,r.q[T]),rota, b);
                    /*
                    std::cout << red;
                    printRoute(d,s.R[v],v);
                    std::cout << normal;
                    */
                    r.z[T][v] = boost::dynamic_bitset<> (b);
                    
                    oldISR = N;
                    counter = 0;
                    
                    added = true;
                    savings.erase(savings.begin() + index);
                    if (savings.size() > 0) index = savings.size()-1;
                }                
            }else{
                if (r.q[T][i] + r.q[T][j] <= d.Q[v]){
                    clearSavListCW(savCopy,i);
                    clearSavListCW(savCopy,j);
                    
                    boost::dynamic_bitset<> b (d.n,0);
                    b.set(i-1);
                    b.set(j-1);
                    std::deque<int> rt = {i,j};
                    //creating a new route
                    s.R[v] = newRoute(d.e[v] + d.c[0][i] + d.c[i][j] + d.c[j][0], d.a[0][i] + d.s[i] + d.a[i][j] + d.s[j] + d.a[j][0], r.q[T][i] + r.q[T][j], rt, b);
                    /*
                    std::cout << blue;
                    printRoute(d,s.R[v],v);
                    std::cout << normal;
                    */
                    rCreated[v] = true;
                    r.z[T][v] = boost::dynamic_bitset<> (b);
                    r.g[T][v] = true;
                    r.V[T]++;
                    
                    N.set(i-1);
                    N.set(j-1);
                    added = true;
                }
            }                    
        }
        
        //inicio do trecho para interrupcao de looping infinito
        if (savings.size() != lstSize) lstSize = savings.size();
        else if (oldISR == N) counter++;
        if (counter >= lstSize*d.V*d.n) break;
        if (index > 0) index--;
        else if (savings.size() > 0) index = savings.size()-1;
        //else std::cout << RED << "\nLISTA VAZIA! COMO TRATAR?";        
    }
    
    //std::cout << "\n\tN = " << N;
    oldISR.flip();
    if (oldISR.any()){            
        //std::cout << "\n\tEstá sem rota: ";
        for (int i = 1; i <= d.n; i++){
            if (!N.test(i-1)){
                //std::cout << i << "(" << r.q[T][i] << ") ";
                bool added = false;
                for (int v = 1; v <= d.V && !added; v++){
                    if (s.R[v].rLoad + r.q[T][i] <= d.Q[v] && s.R[v].rTTT - d.a[s.R[v].rLST][0] + d.a[s.R[v].rLST][i] + d.s[i] + d.a[i][0] <= d.H) {
                        std::deque<int> rt (s.R[v].rOrder);
                        rt.push_back(i);
                        boost::dynamic_bitset<> b (s.R[v].rComp);
                        b.set(i-1);
                        s.R[v] = newRoute(d.e[v] + checkRoutCost(d,rt), checkRoutTTT(d,(d,rt)), checkLoad(rt,r.q[T]), rt, b);
                        
                        r.z[T][v] = boost::dynamic_bitset<> (s.R[v].rComp);
                        r.g[T][v] = true;
                        N.set(i-1);
                        added = true;                      
                    }                    
                }                
            }            
        }
        
        if (!N.all()){
            failed = true;
            return newVRS(d.n,d.V);            
        }        
    }
    
    for (int v = 1; v <= d.V; v++){
        s.c += s.R[v].rCost;
        for (int i = 1; i <= d.n; i++)
            if (r.g[T][v])
                if (s.R[v].rComp.test(i-1)) s.cpp[i] = v;
                
        updatePredSucc(d.n, s.R[v].rOrder, s.ps);      
    }
    
    s.anyRoute = true;
        
    //std::cout << "\n\tN = " << N;    
    //if (!checkVRSCost(d,r,s,T,"PCW")) failed = true;
    
    return s;
}
//Sequential Clarke-Wright procedure
vrs scw(data d, routing &r, int T, bool &failed){
    //std::cout << "\n\t\t\tSCW t = " << T << " >> ";
    vrs s = newVRS(d.n, d.V); //nova solucao
    
    boost::dynamic_bitset<> N (d.n, 0);
    for (int i = 1; i <= d.n; i++) if (r.q[T][i] > 0) N.set(i-1);
    if (N.none()) return s;
    
    
    std::vector<cws> savings;
    createCWSavingList(d,r,savings,T);
    std::vector<cws> savCopy (savings); //corrigindo containers
    
    if (savings.size() == 0) {
        failed = true;
        return s;
    }
    
    //std::cout << blue << "\n\tN = " << N << normal;
    //for (int i = 1; i <= d.n; i++) std::cout << green << "\n\ti = " << i << " >> q = " << r.q[T][i] << normal;
    
    bool fail1 = false; //lista de savings vazia e nao criou todas as rotas
    bool fail2 = false; //numero de veiculos insuficientes
    int v = d.V;
    while(N.any()){
        cws best_cw;
        best_cw.save = -1;
        int e = 0;
        for(int c = 0; c < savings.size(); c++){
            cws cw = savings[c];
            int i = cw.stop.first;
            int j = cw.stop.second;
            if(cw.save > best_cw.save  && r.q[T][i] + r.q[T][j] <= d.Q[v] && d.a[0][i] + d.s[i] + d.a[i][j] + d.s[j] + d.a[j][0] <= d.H){
                best_cw = cw;
                e = c;                
            }            
        }
        
        if (best_cw.save == -1) {
            failed = true;
            return newVRS(d.n, d.V);
        }
        
        savings.erase(savings.begin() + e);
        
        int start = best_cw.stop.first;
        int end = best_cw.stop.second;
        int load = r.q[T][start] + r.q[T][end];
        int ttt = d.a[0][start] + d.s[start] + d.a[start][end] + d.s[end] + d.a[end][0];
        std::deque<int> rota = {start,end};
        boost::dynamic_bitset<> n (d.n,0);
        n.set(start-1);
        n.set(end-1);
        N.reset(start-1);
        N.reset(end-1);
        bool full = false;
        if (load == d.Q[v]) full = true;
        
        while(!full){
            cws nxtcw;
            int choose = 0;
            double sav = 0;
            
            for(int c = 0; c < savings.size(); c++){
                cws cw = savings[c];
                int i = cw.stop.first;
                int j = cw.stop.second;
                if (cw.save > sav){
                    if (i == start && load + r.q[T][j] <= d.Q[v] && n[j-1] == 0 && ttt - d.a[0][i] + d.s[j] + d.a[0][j] + d.a[j][i] <= d.H){
                        sav = cw.save;
                        nxtcw = cw;                        
                        choose = 1;                        
                    }
                    else if (i == end && load + r.q[T][j] <= d.Q[v] && n[j-1] == 0 && ttt - d.a[i][0] + d.s[j] + d.a[i][j] + d.a[j][0] <= d.H){
                        sav = cw.save;
                        nxtcw = cw;                        
                        choose = 2;                        
                    }
                    else if (j == start && load + r.q[T][i] <= d.Q[v] && n[i-1] == 0 && ttt - d.a[0][j] + d.s[i] + d.a[0][i] + d.a[i][j] <= d.H){
                        sav = cw.save;
                        nxtcw = cw;                        
                        choose = 3;                        
                    }
                    else if (j == end && load + r.q[T][i] <= d.Q[v] && n[i-1] == 0 && ttt - d.a[j][0] + d.s[i] + d.a[j][i] + d.a[i][0] <= d.H){
                        sav = cw.save;
                        nxtcw = cw;                        
                        choose = 4;                        
                    }                    
                }                
            }
            
            int exclude_i = 0; //eliminar todos os savings com ele            
            int i = nxtcw.stop.first;
            int j = nxtcw.stop.second;
            
            if (choose == 1){
                load += r.q[T][j];
                ttt = ttt - d.a[0][i] + d.s[j] + d.a[0][j] + d.a[j][i];
                if (load == d.Q[v] || ttt == d.H) full = true;
                exclude_i = start;
                start = j;
                N.reset(start-1);
                n.set(j-1);
                rota.push_front(j);                
            }
            else if (choose == 2){
                load += r.q[T][j];
                ttt = ttt - d.a[i][0] + d.s[j] + d.a[i][j] + d.a[j][0];
                if (load == d.Q[v] || ttt == d.H) full = true;
                exclude_i = end;
                end = j;
                N.reset(end-1);
                n.set(j-1);
                rota.push_back(j);                
            }
            else if (choose == 3){
                load += r.q[T][i];
                ttt = ttt - d.a[0][j] + d.s[i] + d.a[0][i] + d.a[i][j];
                if (load == d.Q[v] || ttt == d.H) full = true;
                exclude_i = start;
                start = i;
                N.reset(start-1);
                n.set(i-1);
                rota.push_front(i);                
            }
            else if (choose == 4){
                load += r.q[T][i];
                ttt = ttt - d.a[j][0] + d.s[i] + d.a[j][i] + d.a[i][0];
                if (load == d.Q[v] || ttt == d.H) full = true;
                exclude_i = end;
                end = i;
                N.reset(end-1);
                n.set(i-1);
                rota.push_back(i);                
            }
            else if (choose == 0) full = true;
            if (exclude_i > 0) clearSavListCW(savings,exclude_i);
        }
        
        clearSavListCW(savings,start);
        clearSavListCW(savings,end);
        
        s.R[v] = newRoute(d.e[v] + checkRoutCost(d,rota), ttt, load, rota, n);
        
        s.c += s.R[v].rCost;
        if (s.R[v].rOrder.size() > 0) updatePredSucc(d.n, s.R[v].rOrder,s.ps);
        for (int p = 0; p < s.R[v].rOrder.size(); p++) s.cpp[s.R[v].rOrder[p]] = v;
        if (N.any()) v--;
        
        //interrompendo porque a lista de saving esvaziou
        if (savings.size() == 0){
            if (v >= 1) fail1 = true;
            break;            
        }
        
        //interrompendo porque os clientes sem rota nao couberam em nenhuma
        if(N.any() && v < 1) {
            fail2 = true;
            break;            
        }        
    }
    
    //tratamento para clientes sem rota
    if (fail1) {
        //std::cout << "\nFalha 1: Savings vazio";
        greedyRoute(d,r.q[T],v,N,s);
        //std::cout << "\nSem rota após greedy, |N| = " << N.count() << " " << N;
        if (N.any()) failed = true;        
    }
    
    if (fail2) failed = true;
    
    //getchar();
    if (N.any())
        failed = true;
    else{
        failed = false;
        //atualizando routing r
        r.V[T] = 0;
        for (int v = 1; v <= d.V; v++){
            //std::cout << BLUE << "\n\t\tv = " << v << " >> " << s.R[v].rComp << normal;
            if (s.R[v].rComp.any()){
                r.g[T][v] = true;
                r.z[T][v] = boost::dynamic_bitset<> (s.R[v].rComp);
                r.V[T]++;                
            }else{
                r.g[T][v] = false;
                r.z[T][v] = boost::dynamic_bitset<> (d.n, 0);                
            }            
        }
        s.anyRoute = true;        
    }
    
    //printVRS(d,r,s,T);        
    //if (!checkVRSCost(d,r,s,T,"SCW")) failed = true;
    
    return s;
}
vrs S0(data d, int t, routing &r, bool &failed){
    //std::cout << "\n\t\t\tinitial VRP solution procedure[t=" << t << "]";
    steady_clock::time_point tInitial = steady_clock::now(); //retorna o ponto de agora no tempo
    steady_clock::time_point tFinal = steady_clock::now(); //retorna o ponto de agora no tempo
    duration<double> diffTime = duration_cast<duration<double> > (tFinal-tInitial);

    //(1) PVCW, (2) SCW, (3) LEX GREEDY, (4) MAX GREEDY
    std::vector<int> success (5,1);
    success[0] = 4;

    //enquanto ainda existe algum método de solução inicial
    std::uniform_int_distribution<> IS(1,4);

    while(success[0] > 0){
        vrs s0 = newVRS(d.n,d.V);
        routing rl = cpRout(r);
        int is = IS(gen);
        while(success[is] == 0) is = IS(gen);
        bool melhora = false;
        //std::cout << PURPLE << "\n\tis = " << is;

        if (is == 1) //(1) PVCW
            s0 = pcw(d,rl,t,failed);
        else if (is == 2) //(2) SCW
            s0 = scw(d,rl,t,failed);
        else if (is == 3) //(3) LEX GREEDY
            s0 = lgr(d,rl,t,failed);
        else if (is == 4) //(4) LEX GREEDY
            s0 = mgr(d,rl,t,failed);

        if (!failed){
            tFinal = steady_clock::now(); //retorna o ponto de agora no tempo
            diffTime = duration_cast<duration<double> > (tFinal-tInitial);
            //std::cout << " >> f(s0) = " << s0.c << " >> t(s) = " << diffTime.count() << normal;
            /*
            std::cout << failed;
            std::cout << YELLOW;
            printVRS(d,rl,s0,t);
            std::cout << normal;
            //*/
            r = cpRout(rl);
            sini[0]++;
            sini[is]++;
            return s0;
            break;
        }else{
            success[is] = 0;
            success[0]--;

            if (success[0] == 0){
                failed = true;
                //std::cout << failed;
                return newVRS(d.n,d.V);
                break;
            }else failed = false;
        }
    }    
}

//INTER-ROUTE HEURISTICS
bool soloCustAlloc(data d, vrs &vrsol, routing &r, int t, int v1){
//corrige rotas únicas no VLNS
    //std::cout << "\n\t\t\tSolo Customer Optmization >> ";
    int i = vrsol.R[v1].rOrder[0];
    bool melhorou = false;
    for (int v2 = 1; v2 <= d.V && !melhorou; v2++){
        if (v2 == v1) continue;
        int nTTT = d.a[0][i] + d.s[i] + d.a[0][i];
        if (!r.g[t][v2] && d.e[v2] < d.e[v1] && r.q[t][i] <= d.Q[v2] && nTTT <= d.H){
            melhorou = true;
            double sc = vrsol.c - vrsol.R[v1].rCost;
            //criando nova rota
            route r2 = newRoute(d.n);
            r2.rCost = d.e[v2] + d.c[0][i] + d.c[0][i];
            r2.rLoad = r.q[t][i];
            r2.rTTT = nTTT;
            r2.rOrder.push_back(i);
            r2.rFST = i;
            r2.rLST = i;
            r2.rComp.set(i-1);
            vrsol.R[v2] = r2;

            //atualizando rota desativada
            vrsol.R[v1] = newRoute(d.n);

            //atualizando a solucao
            vrsol.c = sc + vrsol.R[v2].rCost;
            vrsol.cpp[i] = v2;
            updatePredSucc(d.n,vrsol.R[v2].rOrder,vrsol.ps);

            r.g[t][v1] = false;
            r.g[t][v2] = true;
            r.z[t][v1].reset(i-1);
            r.z[t][v2].set(i-1);
        }
    }
    //std::cout << melhorou << " ";
    return melhorou;
}
//vlns
arc newArc(int i, int j, int w, int t, int n){
    arc a;
    a.i = i;
    a.j = j;
    a.w = w;
    a.t = t;
    a.n = n;
    return a;
}
pos newNodePos(int pred, int succ, int indx){
    pos np;
    np.pred = pred;
    np.succ = succ;
    np.indx = indx;
    return np;
}
net improvementGraph(data d, routing r, vrs vrp, int T){
    //std::cout << "\n\t\t>Improvement Graph";

    //retrieving the customer visitation
    boost::dynamic_bitset <> z (d.n, 0);
    for (int v = 1; v <= d.V; v++)
        for (int i = 1; i <= d.n; i++)
            if (r.z[T][v].test(i-1))
                z.set(i-1);

    //graph initialization
    net graph;
    graph.nNods = d.n + d.V + 1;
    graph.nArcs = 0;
    graph.maxCost = -999999;
    graph.src = 0;
    graph.snk = 0;
    graph.ftd = std::vector<std::vector<int> > (graph.nNods);

    //retrieving solution information
    std::vector<route> R (vrp.R);
    std::vector<int> cpp (vrp.cpp);

    //(1) SUBSTITUICAO
    for (int i = 1; i <= d.n; i++){
        if (z.test(i-1)){
            for (int v = 1; v <= d.V; v++){
                if (cpp[i] != v && r.g[T][v]){ //rota deve ser diferente da da que ele se encontra
                    //retrieving route v information
                    std::deque<int> custOrder (R[v].rOrder);
                    double currCost = R[v].rCost;
                    double negCost = (-1)*currCost;
                    int currTT = R[v].rTTT;
                    int rFST = R[v].rFST;
                    int rLST = R[v].rLST;
                    int load = R[v].rLoad;

                    //so existe substituicao se a rota contem algum cliente
                    if (custOrder.size() > 0){
                        for (int c = 0; c < custOrder.size(); c++){
                            int j = custOrder[c];
                            if (load - r.q[T][j] + r.q[T][i] <= d.Q[v]){//viabilidade de carga
                                if (custOrder.size() == 1){
                                    if (currTT - d.a[0][j] - d.s[j] - d.a[j][0] + d.a[0][i] + d.s[i] + d.a[i][0] <= d.H){
                                        int w = currCost - d.c[0][j] - d.c[j][0] + d.c[0][i] + d.c[i][0] + negCost;
                                        graph.arcs.push_back(newArc(i,j,w,1,graph.nArcs));
                                        graph.ftd[i].push_back(graph.nArcs);
                                        graph.nArcs++;
                                        graph.nPos.push_back(newNodePos(0,0,0));

                                        if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo
                                    }
                                }else if (custOrder.size() > 1){
                                    if (c == 0){
                                        int k = custOrder[c+1];
                                        if (currTT - d.a[0][j] - d.s[j] - d.a[j][k] + d.a[0][i] + d.s[i] + d.a[i][k] <= d.H){
                                            int w = currCost - d.c[0][j] - d.c[j][k] + d.c[0][i] + d.c[i][k] + negCost;;
                                            graph.arcs.push_back(newArc(i,j,w,1,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(0,k,c));
                                        }
                                    }else if (c == custOrder.size()-1){
                                        int h = custOrder[c-1];
                                        if (currTT - d.a[h][j] - d.s[j] - d.a[j][0] + d.a[h][i] + d.s[i] + d.a[i][0] <= d.H){
                                            int w = currCost - d.c[h][j] - d.c[j][0] + d.c[h][i] + d.c[i][0] + negCost;
                                            graph.arcs.push_back(newArc(i,j,w,1,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(h,0,c));
                                        }
                                    }else{
                                        int h = custOrder[c-1];
                                        int k = custOrder[c+1];
                                        if (currTT - d.a[h][j] - d.s[j] - d.a[j][k] + d.a[h][i] + d.s[i] + d.a[i][k] <= d.H){
                                            int w = currCost - d.c[h][j] - d.c[j][k] + d.c[h][i] + d.c[i][k] + negCost;
                                            graph.arcs.push_back(newArc(i,j,w,1,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(h,k,c));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //(2) ALOCACAO
    for (int v = 1; v <= d.V; v++){
        if (r.g[T][v]){
            if (R[v].rLoad < d.Q[v]){
                //recupero a informacao do veiculo
                //retrieving route v information
                std::deque<int> custOrder(R[v].rOrder);
                double currCost = R[v].rCost;
                double negCost = (-1)*currCost;
                int currTT = R[v].rTTT;
                int rFST = R[v].rFST;
                int rLST = R[v].rLST;
                int load = R[v].rLoad;
                for (int i = 1; i <= d.n; i++){
                    if (z.test(i-1) && cpp[i] != v){                         //rota deve ser diferente da que ele se encontra
                        if (load + r.q[T][i] <= d.Q[v]){      //existe viabilidade de carga
                            if (custOrder.size() > 1){
                                for (int c = 0; c < custOrder.size(); c++){
                                    double nCost = negCost;     //old cost
                                    int j = custOrder[c];               //customer to be possibly removed from route v
                                    int h = 0;
                                    if (c > 0) h = custOrder[c-1];

                                    if (c == 0){  //making i the first customer at route v
                                        if (currTT - d.a[0][j] + d.a[0][i] + d.a[i][j] + d.s[i] <= d.H){
                                            nCost += (currCost - d.c[0][j] + d.c[0][i] + d.c[i][j]); //adding the new cost
                                            graph.arcs.push_back(newArc(i,d.n + v,nCost,2,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (nCost > graph.maxCost) graph.maxCost = nCost; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(h,j,c));
                                        }
                                    }
                                    else if (c == custOrder.size()-1){
                                        double nCost1 = nCost;
                                        double nCost2 = nCost;
                                        if (currTT - d.a[h][j] + d.a[h][i] + d.a[i][j] + d.s[i] <= d.H){
                                            nCost1 += (currCost - d.c[h][j] + d.c[h][i] + d.c[i][j]); //adding the new cost

                                            graph.arcs.push_back(newArc(i,d.n + v,nCost1,2,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (nCost1 > graph.maxCost) graph.maxCost = nCost1; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(h,j,c));
                                        }

                                        if (currTT - d.a[j][0] + d.a[j][i] + d.a[i][0] + d.s[i] <= d.H){
                                            nCost2 += (currCost - d.c[j][0] + d.c[j][i] + d.c[i][0]); //adding the new cost
                                            graph.arcs.push_back(newArc(i,d.n + v,nCost2,2,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (nCost2 > graph.maxCost) graph.maxCost = nCost2; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(j,0,c+1));
                                        }
                                    }
                                    else{
                                        if (currTT - d.a[h][j] + d.a[h][i] + d.a[i][j] + d.s[i] <= d.H){
                                            nCost += (currCost - d.c[h][j] + d.c[h][i] + d.c[i][j]); //adding the new cost
                                            graph.arcs.push_back(newArc(i,d.n + v,nCost,2,graph.nArcs));
                                            graph.ftd[i].push_back(graph.nArcs);
                                            graph.nArcs++;

                                            if (nCost > graph.maxCost) graph.maxCost = nCost; //atualizando o arco de maior custo

                                            graph.nPos.push_back(newNodePos(h,j,c));
                                        }
                                    }
                                }
                            }
                            else if (custOrder.size() == 1){     //rota de um único cliente
                                int j = rFST;
                                //i eh adicionado ao inicio
                                if (currTT - d.a[0][j] + d.a[0][i] + d.a[i][j] + d.s[i] <= d.H){
                                    int w = currCost - d.c[0][j] + d.c[0][i] + d.c[i][j] + negCost;
                                    graph.arcs.push_back(newArc(i,d.n + v,w,2,graph.nArcs));
                                    graph.ftd[i].push_back(graph.nArcs);
                                    graph.nArcs++;

                                    if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                                    graph.nPos.push_back(newNodePos(0,j,0));
                                }

                                //i eh adicionado ao fim
                                if (currTT - d.a[j][0] + d.a[j][i] + d.a[i][0] + d.s[i] <= d.H){
                                    int w = currCost - d.c[0][j] + d.c[0][i] + d.c[i][j] + negCost;
                                    graph.arcs.push_back(newArc(i,d.n + v,w,2,graph.nArcs));
                                    graph.ftd[i].push_back(graph.nArcs);
                                    graph.nArcs++;

                                    if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                                    graph.nPos.push_back(newNodePos(j,0,1));
                                }
                            }
                        }
                    }
                }
            }
        }
        else{
            //o veiculo atualmente nao eh utilizado
            for (int i = 1; i <= d.n; i++){
                if (z.test(i-1)){
                    if (r.q[T][i] <= d.Q[v]){
                        int w = d.e[v] + d.c[0][i] + d.c[i][0];
                        graph.arcs.push_back(newArc(i,d.n + v,w,2,graph.nArcs));
                        graph.ftd[i].push_back(graph.nArcs);
                        graph.nArcs++;

                        if (w > graph.maxCost) graph.maxCost = w; //atualizando o arco de maior custo

                        graph.nPos.push_back(newNodePos(0,0,0));
                    }
                }
            }
        }
    }

    //getchar();

    //(3) REMOCAO
    //std::cout << "\n\t\t\t>>(3) arcos de remocao";
    for (int v = 1; v <= d.V; v++){
        if (r.g[T][v]){
            std::deque<int> custOrder(R[v].rOrder); //retrieving the visitation order of customers
            double currCost = R[v].rCost;            // and the current cost
            double negCost = (-1)*currCost;     //old cost to be subtracted
            if (custOrder.size() > 1){
                for (int c = 0; c < custOrder.size(); c++){
                    int j = custOrder[c];               //customer to be possibly removed from route v
                    double nCost = negCost;
                    int h = 0;
                    if (c > 0) h = custOrder[c-1];
                    int k = 0;
                    if (c < custOrder.size()-1) k = custOrder[c+1];
                    if (c == 0)//removing the first customer from the route v
                        nCost += (currCost - d.c[h][j] - d.c[j][k] + d.c[h][k]);
                    else if (c == custOrder.size()-1)
                        nCost += (currCost - d.c[h][j] - d.c[j][0] + d.c[h][0]);
                    else
                        nCost += (currCost - d.c[h][j] - d.c[j][k] + d.c[h][k]);

                    if (nCost > graph.maxCost) graph.maxCost = nCost; //atualizando o arco de maior custo

                    graph.arcs.push_back(newArc(0,j,nCost,3,graph.nArcs));
                    graph.ftd[0].push_back(graph.nArcs);
                    graph.nArcs++;
                    graph.nPos.push_back(newNodePos(h,k,c));
                }
            }else if (custOrder.size() == 1){
                int j = custOrder[0];               //customer to be possibly removed from route v

                if (negCost > graph.maxCost) graph.maxCost = negCost; //atualizando o arco de maior custo

                graph.arcs.push_back(newArc(0,j,negCost,3,graph.nArcs));
                graph.ftd[0].push_back(graph.nArcs);
                graph.nArcs++;
                graph.nPos.push_back(newNodePos(0,0,0));
            }
        }
    }


    //(4) Complementar de adicao
    for (int v = 1; v <= d.V; v++){
        if (0 > graph.maxCost) graph.maxCost = 0; //atualizando o arco de maior custo

        graph.arcs.push_back(newArc(d.n + v,0,0,4,graph.nArcs));
        graph.ftd[d.n + v].push_back(graph.nArcs);
        graph.nArcs++;
        graph.nPos.push_back(newNodePos(0,0,0));
    }
    return graph;
}
std::vector<int> labelCorrecting(data d, net g){
    //std::cout << "\n\t\t>Label Correcting";

    //inicializacao
    std::vector<int> pred (g.nNods, 0);
    std::vector<double> dist (g.nNods, 9999999);
    std::vector<int> wasInTheList (g.nNods, 0);
    std::vector<int> minArcList (g.nNods, 0);
    boost::dynamic_bitset<> inTheList (g.nNods, 0);


    pred[0] = 0;
    dist[0] = 0.0;

    std::deque<int> LIST;
    LIST.push_back(0);
    inTheList.set(0);
    //wasInTheList[0]++;

    bool ctn = true;
    int i,i0;

    //exploring the nodes and their successors
    //std::cout << "\n\t\t\t>exploring the nodes and their successors..." << std::endl;
    while (LIST.size() > 0 && ctn){
        i = LIST[0];
        wasInTheList[i]++;

        if (wasInTheList[i] > g.nNods - 1){
            i0 = pred[i];
            ctn = false;
            break;
        }
        //acompanhando as alteracoes de gustavo
        LIST.pop_front();
        inTheList.reset(i);

        for (int ca = 0; ca < g.ftd[i].size(); ca++){
            arc a = g.arcs[g.ftd[i][ca]];
            if(dist[a.j] > dist[i] + a.w){
                dist[a.j] = dist[i] + a.w;      //custo acumulado
                pred[a.j] = i;                  //predecessor
                minArcList[a.j] = g.ftd[i][ca]; //arco utilizado, pois existem arcos paralelos
                if (!inTheList.test(a.j)){
                    inTheList.set(a.j);
                    LIST.push_back(a.j);
                }
            }
        }
    }

    std::vector<int> cycle;
    boost::dynamic_bitset<> cBit (g.nNods, 0);
    std::vector<int> inTheCycle (g.nNods, 0);
    if (!ctn){
        recoverCycle(i, pred, minArcList, cycle, cBit);
        for (int ia = 0; ia < cycle.size(); ia++){
            arc a = g.arcs[cycle[ia]];
            inTheCycle[a.i]++;
            inTheCycle[a.j]++;
            if (inTheCycle[a.i] > 2){
                cycle.clear();
                cBit.reset();
                i = a.i;
                recoverCycle(i, pred, minArcList, cycle, cBit);
            }else if (inTheCycle[a.j] > 2){
                i = a.j;
                recoverCycle(i, pred, minArcList, cycle, cBit);
            }
        }
    }
    return cycle;
}
void recoverCycle(int i, std::vector<int> pred, std::vector<int> minArcList, std::vector<int> &cycle, boost::dynamic_bitset<> &cBit){
    //std::cout << "\n\t\t>recovering the cycle..." << std::endl;
    int rFST = i;
    cycle.push_back(minArcList[i]);
    cBit.set(i);
    i = pred[i];
    while (i != rFST){
        if (!cBit.test(i)){
            cBit.set(i);
            cycle.push_back(minArcList[i]);
            if (pred[i] == i) {
                cycle.clear();
                cBit.reset();
                break;
            }
            i = pred[i];
        }
        else
            break;
    }
}
bool cycleValidation(data d, std::vector<int> cpp, net g, std::vector<int> &cycle){
    //std::cout << "\n\t\t>validating the cycle: ";
    std::vector<int> routes (d.V+1, 0);
    bool valid = true;

    std::vector<int> inTheCycleI (g.nNods, 0);
    std::vector<int> inTheCycleJ (g.nNods, 0);

    if (cycle.size() > 1){
        for(int ia = 0; ia < cycle.size() && valid; ia++){
            arc a = g.arcs[cycle[ia]];

            //std::cout << lgray << "\n\t\t\t>> #arc: " << a.n << " -> i: " << a.i <<  " -> j: " << a.j <<  " -> w: " << a.w << " -> t: " << a.t << normal;
            inTheCycleI[a.i]++;
            inTheCycleJ[a.j]++;
            if (inTheCycleI[a.i] > 1 || inTheCycleJ[a.j] > 1) valid = false;

            //se alguma particao for utilizada mais do que 1 vez, então o ciclo é inválido
            if (a.i > 0 && a.i <= d.n){ //vejo a rota que um cliente está
                int v = cpp[a.i];       //a origem do arco eh um cliente
                routes[v]++;            //sinaliza-se que a particao/rota/veiculo foi usada
                if (routes[v] > 1)      //se foi usada mais de uma vez, o ciclo eh invalido
                    valid = false;
            }else if (a.i > d.n && a.i <= d.n + d.V){ // vejo a rota que recebeu ou perdeu um cliente
                int v = a.i - d.n;      //a origem do arco eh uma particao
                routes[v]++;            //sinaliza-se que a particao/rota/veiculo foi usada
                if (routes[v] > 1)      //se foi usada mais de uma vez, o ciclo eh invalido
                    valid = false;
            }
        }
    }
    else
        valid = false;

    if (!valid) cycle.clear();

    return valid;
}
void buildSolution(data d, std::vector<int> &cycle, net g, routing &rl, vrs &sol, int T){
    //std::cout << "\n\t\t>Building the new solution: ";
    std::vector<route> R (sol.R);
    std::vector<int> cpp (sol.cpp);
    double cost = sol.c;
    std::vector<std::pair<int,int> > ps (sol.ps);
    routing r = cpRout(rl);

    for (int ic = 0; ic < cycle.size(); ic++){
        arc a = g.arcs[cycle[ic]];
        if (a.t == 1){
            //substituicao i->j
            int v = sol.cpp[a.j];
            cpp[a.i] = v;
            cost += a.w;

            /*
            std::cout << purple << "\n" << cost;
            std::cout << blue;
            printRoute(d,R[v],v);
            std::cout << normal;
            */

            r.z[T][v].set(a.i-1);
            r.z[T][v].reset(a.j-1);
            r.z[T][sol.cpp[a.i]].reset(a.i-1);

            int l = g.nPos[a.n].indx;   //localization
            int p = g.nPos[a.n].pred;   //predecessor
            int s = g.nPos[a.n].succ;   //successor
            if (s == 0) s = d.n+1;

            ps[a.i] = std::pair<int,int> (p,s);
            if (p > 0) ps[p].second = a.i;
            if (s < d.n + 1) ps[s].first = a.i;

            //R[v].rCost = R[v].rCost - d.c[p][a.j] - d.c[a.j][s] + d.c[p][a.i] + d.c[a.i][s];
            R[v].rLoad = R[v].rLoad - r.q[T][a.j] + r.q[T][a.i];
            if (l == 0) R[v].rFST = a.i;
            if (l == R[v].rOrder.size()-1) R[v].rLST = a.i;
            R[v].rOrder[l] = a.i;
            R[v].rCost = d.e[v] + checkRoutCost(d,R[v].rOrder);
            R[v].rComp.reset(a.j-1);
            R[v].rComp.set(a.i-1);
            //R[v].rTTT = R[v].rTTT - d.a[p][a.j] - d.a[a.j][s] + d.a[p][a.i] + d.a[a.i][s] - d.s[a.j] + d.s[a.i];
            R[v].rTTT = checkRoutTTT(d,R[v].rOrder);

            //teste
            updatePredSucc(d.n,R[v].rOrder,ps);
            /*
            std::cout << red;
            printRoute(d,R[v],v);
            std::cout << normal;
            */
        }
        else if (a.t == 2){
            //alocacao i->v
            int v = a.j - d.n;
            cpp[a.i] = v;
            cost += a.w;

            /*
            std::cout << purple << "\n" << cost;
            std::cout << blue;
            printRoute(d,R[v],v);
            std::cout << normal;
            */

            r.z[T][v].set(a.i-1);
            r.z[T][sol.cpp[a.i]].reset(a.i-1);

            int l = g.nPos[a.n].indx;   //localization/position
            int p = g.nPos[a.n].pred;   //predecessor
            int s = g.nPos[a.n].succ;   //successor
            if (s == 0) s = d.n+1;

            ps[a.i] = std::pair<int,int> (p,s);
            if (p > 0) ps[p].second = a.i;
            if (s < d.n + 1) ps[s].first = a.i;

            R[v].rCost = R[v].rCost - d.c[p][s] + d.c[p][a.i] + d.c[a.i][s];
            R[v].rLoad = R[v].rLoad + r.q[T][a.i];

            if (l == 0) R[v].rFST = a.i;
            if (l == R[v].rOrder.size()-1) R[v].rLST = a.i;

            R[v].rOrder.emplace(R[v].rOrder.begin()+l,a.i);
            //R[v].rCost = d.e[v] + checkRoutCost(d,R[v].rOrder);
            R[v].rComp.set(a.i-1);
            R[v].rTTT = R[v].rTTT - d.a[p][s] + d.a[p][a.i] + d.a[a.i][s]  + d.s[a.i];
            //R[v].rTTT = checkRoutTTT(d,R[v].rOrder);
            if (!r.g[T][v]){
                r.g[T][v] = true;
                R[v].rCost = d.e[v] + d.c[p][a.i] + d.c[a.i][s];
                R[v].rTTT = d.a[p][a.i] + d.a[a.i][s] + d.s[a.i];
            }
            //teste
            updatePredSucc(d.n,R[v].rOrder,ps);
            /*
            std::cout << red;
            printRoute(d,R[v],v);
            std::cout << normal;
            */
        }
        else if (a.t == 3){
            //remocao
            int v = sol.cpp[a.j];

            /*
            std::cout << blue;
            printRoute(d,R[v],v);
            std::cout << normal;
            */

            cost += a.w;
            //std::cout << purple << "\n" << cost;
            r.z[T][v].reset(a.j-1);

            int l = g.nPos[a.n].indx;   //localization/position
            int p = g.nPos[a.n].pred;   //predecessor
            int s = g.nPos[a.n].succ;   //successor
            if (s == 0) s = d.n+1;

            //na duvida da necessidade desta operacao
            if (p > 0) ps[p].second = s;
            if (s < d.n + 1) ps[s].first = p;

            R[v].rCost = R[v].rCost - d.c[p][a.j] - d.c[a.j][s] + d.c[p][s];
            R[v].rLoad = R[v].rLoad - r.q[T][a.j];
            R[v].rComp.reset(a.j-1);
            R[v].rTTT = R[v].rTTT - d.a[p][a.j] - d.a[a.j][s] + d.a[p][s] - d.s[a.j];

            // erasing the lth element
            R[v].rOrder.erase(R[v].rOrder.begin()+l);
            if (R[v].rOrder.size() == 0){
                r.g[T][v] = false;
                r.z[T][v] = boost::dynamic_bitset<> (R[v].rComp);
                R[v].rCost -= d.e[v];
                //r.z[T][v] = boost::dynamic_bitset<> (d.n,0);
                //R[v] = newRoute(d.n);
            }

            if (l == 0) R[v].rFST = s;
            else if (l == R[v].rOrder.size()-1) R[v].rLST = p;
            //teste
            if (R[v].rOrder.size() > 0) updatePredSucc(d.n,R[v].rOrder,ps);
            /*
            std::cout << red;
            printRoute(d,R[v],v);
            std::cout << normal;
            */
        }
        else if (a.t == 4){
            //complemento de adicao
        }
    }

    bool melhorou = true;

    if (melhorou){
        sol.R = std::vector<route>(R);
        sol.cpp = std::vector<int> (cpp);
        sol.c = cost;
        sol.ps = std::vector<std::pair<int,int> > (ps);
        rl = cpRout(r);
        //if (!checkVRSCost(d,rl,sol,T,"VLNS")) melhorou = false;
    }
}
net clearGraph(){
    net g;
    g.nArcs = 0;
    g.nNods = 0;
    g.maxCost = 999999;
    g.src = 0;
    g.snk = 0;
    g.arcs = std::vector<arc> (0);
    g.nPos = std::vector<pos> (0);
    g.ftd = std::vector<std::vector<int> > (0);
    return g;
}
bool vlns(data d, routing &r, int T, vrs &s){
    //std::cout << "\n\t\t\tVLNS[t=" << T << "] >> ";
    bool melhorou = false;
    bool melhora = true;

    while(melhora){
        vrs sl = cpyVRS(s);
        net g = improvementGraph(d, r, sl, T);

        bool valid = true;
        while (valid){
            std::vector<int> cycle = labelCorrecting(d, g);

            if (cycle.size() > 1)
                valid = cycleValidation(d, sl.cpp, g, cycle);
            else
                valid = false;

            if (valid){
                buildSolution(d, cycle, g, r, sl, T);
                g = clearGraph();
                g = improvementGraph(d, r, sl, T);
                melhorou = true;
            }else
                melhora = false;
        }

        bool solo = false;
        for (int v = 1; v <= d.V; v++){
            if (r.g[T][v] && sl.R[v].rOrder.size() == 1){
                solo = soloCustAlloc(d,sl,r,T,v);
                if(solo) melhora = true;
            }
        }

        s = cpyVRS(sl);
    }
    //std::cout << melhorou;
    return melhorou;
}

//K clientes adjacentes são transferidos para o final de uma outra rota
bool shift(data d, vrs &s, routing &r, int t, int K){
    //std::cout << "\n\t\t\tSHIFT(K=" << K << ",t="<< t << ") >> ";
    //sinaliza melhora
    bool melhora = false;
    bool melhorou = false;

    do {
        //armazena a melhor alteracao
        int bV1, bV2;
        std::deque<int> bR1, bR2;
        double bestDiff = 0;

        melhora = false;
        for (int v1 = 1; v1 <= d.V; v1++){
            if (s.R[v1].rOrder.size() < K) continue;
            for (int pos1 = 0; pos1 < s.R[v1].rOrder.size() - K; pos1++){
                //current vector of transferred customers
                std::vector<int> currK (K,0);
                int demK = 0;
                for (int k = 0; k < K; k++) {
                    currK[k] = s.R[v1].rOrder[pos1 + k];
                    demK += r.q[t][currK[k]];
                }

                for (int v2 = 1; v2 <= d.V; v2++){
                    //teste da demanda
                    int newLoad = demK + s.R[v2].rLoad;
                    if (v1 != v2 && newLoad <= d.Q[v2]){
                        std::deque<int> r2 (s.R[v2].rOrder);
                        for (int k = 0; k < K; k++) r2.push_back(currK[k]);
                        //teste do tempo de viagem, adicionando string ao final
                        if (checkRoutTTT(d,r2) <= d.H){
                            std::deque<int> r1 (s.R[v1].rOrder);

                            for (int k = 0; k < K; k++){
                                int i = currK[k];
                                r1.erase(r1.begin() + pos1);
                            }

                            int cr2 = d.e[v2] + checkRoutCost(d,r2);
                            int cr1 = 0;
                            if (r1.size() > 0) cr1 = d.e[v1] + checkRoutCost(d,r1);

                            double diffCosts = ((cr1 + cr2) - (s.R[v1].rCost + s.R[v2].rCost));
                            if (diffCosts < bestDiff){
                                melhora = true;
                                melhorou = true;
                                bestDiff = diffCosts;

                                bV1 = v1;
                                bV2 = v2;

                                bR1 = std::deque<int> (r1);
                                bR2 = std::deque<int> (r2);
                            }
                        }
                    }
                }
            }
        }

        if (melhora){                                    
            //atualizando solução
            double oldV1Cost = s.R[bV1].rCost;
            double oldV2Cost = s.R[bV2].rCost;

            s.c += bestDiff;

            boost::dynamic_bitset<> btR2(d.n,0);
            for (int pos = 0; pos < bR2.size(); pos++) {
                s.cpp[bR2[pos]] = bV2;
                btR2.set(bR2[pos]-1);
            }
            s.R[bV2] = newRoute(d.e[bV2] + checkRoutCost(d,bR2), checkRoutTTT(d,bR2), checkLoad(bR2,r.q[t]), bR2, btR2);
            r.z[t][bV2] = boost::dynamic_bitset<> (btR2);
            r.g[t][bV2] = true;

            if (bR1.size() > 0){
                boost::dynamic_bitset<> btR1(d.n,0);
                for (int pos = 0; pos < bR1.size(); pos++) {
                    s.cpp[bR1[pos]] = bV1;
                    btR1.set(bR1[pos]-1);
                }
                s.R[bV1] = newRoute(d.e[bV1] + checkRoutCost(d,bR1), checkRoutTTT(d,bR1), checkLoad(bR1,r.q[t]), bR1, btR1);
                r.z[t][bV1] = boost::dynamic_bitset<> (btR1);
                r.g[t][bV1] = true;
            }else{
                s.R[bV1] = newRoute(d.n);
                r.z[t][bV1] = boost::dynamic_bitset<> (d.n,0);
                r.V[t]--;
                r.g[t][bV1] = false;
            }

            if (s.R[bV1].rOrder.size() >= 1) updatePredSucc(d.n,s.R[bV1].rOrder,s.ps);
            if (s.R[bV2].rOrder.size() >= 1) updatePredSucc(d.n,s.R[bV2].rOrder,s.ps);                                    
        }
    }while(melhora);

    //if (melhorou && !checkVRSCost(d,r,s,t,"SHIFT(K)")) melhorou = false;

    //std::cout << melhorou;
    return melhorou;
}
bool cross(data d, vrs &s, routing &r, int t){
    //std::cout << "\n\t\t\tCROSS(t="<< t << ") >> ";
    //armazena a melhor alteracao
    bool melhora = false;
    bool melhorou = false;

    do{
        int bV1, bV2;
        //std::vector<std::pair<int,int> > bestPS (s.ps);
        std::vector<std::pair<int,int> > bestPS;
        route bestV1 = newRoute(d.n);
        route bestV2 = newRoute(d.n);
        routing bestR;
        std::vector<int> bstV1;
        std::vector<int> bstV2;
        double bestDiff = 0;

        melhora = false;
        for (int v1 = 1; v1 <= d.V; v1++){
            if (s.R[v1].rOrder.size() <= 2) continue;

            int ss1 = std::ceil(s.R[v1].rOrder.size()/2);

            std::deque<int> r1;
            boost::dynamic_bitset<> B1(d.n,0);
            int d1 = 0;

            for (int pos = 0; pos < ss1; pos++){
                int i = s.R[v1].rOrder[pos];
                r1.push_back(i);
                B1.set(i-1);
                d1 += r.q[t][i];
            }

            std::vector<int> nV2;
            std::deque<int> r2;
            boost::dynamic_bitset<> B2(d.n,0);
            int d2 = 0;
            for (int pos = ss1; pos < s.R[v1].rOrder.size(); pos++){
                int i = s.R[v1].rOrder[pos];
                r2.push_back(i);
                B2.set(i-1);
                d2 += r.q[t][i];
                nV2.push_back(i);
            }

            for (int v2 = 1; v2 <= d.V; v2++){
                if (s.R[v2].rOrder.size() <= 2 || v1 == v2) continue;

                int ss2 = std::ceil(s.R[v2].rOrder.size()/2);

                std::vector<int> nV1;
                for (int pos = 0; pos < ss2; pos++){
                    int i = s.R[v2].rOrder[pos];
                    r1.push_back(i);
                    B1.set(i-1);
                    d1 += r.q[t][i];
                    if (d1 > d.Q[v1]) break;
                    nV1.push_back(i);
                }

                if (d1 > d.Q[v1]) continue;

                for (int pos = s.R[v2].rOrder.size() - 1; pos >= ss2; pos--){
                    int i = s.R[v2].rOrder[pos];
                    r2.push_front(i);
                    B2.set(i-1);
                    d2 += r.q[t][i];
                    if (d2 > d.Q[v2]) break;
                }

                if (d2 > d.Q[v2]) continue;

                int cr1 = checkRoutCost(d,r1);
                //int ttt1 = cr1 + d.s[1]*r1.size();
                int ttt1 = checkRoutTTT(d,r1);

                if (ttt1 > d.H) continue;
                cr1 += d.e[v1];

                int cr2 = checkRoutCost(d,r2);
                //int ttt2 = cr2 + d.s[1]*r2.size();
                int ttt2 = checkRoutTTT(d,r2);
                if (ttt2 > d.H) continue;
                cr2 += d.e[v2];

                double diffCosts = ((cr1 + cr2) - (s.R[v1].rCost + s.R[v2].rCost));
                //if (ttt1 <= d.H && ttt2 <= d.H && diffCosts < 0){
                if (ttt1 <= d.H && ttt2 <= d.H && diffCosts < bestDiff){
                    melhora = true;
                    melhorou = true;

                    bV1 = v1;
                    bV2 = v2;

                    bestV1 = newRoute(cr1, ttt1, d1, r1, B1);
                    bestV2 = newRoute(cr2, ttt2, d2, r2, B2);

                    bestDiff = diffCosts;

                    bstV1 = std::vector<int>(nV1);
                    bstV2 = std::vector<int>(nV2);

                    bestPS = std::vector<std::pair<int,int> > (s.ps);
                    if (bestV1.rOrder.size() > 0) updatePredSucc(d.n,bestV1.rOrder,bestPS);
                    if (bestV2.rOrder.size() > 0) updatePredSucc(d.n,bestV2.rOrder,bestPS);
                }
            }
        }

        if (melhora){
            //atualizando solução
            double oldV1Cost = s.R[bV1].rCost;
            double oldV2Cost = s.R[bV2].rCost;
            /*
            std::cout << yellow << "\n\told cost = " << s.c;
            std::cout << red;
            printRoute(d,s.R[bV1],bV1);
            printRoute(d,s.R[bV2],bV2);
            std::cout << green;
            printRoute(d,bestV1,bV1);
            printRoute(d,bestV2,bV2);
            std::cout << normal;
            */

            s.R[bV1] = cpyRoute(bestV1);
            s.R[bV2] = cpyRoute(bestV2);
            s.ps = std::vector<std::pair<int,int> > (bestPS);
            s.c += bestDiff;


            for (int k = 0; k < bstV1.size(); k++) s.cpp[bstV1[k]] = bV1;
            for (int k = 0; k < bstV2.size(); k++) s.cpp[bstV2[k]] = bV2;

            r.z[t][bV1] = boost::dynamic_bitset<> (s.R[bV1].rComp);
            r.z[t][bV2] = boost::dynamic_bitset<> (s.R[bV2].rComp);

            //std::cout << yellow << "\n\tnew cost = " << s.c;
        }
    }while(melhora);

    //if (melhorou && !checkVRSCost(d,r,s,t,"CROSS")) melhorou = false;

    //std::cout << melhorou;
    return melhorou;
}
bool shift20(data d, vrs &s, routing &r, int t){
    //std::cout << "\n\t\t\tSHIFT20(t="<< t << ") >> ";
    //sinaliza melhoria
    bool melhora = false;
    bool melhorou = false;

    do{
        //armazena a melhor alteracao
        double bestDiff = 0;
        int bV1, bV2;
        std::deque<int> bR1, bR2;                

        melhora = false;
        for (int v1 = 1; v1 <= d.V; v1++){
            if (s.R[v1].rOrder.size() < 2) continue;
            for (int pos1 = 0; pos1 < s.R[v1].rOrder.size() - 1; pos1++){
                int i = s.R[v1].rOrder[pos1];
                int j = s.R[v1].rOrder[pos1 + 1];
                int loadIJ = r.q[t][i] + r.q[t][j];
                
                std::deque<int> r1 (s.R[v1].rOrder);
                r1.erase(r1.begin() + pos1);
                r1.erase(r1.begin() + pos1);
                if (checkRoutTTT(d,r1) > d.H) break;

                for (int v2 = 1; v2 <= d.V; v2++){
                    if (v1 != v2 && (loadIJ + s.R[v2].rLoad <= d.Q[v2])){
                        for (int pos2 = 0; pos2 < s.R[v2].rOrder.size(); pos2++){
                            //teste da viabilidade de tempo de viagem
                            std::deque<int> r2 (s.R[v2].rOrder);
                            r2.emplace(r2.begin() + pos2,i);
                            r2.emplace(r2.begin() + pos2 + 1,j);
                            
                            if (checkRoutTTT(d,r2) > d.H) break;
                                                                            
                            double diffCosts = 0;
                            
                            if (r1.size() > 0)
                                diffCosts = ((d.e[v1] + checkRoutCost(d,r1) + d.e[v2] + checkRoutCost(d,r2)) - (s.R[v1].rCost + s.R[v2].rCost));
                            else 
                                diffCosts = (d.e[v2] + checkRoutCost(d,r2) - (s.R[v1].rCost + s.R[v2].rCost));
                            
                            if (diffCosts < bestDiff) {
                                //std::cout << "\n\tMelhoria: (" << i << "," << j << ") alocados do veiculo v1 = " << v1 << " para v2 = " << v2 << ". Melhoria = " << diffCosts;

                                melhora = true;
                                melhorou = true;

                                bV1 = v1;
                                bV2 = v2;

                                bestDiff = diffCosts;
                                
                                bR1 = std::deque<int> (r1);
                                bR2 = std::deque<int> (r2);
                            }

                        }
                    }
                }
            }
        }

        if (melhora){
            //atualizando solução
            double oldV1Cost = s.R[bV1].rCost;
            double oldV2Cost = s.R[bV2].rCost;
            
            //std::cout << yellow << "\n\told cost = " << s.c;
            
            s.c += bestDiff;
            
            /*                  
            std::cout << red;
            printRoute(d,s.R[bV1],bV1);
            printRoute(d,s.R[bV2],bV2);
            //*/
            
            //atualizando a rota que recebeu os dois clientes
            boost::dynamic_bitset<> btR2(d.n,0);
            for (int pos = 0; pos < bR2.size(); pos++) {
                s.cpp[bR2[pos]] = bV2;
                btR2.set(bR2[pos]-1);
            }
            
            s.R[bV2] = newRoute(d.e[bV2] + checkRoutCost(d,bR2), checkRoutTTT(d,bR2), checkLoad(bR2,r.q[t]), bR2, btR2);
            r.z[t][bV2] = boost::dynamic_bitset<> (btR2);
            r.g[t][bV2] = true;
            
            //atualizando a rota que perdeu os dois clientes            
            if (bR1.size() > 0){
                boost::dynamic_bitset<> btR1(d.n,0);
                for (int pos = 0; pos < bR1.size(); pos++) {
                    s.cpp[bR1[pos]] = bV1;
                    btR1.set(bR1[pos]-1);
                }
                s.R[bV1] = newRoute(d.e[bV1] + checkRoutCost(d,bR1), checkRoutTTT(d,bR1), checkLoad(bR1,r.q[t]), bR1, btR1);
                r.z[t][bV1] = boost::dynamic_bitset<> (btR1);
                r.g[t][bV1] = true;
            }else{
                s.R[bV1] = newRoute(d.n);
                r.z[t][bV1] = boost::dynamic_bitset<> (d.n,0);
                r.V[t]--;
                r.g[t][bV1] = false;
            }

            if (s.R[bV1].rOrder.size() >= 1) updatePredSucc(d.n,s.R[bV1].rOrder,s.ps);
            if (s.R[bV2].rOrder.size() >= 1) updatePredSucc(d.n,s.R[bV2].rOrder,s.ps);  
            
            /*                        
            std::cout << green;
            printRoute(d,s.R[bV1],bV1);
            printRoute(d,s.R[bV2],bV2);
            std::cout << normal;
            //std::cout << yellow << "\n\tnew cost = " << s.c;
            //*/                                    
        }
    }while(melhora);

    //if (melhorou && !checkVRSCost(d,r,s,t,"SHIFT20")) melhorou = false;

    //std::cout << melhorou;
    return melhorou;
}
//dois clientes adjacentes são trocados por um clientes de outra rota
bool swap21(data d, vrs &s, routing &r, int t){
    //std::cout << "\n\t\t\tSWAP21(t="<< t << ") >> ";
    //sinaliza melhoria
    bool melhora = false;
    bool melhorou = false;

    do{
        //armazena a melhor alteracao
        int bV1, bV2, bestI, bestJ, bestK;
        //std::vector<std::pair<int,int> > bestPS (s.ps);
        std::vector<std::pair<int,int> > bestPS;
        route bestV1, bestV2;
        routing bestR;
        double bestDiff = 0;

        melhora = false;
        for (int v1 = 1; v1 <= d.V; v1++){
            if (s.R[v1].rOrder.size() < 2) continue;

            for (int pos1 = 0; pos1 < s.R[v1].rOrder.size() - 1; pos1++){
                int i = s.R[v1].rOrder[pos1];
                int j = s.R[v1].rOrder[pos1 + 1];
                int loadIJ = r.q[t][i] + r.q[t][j];

                for (int v2 = 1; v2 <= d.V; v2++){
                    if (s.R[v2].rOrder.size() < 2 || v1 == v2) continue;

                    for (int pos2 = 0; pos2 < s.R[v2].rOrder.size(); pos2++){
                        int k = s.R[v2].rOrder[pos2];
                        //se a troca entre (i,j) e k eh viavel, entao tentaremos ver se melhora a solução
                        if ((s.R[v2].rLoad - r.q[t][k] + loadIJ <= d.Q[v2]) && (s.R[v1].rLoad + r.q[t][k] - loadIJ <= d.Q[v1])){
                            //teste da viabilidade de tempo de viagem para R2
                            std::deque<int> r2 (s.R[v2].rOrder);
                            r2.erase(r2.begin() + pos2);
                            r2.emplace(r2.begin() + pos2,i);
                            r2.emplace(r2.begin() + pos2 + 1,j);

                            int cr2 = checkRoutCost(d,r2);
                            //int ttt2 = cr2 + d.s[1]*r2.size();
                            int ttt2 = checkRoutTTT(d,r2);
                            if (ttt2 > d.H) break;

                            //teste da viabilidade de tempo de viagem para R1
                            std::deque<int> r1 (s.R[v1].rOrder);
                            r1.erase(r1.begin() + pos1);
                            r1.erase(r1.begin() + pos1);
                            r1.emplace(r1.begin() + pos1,k);

                            int cr1 = checkRoutCost(d,r1);
                            //int ttt1 = cr1 + d.s[1]*r1.size();
                            int ttt1 = checkRoutTTT(d,r1);
                            if (ttt1 > d.H) continue;

                            cr2 += d.e[v2];
                            cr1 += d.e[v1];

                            //se os tempos sao viaveis e a diferenca na troca reduz o custo
                            double diffCosts = ((cr1 + cr2) - (s.R[v1].rCost + s.R[v2].rCost));
                            //if (ttt1 <= d.H && ttt2 <= d.H && diffCosts < 0){
                            if (ttt1 <= d.H && ttt2 <= d.H && diffCosts < bestDiff){
                                //std::cout << "\n\tMelhoria: (" << i << "," << j << ") trocados do veiculo v1 = " << v1 << " por " << k << " do veiculo v2 = " << v2;

                                melhora = true;
                                melhorou = true;

                                bV1 = v1;
                                bV2 = v2;

                                bestI = i;
                                bestJ = j;
                                bestK = k;

                                bestR = cpRout(r);
                                //atualizando a rota que recebe os clientes adjacentes
                                bestR.z[t][v2] = boost::dynamic_bitset<> (s.R[v2].rComp);
                                bestR.z[t][v2].set(i-1);
                                bestR.z[t][v2].set(j-1);
                                bestR.z[t][v2].reset(k-1);
                                bestV2 = newRoute(cr2, ttt2, s.R[v2].rLoad + loadIJ - r.q[t][k], r2, bestR.z[t][v2]);

                                //atualizando a rota que perde os clientes adjacentes
                                bestR.z[t][v1] = boost::dynamic_bitset<> (s.R[v1].rComp);
                                bestR.z[t][v1].reset(i-1);
                                bestR.z[t][v1].reset(j-1);
                                bestR.z[t][v1].set(k-1);
                                bestV1 = newRoute(cr1, ttt1, s.R[v1].rLoad - loadIJ + r.q[t][k], r1, bestR.z[t][v1]);

                                bestDiff = diffCosts;

                                bestPS = std::vector<std::pair<int,int> > (s.ps);
                                if (bestV1.rOrder.size() >= 1) updatePredSucc(d.n,bestV1.rOrder,bestPS);
                                updatePredSucc(d.n,bestV2.rOrder,bestPS);

                                //getchar();
                            }
                        }
                    }
                }
            }
        }

        if (melhora){
            //atualizando solução
            double oldV1Cost = s.R[bV1].rCost;
            double oldV2Cost = s.R[bV2].rCost;

            /*
            std::cout << yellow << "\n\told cost = " << s.c;
            std::cout << red;
            printRoute(d,s.R[bV1],bV1);
            printRoute(d,s.R[bV2],bV2);
            std::cout << green;
            printRoute(d,bestV1,bV1);
            printRoute(d,bestV2,bV2);
            std::cout << normal;
            */

            s.R[bV1] = cpyRoute(bestV1);
            s.R[bV2] = cpyRoute(bestV2);
            s.ps = std::vector<std::pair<int,int> > (bestPS);
            s.c += bestDiff;
            s.cpp[bestI] = bV2;
            s.cpp[bestJ] = bV2;
            s.cpp[bestK] = bV1;

            r = cpRout(bestR);

            //std::cout << yellow << "\n\tnew cost = " << s.c;
        }
    }while(melhora);

    //if (melhorou && !checkVRSCost(d,r,s,t,"SWAP21")) melhorou = false;
    //std::cout << melhorou;
    return melhorou;
}
//dois clientes adjacentes são trocados por outros dois adjacentes de outra rota
bool swap22(data d, vrs &s, routing &r, int t){
    //std::cout << "\n\t\t\tSWAP22(t="<< t << ") >> ";
    //sinaliza melhoria
    bool melhora = false;
    bool melhorou = false;

    do{
        //armazena a melhor alteracao
        int bV1, bV2, bestI, bestJ, bestK, bestL;
        //std::vector<std::pair<int,int> > bestPS (s.ps);
        std::vector<std::pair<int,int> > bestPS;
        route bestV1, bestV2;
        routing bestR;
        double bestDiff = 0;

        melhora = false;
        for (int v1 = 1; v1 <= d.V; v1++){
            if (s.R[v1].rOrder.size() < 2) continue;

            for (int pos1 = 0; pos1 < s.R[v1].rOrder.size() - 1; pos1++){
                int i = s.R[v1].rOrder[pos1];
                int j = s.R[v1].rOrder[pos1 + 1];
                int loadIJ = r.q[t][i] + r.q[t][j];

                for (int v2 = 1; v2 <= d.V; v2++){
                    if (s.R[v2].rOrder.size() < 2 || v1 == v2) continue;

                    for (int pos2 = 0; pos2 < s.R[v2].rOrder.size() - 1; pos2++){
                        int k = s.R[v2].rOrder[pos2];
                        int l = s.R[v2].rOrder[pos2 + 1];
                        int loadKL = r.q[t][k] + r.q[t][l];

                        //se a troca entre (i,j) e (k,l) eh viavel, entao tentaremos ver se melhora a solução
                        if ((s.R[v2].rLoad - loadKL + loadIJ <= d.Q[v2]) && (s.R[v1].rLoad + loadKL - loadIJ <= d.Q[v1])){
                            //teste da viabilidade de tempo de viagem para R2
                            std::deque<int> r2 (s.R[v2].rOrder);
                            r2.erase(r2.begin() + pos2);
                            r2.erase(r2.begin() + pos2);
                            r2.emplace(r2.begin()+pos2,i);
                            r2.emplace(r2.begin()+pos2 + 1,j);

                            int cr2 = checkRoutCost(d,r2);
                            //int ttt2 = cr2 + d.s[1]*r2.size();
                            int ttt2 = checkRoutTTT(d,r2);
                            if (ttt2 > d.H) continue;

                            //teste da viabilidade de tempo de viagem para R1
                            std::deque<int> r1 (s.R[v1].rOrder);
                            r1.erase(r1.begin() + pos1);
                            r1.erase(r1.begin() + pos1);
                            r1.emplace(r1.begin()+pos1,k);
                            r1.emplace(r1.begin()+pos1 + 1,l);

                            int cr1 = checkRoutCost(d,r1);
                            //int ttt1 = cr1 + d.s[1]*r1.size();
                            int ttt1 = checkRoutTTT(d,r1);
                            if (ttt1 > d.H) continue;

                            cr2 += d.e[v2];
                            cr1 += d.e[v1];

                            //se os tempos sao viaveis e a diferenca na troca reduz o custo
                            double diffCosts = ((cr1 + cr2) - (s.R[v1].rCost + s.R[v2].rCost));
                            //if (ttt1 <= d.H && ttt2 <= d.H &&  diffCosts < 0){
                            if (ttt1 <= d.H && ttt2 <= d.H &&  diffCosts < bestDiff){
                                //std::cout << "\n\tMelhoria: (" << i << "," << j << ") trocados do veiculo v1 = " << v1 << " por (" << k << "," << l << ") do veiculo v2 = " << v2;

                                melhora = true;
                                melhorou = true;

                                bV1 = v1;
                                bV2 = v2;

                                bestI = i;
                                bestJ = j;
                                bestK = k;
                                bestL = l;

                                bestR = cpRout(r);
                                //atualizando a rota que recebe os clientes adjacentes
                                bestR.z[t][v2] = boost::dynamic_bitset<> (s.R[v2].rComp);
                                bestR.z[t][v2].set(i-1);
                                bestR.z[t][v2].set(j-1);
                                bestR.z[t][v2].reset(k-1);
                                bestR.z[t][v2].reset(l-1);
                                bestV2 = newRoute(cr2, ttt2, s.R[v2].rLoad + loadIJ - loadKL, r2, bestR.z[t][v2]);

                                //atualizando a rota que perde os clientes adjacentes
                                bestR.z[t][v1] = boost::dynamic_bitset<> (s.R[v1].rComp);
                                bestR.z[t][v1].reset(i-1);
                                bestR.z[t][v1].reset(j-1);
                                bestR.z[t][v1].set(k-1);
                                bestR.z[t][v1].set(l-1);
                                bestV1 = newRoute(cr1, ttt1, s.R[v1].rLoad - loadIJ + loadKL, r1, bestR.z[t][v1]);

                                bestDiff = diffCosts;

                                bestPS = std::vector<std::pair<int,int> > (s.ps);
                                if (bestV1.rOrder.size() > 0) updatePredSucc(d.n,bestV1.rOrder,bestPS);
                                if (bestV2.rOrder.size() > 0) updatePredSucc(d.n,bestV2.rOrder,bestPS);

                                //getchar();
                            }
                        }
                    }
                }
            }
        }

        if (melhora){
            //atualizando solução
            double oldV1Cost = s.R[bV1].rCost;
            double oldV2Cost = s.R[bV2].rCost;

            /*
            std::cout << yellow << "\n\told cost = " << s.c;
            std::cout << red;
            printRoute(d,s.R[bV1],bV1);
            printRoute(d,s.R[bV2],bV2);
            std::cout << green;
            printRoute(d,bestV1,bV1);
            printRoute(d,bestV2,bV2);
            std::cout << normal;
            */

            s.R[bV1] = cpyRoute(bestV1);
            s.R[bV2] = cpyRoute(bestV2);
            s.ps = std::vector<std::pair<int,int> > (bestPS);
            s.c += bestDiff;
            s.cpp[bestI] = bV2;
            s.cpp[bestJ] = bV2;
            s.cpp[bestK] = bV1;
            s.cpp[bestL] = bV1;

            r = cpRout(bestR);

            //std::cout << yellow << "\n\tnew cost = " << s.c;
        }
    }while(melhora);

    //if (melhorou && !checkVRSCost(d,r,s,t,"SWAP22")) melhorou = false;
    //std::cout << melhorou;
    return melhorou;
}
bool interRouteOpt(data d, routing &r, vrs &s, int t){
    //std::cout << "\n\t\t\tinter-Route procedure[t=" << t << "] >> ";    
    /*
    std::cout << YELLOW;
    printVRS(d,r,s,t);
    std::cout << normal;
    //*/
    double currSolValue = s.c;
    bool melhorou = false;
    vrs sl = cpyVRS(s);
    routing rl = cpRout(r);

    if (r.V[t] == 1)
        melhorou = intraRouteOpt(d, std::vector<bool> (d.V + 1, true), sl);
    else if (r.V[t] > 1){
        // sinaliza se a busca trouxe melhoria
        //(1) VLNS + Solo, (2) Shift20, (3) Swap21, (4) Swap22, (5) Cross, (6,7,8) ShiftK, K = 3,4,5
        std::vector<int> success (9,1);
        success[0] = 8;

        //enquanto ainda existe alguma busca capaz de melhorar a rota v
        std::uniform_int_distribution<> LS(1,8);
        while(success[0] > 0){
            int ls = LS(gen);
            while(success[ls] == 0) ls = LS(gen);
            bool melhora = false;

            std::vector<route> rota (sl.R);

            if (ls == 1) //(1) VLNS + Solo
                melhora = vlns(d,rl,t,sl);            
            else if (ls == 2) //(2) Shift20
                melhora = shift20(d,sl,rl,t);            
            else if (ls == 3) //(3) Swap21
                melhora = swap21(d,sl,rl,t);
            else if (ls == 4) //(4) Swap22
                melhora = swap22(d,sl,rl,t);
            else if (ls == 5) //(5) Cross
                melhora = cross(d,sl,rl,t);
            else if (ls == 6) //(6) Shift, K = 3
                melhora = shift(d,sl,rl,t,3);
            else if (ls == 7) //(7) Shift, K = 4
                melhora = shift(d,sl,rl,t,4);
            else if (ls == 8) //(8) Shift, K = 5
                melhora = shift(d,sl,rl,t,5);

            if (melhora){
                //std::cout << yellow << "\n\t\t\tnew VRS cost = " << sl.c << normal;
                /*
                std::cout << BLUE;
                printVRS(d,r,sl,t);
                std::cout << normal;
                //*/
                std::vector<bool> mudou (d.V + 1, false);
                for (int v = 1; v <= d.V; v++) if (rota[v].rCost != sl.R[v].rCost) mudou[v] = true;
                melhora = intraRouteOpt(d,mudou,sl);
                /*
                std::cout << GREEN;
                printVRS(d,r,sl,t);
                std::cout << normal;
                //*/
                //std::cout << yellow << "\n\t\t\tnew VRS cost = " << sl.c << normal;
                //*/
                //std::cout << 1 << " ";
                melhorou = true;
                s = cpyVRS(sl);
                r = cpRout(rl);
                iers[ls]++;
                iers[0]++;
            }else{
                success[ls] = 0;
                success[0]--;
            }
        }
    }

    if (melhorou){
        s = cpyVRS(sl);
        r = cpRout(rl);
        /*
        std::cout << PURPLE;
        printVRS(d,r,s,t);
        std::cout << normal;
        //*/
        return true;
    }else
        return false;
}

//INTRA-ROUTE HEURISTICS
//One-point Move
double OnePnt(data d, std::vector<std::pair<int,int> > &ps, route &r){
    //std::cout << "\n\t\t\t\tOne-point >> cost before = " << r.rCost << " >> cost after = ";
    bool improved;

    do{
        improved = false;
        int bP1 = -1;
        int bP2 = -1; //best positions to be changed
        double bestCost = r.rCost;
        double currCost = r.rCost;
        int currTTT = r.rTTT;
        int bestTTT = r.rTTT;
        for (int p1 = 0; p1 < r.rOrder.size() - 2; p1++){
            int i = r.rOrder[p1];
            for (int p2 = p1 + 1; p2 < r.rOrder.size() - 1; p2++){
                int j = r.rOrder[p2];
                double newCost = 999999;
                double newTTT = 999999;
                if (p1 == p2-1){
                    newCost = currCost - (d.c[ps[i].first][i] + d.c[j][ps[j].second]) + (d.c[ps[i].first][j] + d.c[i][ps[j].second]);
                    newTTT = currTTT - (d.a[ps[i].first][i] + d.a[j][ps[j].second]) + (d.a[ps[i].first][j] + d.a[i][ps[j].second]);
                }else{
                    newCost = currCost - (d.c[ps[i].first][i] + d.c[i][ps[i].second] + d.c[j][ps[j].second]) + (d.c[ps[i].first][ps[i].second] + d.c[j][i] + d.c[i][ps[j].second]);
                    newTTT = currTTT - (d.a[ps[i].first][i] + d.a[i][ps[i].second] + d.a[j][ps[j].second]) + (d.a[ps[i].first][ps[i].second] + d.a[j][i] + d.a[i][ps[j].second]);
                }

                if (newCost < bestCost && newTTT <= d.H){
                    bestCost = newCost;
                    bestTTT = newTTT;
                    bP1 = p1;
                    bP2 = p2;
                    improved = true;
                }
            }
        }

        if (improved){
            //houve melhora
            int i = r.rOrder[bP1];
            int j = r.rOrder[bP2];

            //atualizando
            r.rOrder.emplace(r.rOrder.begin() + bP2 + 1,i);
            r.rOrder.erase(r.rOrder.begin() + bP1);

            r.rCost = bestCost;
            r.rTTT = bestTTT;
            r.rFST = r.rOrder.front();
            r.rLST = r.rOrder.back();
            updatePredSucc(d.n, r.rOrder, ps);
        }
    }while(improved);

    //std::cout << r.rCost;
    return r.rCost;
}
//Two-point Move
double TwoPnt(data d, std::vector<std::pair<int,int> > &ps, route &r){
    //std::cout << "\n\t\t\t\tTwo-point >> cost before = " << r.rCost << " >> cost after = ";
    bool improved;

    std::vector<boost::dynamic_bitset<> > notAllowedTrade (d.n + 1, boost::dynamic_bitset<> (d.n,0));

    do{
        improved = false;
        int bP1 = -1;
        int bP2 = -1; //best positions to be changed
        double bestCost = r.rCost;
        double currCost = r.rCost;
        int currTTT = r.rTTT;
        int bestTTT = r.rTTT;
        for (int p1 = 0; p1 < r.rOrder.size() - 2; p1++){
            int i = r.rOrder[p1];
            for (int p2 = p1 + 1; p2 < r.rOrder.size() - 1; p2++){
                int j = r.rOrder[p2];
                if (!notAllowedTrade[i].test(j-1)){
                    double newCost = 999999;
                    double newTTT = 999999;
                    if (p1 == p2-1){
                        newCost = currCost - (d.c[ps[i].first][i] + d.c[j][ps[j].second]) + (d.c[ps[i].first][j] + d.c[i][ps[j].second]);
                        newTTT = currTTT - (d.a[ps[i].first][i] + d.a[j][ps[j].second]) + (d.a[ps[i].first][j] + d.a[i][ps[j].second]);
                    }
                    else{
                        newCost = currCost - (d.c[ps[i].first][i] + d.c[i][ps[i].second] + d.c[ps[j].first][j] + d.c[j][ps[j].second]) + (d.c[ps[i].first][j] + d.c[j][ps[i].second] + d.c[ps[j].first][i] + d.c[i][ps[j].second]);
                        newTTT = currTTT - (d.a[ps[i].first][i] + d.a[i][ps[i].second] + d.a[ps[j].first][j] + d.a[j][ps[j].second]) + (d.a[ps[i].first][j] + d.a[j][ps[i].second] + d.a[ps[j].first][i] + d.a[i][ps[j].second]);
                    }

                    if (newCost < bestCost && newTTT <= d.H){
                        bestCost = newCost;
                        bestTTT = newTTT;
                        bP1 = p1;
                        bP2 = p2;
                        improved = true;
                    }
                }
            }
        }

        if (improved){
            //houve melhora
            int i = r.rOrder[bP1];
            int j = r.rOrder[bP2];

            notAllowedTrade[i].set(j-1);
            notAllowedTrade[j].set(i-1);

            //atualizando
            int aux = r.rOrder[bP1];
            r.rOrder[bP1] = r.rOrder[bP2];
            r.rOrder[bP2] = aux;

            r.rCost = bestCost;
            r.rTTT = bestTTT;
            r.rFST = r.rOrder.front();
            r.rLST = r.rOrder.back();
            updatePredSucc(d.n, r.rOrder, ps);
        }
        //getchar();
    }while(improved);

    //std::cout << r.rCost;
    return r.rCost;
}
//Three-point Move
double ThrPnt(data d, std::vector<std::pair<int,int> > &ps, route &r, int v){
    //std::cout << "\n\t\t\t\tThree-point >> cost before = " << r.rCost << " >> cost after = ";
    std::deque<int> rS = std::deque<int> (r.rOrder);
    double bCost = r.rCost;
    double bTTT = r.rTTT;
    std::deque<int> bR = std::deque<int> (r.rOrder);
    bool continua = false;
    do{
        bool melhoria = false;
        for (int p1 = 0; p1 < rS.size() - 2; p1++){
            int i = rS[p1];
            int p2 = p1 + 1;
            int j = rS[p2];
            for (int p3 = p1 + 2; p3 < rS.size(); p3++){
                std::deque<int> r_ = std::deque<int> (rS);
                int k = rS[p3];

                if (p3 == p1 + 2){
                    r_[p1] = k;
                    r_[p2] = i;
                    r_[p3] = j;
                }else if (p3 > p1 + 2){
                    for (int pAux = p2; pAux < p3 - 1; pAux++) r_[pAux] = r_[pAux+1];
                    r_[p1] = k;
                    r_[p3] = j;
                    r_[p3-1] = i;
                }

                double newCost = d.e[v] + checkRoutCost(d,r_);
                int newTTT = (newCost - d.e[v] + r_.size()*d.s[1]);

                //se o custo da nova rota eh melhor e ela respeita o limite de viagem
                if (newCost < bCost && newTTT <= d.H){
                    bCost = newCost;
                    bTTT = newTTT;
                    bR = std::deque<int> (r_);
                    melhoria = true;
                }
            }
        }

        if (melhoria) {
            rS = std::deque<int> (bR);
            continua = false;
        }
    }while(continua);

    r.rOrder = std::deque<int> (bR);
    r.rTTT = bTTT;
    r.rFST = r.rOrder.front();
    r.rLST = r.rOrder.back();
    r.rCost = bCost;
    updatePredSucc(d.n, r.rOrder, ps);

    //std::cout << r.rCost;
    return r.rCost;
}
//Two-opt move
double TwoOpt(data d, std::vector<std::pair<int,int> > &ps, route &r){
    //std::cout << "\n\t\t\t\tTwo-opt >> cost before = " << r.rCost << " >> cost after = ";
    bool improved;
    do{
        improved = false;
        int bP1 = -1;
        int bP2 = -1; //best positions to be changed
        double bestCost = r.rCost;
        double currCost = r.rCost;
        int currTTT = r.rTTT;
        int bestTTT = r.rTTT;

        for (int p1 = 0; p1 < r.rOrder.size() - 2; p1++){
            for (int p2 = p1 + 1; p2 < r.rOrder.size() - 1; p2++){
                int pred = ps[r.rOrder[p1]].first;
                int succ = ps[r.rOrder[p2]].second;
                int i = r.rOrder[p1];
                int j = r.rOrder[p2];

                int newTTT = r.rTTT - d.a[pred][i] - d.a[j][succ] + d.a[pred][j] + d.a[i][succ];
                double newCost = r.rCost - d.c[pred][i] - d.c[j][succ] + d.c[pred][j] + d.c[i][succ];

                if (newCost < bestCost && newTTT <= d.H){
                    bestCost = newCost;
                    bestTTT = newTTT;
                    bP1 = p1;
                    bP2 = p2;
                    improved = true;
                }
            }
        }

        if (improved){
            //houve melhora
            int i = r.rOrder[bP1];
            int j = r.rOrder[bP2];

            std::deque<int> rOrder;
            for (int pp = 0; pp < bP1; pp++)
                rOrder.push_back(r.rOrder[pp]);

            for (int pp = bP2; pp >= bP1; pp--)
                rOrder.push_back(r.rOrder[pp]);

            for (int pp = bP2 + 1; pp < r.rOrder.size(); pp++)
                rOrder.push_back(r.rOrder[pp]);

            r.rOrder = rOrder;
            r.rTTT = bestTTT;
            r.rCost = bestCost;
            r.rFST = r.rOrder.front();
            r.rLST = r.rOrder.back();
            updatePredSucc(d.n, r.rOrder, ps);
        }
    }while(improved);

    //std::cout << r.rCost;
    return r.rCost;
}
//Or-opt move
double OrOpt(data d, std::vector<std::pair<int,int> > &ps, route &r, int Or, int v){
    //std::cout << "\n\t\t\t\tOr-opt, Or = " << Or << " >> cost before = " << r.rCost << " >> cost after = ";
    std::deque<int> rS = std::deque<int> (r.rOrder);
    double bCost = r.rCost;
    double bTTT = r.rTTT;
    std::deque<int> bR = std::deque<int> (r.rOrder);
    bool continua = false;
    do{
        bool melhoria = false;
        for (int p1 = 0; p1 < rS.size() - Or; p1++){
            std::deque<int> seq;
            for (int ind = p1; ind < p1 + Or; ind++) seq.push_back(rS[ind]);

            for (int p2 = p1 + Or; p2 < rS.size(); p2++){
                std::deque<int> seq2 (seq);
                if (p1 > 0){
                    for (int aux = p2; aux > p1 + Or - 1; aux--) seq2.emplace(seq2.begin(),rS[aux]);
                    for (int aux = p1 - 1; aux > -1; aux--) seq2.push_front(rS[aux]);
                    for (int aux = p2 + 1; aux < rS.size(); aux++) seq2.push_back(rS[aux]);
                }else{
                    for (int aux = p2; aux > p1 + Or - 1; aux--) seq2.emplace(seq2.begin(),rS[aux]);
                    for (int aux = p2 + 1; aux < rS.size(); aux++) seq2.push_back(rS[aux]);
                }

                double newCost = d.e[v] + checkRoutCost(d,seq2);
                int newTTT = (newCost - d.e[v] + seq2.size()*d.s[1]);
                //se o custo da nova rota eh melhor e ela respeita o limite de viagem
                if (newCost < bCost && newTTT <= d.H){
                    bCost = newCost;
                    bTTT = newTTT;
                    bR = std::deque<int> (seq2);
                    melhoria = true;
                }
            }
        }

        if (melhoria) {
            rS = std::deque<int> (bR);
            continua = false;
        }
    }while(continua);

    r.rOrder = std::deque<int> (bR);
    r.rTTT = bTTT;
    r.rFST = r.rOrder.front();
    r.rLST = r.rOrder.back();
    r.rCost = bCost;
    updatePredSucc(d.n, r.rOrder, ps);

    //std::cout << r.rCost;
    return r.rCost;
}
bool intraRouteOpt(data d, std::vector<bool> mudou, vrs &s){
    //std::cout << "\n\t\t\tintra-Route procedure";
    double currSolCost = s.c;
    //bool melhora = false;
    bool melhorou = false;
    
    for (int v = 1; v <= d.V; v++){
        if (s.R[v].rOrder.size() <= 2 || !mudou[v]) continue;
        
        //sinaliza se a busca trouxe melhoria
        //(1) OnePnt, (2) TwoPnt, (3) TwoOpt, (4) ThrPnt, (5,6,7) Or-opt, Or = 2,3,4
        std::vector<int> success (8,1);
        success[0] = 7;
        
        //nao aplico Or
        if (s.R[v].rOrder.size() <= 4){
            for (int ls = 5; ls <= 7; ls++) {
                success[ls] = 0; //
                success[0]--;                
            }            
        }
        
        //nao aplico ThrPnt
        if (s.R[v].rOrder.size() <= 3){
            success[4] = 0;
            success[0]--;            
        }
        
        //nao aplico nenhuma
        if (s.R[v].rOrder.size() <= 2) success[0] = 0;
        if (success[0] == 0) continue;
        
        //std::cout << "\n\t\t\t\tv = " << v;
        //enquanto ainda existe alguma busca capaz de melhorar a rota v
        std::uniform_int_distribution<> LS(1,7);
        while(success[0] > 0){
            int ls = LS(gen);
            while(success[ls] == 0) ls = LS(gen);
            
            double currRouteCost = s.R[v].rCost;
            double newRouteCost = s.R[v].rCost;
            
            route rl = cpyRoute(s.R[v]);
            std::vector<std::pair<int,int> > psl (s.ps);
            
            if (ls == 1) //1-opt Move
                newRouteCost = OnePnt(d, psl, rl);
            else if (ls == 2) //2-pnt Move
                newRouteCost = TwoPnt(d, psl, rl);
            else if (ls == 3) //2-opt Move
                newRouteCost = TwoOpt(d, psl, rl);
            else if (ls == 4) //3-pnt Move
                newRouteCost = ThrPnt(d, psl, rl,v);
            else if (ls == 5) //Or-opt Move
                newRouteCost = OrOpt(d, psl, rl, 3, v);
            else if (ls == 6) //Or-opt Move
                newRouteCost = OrOpt(d, psl, rl, 4, v);
            else if (ls == 7) //Or-opt Move
                newRouteCost = OrOpt(d, psl, rl, 5, v);
            
            //VERIFICAÇÃO DE ERRO
            /*
            //double cost = d.e[v] + checkRoutCost(d,s.R[v].rOrder);
            //if (cost != s.R[v].rCost){
                    std::cout << RED << "\n\t\tdeu zica nos custos de roteamento [v=" << v << "]" << " - cost: " << cost << " <> " << " - wrong: " << rl.rCost << normal;

                    std::cout << green;
                    printRoute(d,rl,v);
                    std::cout << red;
                    printRoute(d,s.R[v],v);
                    std::cout << normal;
                    for (int pos = 0; pos < s.R[v].rOrder.size(); pos++){
                        int i = s.R[v].rOrder[pos];
                        int pi = psl[i].first;
                        int si = psl[i].second;
                        std::cout << red << "\n\t\t pred = " << pi << " >> i = " << i << " >> succ = " << si;
                        pi = s.ps[i].first;
                        si = s.ps[i].second;
                        std::cout << green << "\n\t\t pred = " << pi << " >> i = " << i << " >> succ = " << si << normal;
                    }
                    help();
            }
            double cost = d.e[v] + checkRoutCost(d,rl.rOrder);
            if (cost != rl.rCost) newRouteCost = 9999999;
            //*/
            
            //se houve melhoria, atualizo o custo da solução
            if (newRouteCost < currRouteCost){
                s.R[v] = cpyRoute(rl);
                s.ps = std::vector<std::pair<int,int> > (psl);
                s.c -= currRouteCost;
                s.c += newRouteCost;
                //printRoute(d,s.R[v],v);
                iars[ls]++;
                iars[0]++;
            } else{ //caso contrario, removo a busca da lista de candidatas
                success[ls] = 0;
                success[0]--;
            }
        }
    }
    
    if (s.c < currSolCost) melhorou = true;
                 
    double costVRS = 0;
    for (int v = 1; v <= d.V; v++)
        if (s.R[v].rComp.any()) costVRS += (d.e[v] + checkRoutCost(d,s.R[v].rOrder));

    if (costVRS != s.c)  {
        /*
        std::cout << RED << "\n\t\tAlgo errado na INTRAROUTE, custo = " << s.c << " <> custo real = " << costVRS;
        for (int v = 1; v <= d.V; v++)
            if (s.R[v].rComp.any())
                std::cout << "\n\t\t\tv = " << v << " >> custo = " << s.R[v].rCost << " <> real = " << (d.e[v] + checkRoutCost(d,s.R[v].rOrder));
        //*/
        melhorou = false;
        //help();
    }

    return melhorou;
}
bool BuildAndOptimizeVRS(data d, std::vector<vrs> &vrp, routing &r, double &totVRPcost){
    //std::cout << "\n\t\t\tConstruindo e otimizando solução do VRP";
    bool failed = false;
    for (int t = 1; t <= d.T && !failed; t++){
        vrp[t] = S0(d,t,r,failed);                
        if (!failed){
            if (vrp[t].anyRoute) 
                interRouteOpt(d,r,vrp[t],t);                
            
            totVRPcost += vrp[t].c;
            /*
            std::cout << RED;
            printVRS(d,r,vrp[t],t);
            std::cout << normal;
            //*/
        }else{
            totVRPcost += 999999;
            return false;
        }
    }


    if (failed) return false;
    else return true;
}
//end local search routing methods

//perturbacao do PPP
mov newMove(int k, int p, int n){
    mov m;    
    m.k = k;
    m.p = p;
    m.n = n;
    return m;
}
//transferir entregas de produtos k dos periodos que nao tem setup para periodos que tem ou proximos
std::vector<std::vector<std::vector<mov> > > genTransfProdMoves (data d, prodinv p, std::vector<vrs> vrp, int &n){
    //std::cout << RED << "\n\tGERANDO MOVIMENTOS VIÁVEIS" << normal;
    std::vector<std::vector<std::vector<mov> > > moves (d.T + 1, std::vector<std::vector<mov> > (d.T + 1, std::vector<mov> ()));
            
    //so removerei cargas de periodos que nao ocorreram setup produtivo
    std::vector<boost::dynamic_bitset<> > setup (d. P + 1, boost::dynamic_bitset<> (d.T, 0));
    for (int k = 1; k <= d.P; k++){
        //std::cout << "\nk = " << k;
        for (int t = 1; t <= d.T; t++){
            if (p.y[t][k] == 1) setup[k].set(t-1);
            //std::cout << "\t y[" << t << "] = " << p.y[t][k] << " >> " << p.p[t][k] << "/" << d.M_[t][k] << "/" << d.C[k];
        }
    }
    //std::cout << "\n";
    //essas cargas podem ser transferidas para outros periodos, diferentes do corrente t1 e menores do que o proximo onde também ocorre setup
    n = 0;
    for (int k = 1; k <= d.P; k++){
        for (int t1 = 1; t1 <= d.T; t1++){
            //se nao existe setup em t1, nao ha o que transferir
            if (setup[k][t1-1] == 0) continue;
            
            //tl e tll sao os limites inferior e superior para t2, estagios
            int tl = 0;
            int tll = 0;
            /*
            for (int t = 1; t <= d.T; t++) {
                if (t < t1 && tl == 0 && setup[k][t-1] == 1) tl = t;
                if (t > t1 && tll == 0 && setup[k][t-1] == 1) tll = t;                                            
            }
            */
            for (int t = t1 - 1; t >= 1; t--) 
                if (tl == 0 && setup[k][t-1] == 1) tl = t;
            
            for (int t = t1 + 1; t <= d.T; t++) 
                if (tll == 0 && setup[k][t-1] == 1) tll = t;                                            
                        
            if (tl == 0) tl = t1;
            if (tll == 0) tll = t1;            
            //std::cout << "\n\ttl/t1/tll = " << tl << "/" << t1 << "/" << tll;
            
            for (int t2 = 1; t2 <= d.T; t2++){
                if (t2 != t1 && t2 >= tl && t2 <= tll){
                    int o = 99999;
                    o = std::min(int(p.p[t1][k]), int(d.M_[t2][k] - p.p[t2][k]));
                    o = std::min(o,int(d.U[k][0] - p.I[t2][k][0]));
                    //std::cout << "\n\t1 o = " << p;
                    if (o == 0) continue;
                    if (o > 0 && o < 99999) {
                        //std::cout << red << "\n\t(" << t1 << ">>" << t2 << "), k=" << k << ", p = " << o << normal;
                        moves[t1][t2].push_back(newMove(k,o,n));
                        n++;                        
                    }                    
                }                
            }            
        }        
    }
    
    return moves;
}
bool modifyProdPlan(data d, prodinv &p, std::vector<vrs> vrp, cpx &ppp){
    
    //std::cout << BLUE << "\n\tCHANGING DELIVERY PLAN" << normal;	
    int totMov = 0;
    std::vector<std::vector<std::vector<mov> > > transfMoves (genTransfProdMoves (d, p, vrp, totMov));
    //std::cout << red <<  "\n\ttotMov = " << totMov;
    //exit(1);
    
    //getchar();    
    if (totMov > 0){        
        //sinaliza se já teve alteracao neste produto
        std::vector<bool> P (d.P + 1, false);
        //sinaliza se o movimento jah foi usado
        std::vector<bool> usedMove (totMov + 1, false);
                
        //sorteando os movimentos
        std::uniform_int_distribution<> T1T2(1,d.T);
        std::vector<mov> moves;
        std::vector<std::pair<int,int> > periods;
        steady_clock::time_point initialTime = steady_clock::now();
        duration<double> diffTime = duration_cast<duration<double> > (steady_clock::now()-initialTime);
        do{
            int t1 = T1T2(gen);                
            int t2 = T1T2(gen);
            while(t2 == t1) t2 = T1T2(gen);
            //std::cout << "\n\tt1 = " << t1 << " >> t2 = " << t2 << " - |t1.t2| = " << transfMoves[t1][t2].size();
            //getchar();
            if (transfMoves[t1][t2].size() > 0){
                std::uniform_int_distribution<> K1(0,transfMoves[t1][t2].size() - 1);
                mov m = transfMoves[t1][t2][K1(gen)];
                if (usedMove[m.n]) continue;
                moves.push_back(m);
                periods.push_back(std::pair<int,int> (t1,t2));
                //std::cout << blue << "\n\tm = (" << t1 << ">>" << t2 << "), k=" << m.k << ", o = " << m.o << normal;
                usedMove[m.n] = true;
                P[m.k] = true;
            }
            diffTime = duration_cast<duration<double> > (steady_clock::now()-initialTime);
        }while (diffTime <= segundo);    
                        
        std::uniform_int_distribution<> M(0,periods.size() - 1);
        std::vector<bool> K (d.P + 1, false);
        std::vector<bool> T (d.T + 1, false);        
        int cntMov = 0;
        //int maxNumbMov = std::round(moves.size()/3);
        //int maxNumbMov = std::round(d.T/3);
        int maxTest = moves.size();
        
        while(cntMov <= maxPertP){
            //std::cout << GREEN << "\n\taqui 1";
            int mv = M(gen);
            int t1 = periods[mv].first;
            int t2 = periods[mv].second;            
            int k = moves[mv].k;
            int o = moves[mv].p;
            
            //std::cout << blue << "\n\tm1 = (" << t1 << ">>" << t2 << "), k=" << k << ", p = " << o << normal;
                        
            while (T[t1] || T[t2] || K[k]) {
                //std::cout << YELLOW << "\n\t\taqui 2";
                mv = M(gen);
                maxTest--;
                if (maxTest == 0)break;
            }
            if (maxTest == 0) break;
                                                
            T[t1]  = true;
            T[t2]  = true;
            K[k] = true;
            cntMov++;
            
            try{
                //periodo que perde parte do lote de produção
                ppp.p[d.P*(t1-1) + (k-1)].setUB(p.p[t1][k] - o);
                p.p[t1][k] -= o;
                //periodo que ganha parte do lote de produção
                ppp.p[d.P*(t2-1) + (k-1)].setLB(p.p[t2][k] + o);
                p.p[t2][k] += o;
                
				totalOfPerturbations++;
            }
            catch(IloException& ex){
                //std::cerr << RED << "\nError: " << ex << std::endl;
                return false;
            }
        }
        //ppp.cplx.exportModel("ModelWithPerturbations.lp");    
        return true; 
        //std::cout << BLUE << "\n\tEND CHANGING DELIVERY PLAN" << normal;
    }else return false;
    //*/
}
bool correctBoundsProdVars(data d, cpx &ppp){
    try{
        for (int t = 1; t <= d.T; t++)
            for (int k = 1; k <= d.P; k++)
                ppp.p[d.P*(t-1) + (k-1)].setBounds(0, d.M_[t][k]);
            
        //ppp.cplx.exportModel("ModelWithCorrectionss.lp");
        
        return true;
    }catch(IloException& ex){
        //std::cerr << RED << "\nError: " << ex << std::endl;
        return false;
    }    
}

//Tactical ILS
double tils(data d){
    //std::cout << BLUE << "\n\t\tTactical ILS" << normal;
    
    //solution sStar/s0    
    routing rs = newRout(d);
    prodinv ps = newProd(d);
    std::vector<vrs> vrps (d.T + 1, newVRS(d.n,d.V));
    std::vector<double> fos (4,0);
        
    steady_clock::time_point tStarIni = steady_clock::now();
    steady_clock::time_point tInitial = steady_clock::now();
    //resolvendo o nivel tático
    cpx ppp;    
    bool pppFeas = buildAndSolvePPP(d,ppp);
    if (pppFeas){
        pppRecovInfo(d,ppp,rs,ps,fos);        
    }else
        return 999999;
        
    steady_clock::time_point tFinal = steady_clock::now(); //retorna o ponto de agora no tempo
    duration<double> diffTime = duration_cast<duration<double> > (tFinal-tInitial);
    double timePPP = double(diffTime.count());
            
    bool vrpFeas = false;
    tInitial = steady_clock::now();
    if (pppFeas) vrpFeas = BuildAndOptimizeVRS(d,vrps,rs,fos[3]);
    tFinal = steady_clock::now();
    diffTime = duration_cast<duration<double> > (tFinal-tInitial);        
    double timeVRP = double(diffTime.count());    
    double timeStar = timePPP + timeVRP;
    if (vrpFeas) fos[0] += fos[3];        
    
    //armazena o tempo gasto na otimização de cada um dos níveis
    double totVRPtime = timeVRP;
    double totPPPtime = timePPP;
    
    while (!(vrpFeas && pppFeas)){
         tInitial = steady_clock::now();
         pppFeas = buildAndSolvePPP(d,ppp);
         if (pppFeas)
             pppRecovInfo(d,ppp,rs,ps,fos);
         else
             return 999999;
            
        tFinal = steady_clock::now(); //retorna o ponto de agora no tempo
        diffTime = duration_cast<duration<double> > (tFinal-tInitial);
        timePPP = double(diffTime.count());
        
        //VRP
        tInitial = steady_clock::now();
        vrpFeas = false;
        if (pppFeas) vrpFeas = BuildAndOptimizeVRS(d,vrps,rs,fos[3]);
        tFinal = steady_clock::now();
        diffTime = duration_cast<duration<double> > (tFinal-tInitial);
        timeVRP = double(diffTime.count());
        
        timeStar = timePPP + timeVRP;
        if (vrpFeas) fos[0] += fos[3];
        //std::cout << "\n\tInviavel.";
        diffTime = duration_cast<duration<double> > (tFinal-tStarIni);
        if (diffTime > cincomin){
            std::cout << d.name << " >> ";
            return 999999;            
        }
    }
    
    bool setupsFixed = fixSetups(d,true,ps,ppp);
    if (!setupsFixed) exit(1);    
    
    int iter = 0;
    do{
        iter++;
        //solucao corrente (s)
        routing r = cpRout(rs);
        prodinv p = cpProd(ps);
        std::vector<vrs> vrp (vrps);
        std::vector<double> fo (fos);
        double currTime;
        
        //enquanto o número máximo de iterações do ILS não foi atingido
        int iterILS = 0;        
        while(iterILS <= maxIterILS){            
            //solucao vizinha (s')
            routing rl = cpRout(r);
            prodinv pl = cpProd(p);
            std::vector<vrs> vrpl (vrp);
            std::vector<double> fol (fo);
            
            //perturba PPP            
            bool perturbou = modifyProdPlan(d,pl,vrpl, ppp);
                        
			//atualizo delta(i,k,v,t), levo para o PPP
            updateDelta(d,vrpl,pl,ppp);        //AQUI DELTA DEVE CONSIDERAR A MUDANÇA NO PLANO DE PRODUCAO                
            
            //verifico se houve melhora otimizando
            //PPP
            pppFeas = false;
            tInitial = steady_clock::now();
            try{
                //resolvendo        
                IloTimer crono(ppp.cplx.getEnv());
                ppp.crono = &crono;
                ppp.cplx.solve();        
                ppp.crono->stop();

                //recuperando informacoes
                if(ppp.cplx.getStatus() == IloAlgorithm::Optimal || ppp.cplx.getStatus() == IloAlgorithm::Feasible) {
                    getInfoPPP(d,ppp);
                    pppRecovInfo(d,ppp,rl,pl,fol);
                    pppFeas = true;
                }
            }catch(IloException& ex){
                //std::cerr << RED << "\n\t\tError: " << ex << std::endl;
                pppFeas = false;
            }
            tFinal = steady_clock::now(); //retorna o ponto de agora no tempo
            diffTime = duration_cast<duration<double> > (tFinal-tInitial);
            timePPP = double(diffTime.count());
            totPPPtime += timePPP;
            
            //VRP
            tInitial = steady_clock::now();
            vrpFeas = false;
            if (pppFeas) vrpFeas = BuildAndOptimizeVRS(d,vrpl,rl,fol[3]);
            tFinal = steady_clock::now();
            diffTime = duration_cast<duration<double> > (tFinal-tInitial);        
            timeVRP = double(diffTime.count());
            totVRPtime += timeVRP;
            //timeStar = timePPP + timeVRP;                                    
            if (vrpFeas) fol[0] += fol[3];
            
            if(perturbou){
                //std::cout << "\npertubou = " << perturbou;
				bool corrigiu = correctBoundsProdVars(d, ppp);
				if (!corrigiu) help();
			}
                        
            bool modSetup = false;
            if (fol[0] < fo[0] && vrpFeas && pppFeas){
                //salvando o momento que encontrou uma solução melhor
				tFinal = steady_clock::now();
				diffTime = duration_cast<duration<double> > (tFinal-tStarIni);
				currTime = double(diffTime.count());
                
				//atualizando a nova melhor solução
                r = cpRout(rl);
                p = cpProd(pl);
                vrp = std::vector<vrs> (vrpl);
                fo = std::vector<double> (fol);
                				
                //atualizando parametros				
                if (fo[0] < fos[0]){
					iterILS = 0;
                    //atualizando a nova melhor solução
                    rs = cpRout(r);
                    ps = cpProd(p);
                    vrps = std::vector<vrs> (vrp);
                    fos = std::vector<double> (fo);
                    timeStar = currTime;
                    
                    //std::cout << blue << "\n\tfos = " << fos[0] << " = " << fos[1] << " + " << fos[2] << " + " << fos[3] << " >> t(s) = " << timeStar << normal;
                    //getchar();
                    totalOfImprovements++;
                }
                
            }else if (fol[0] < (1 + alpha)*fo[0] && vrpFeas && pppFeas){
                //atualizando a nova melhor solução
                r = cpRout(rl);
                p = cpProd(pl);
                vrp = std::vector<vrs> (vrpl);
                fo = std::vector<double> (fol);
                modSetup = true;
            }            
            else              
                modSetup = true;
                                          
            if (modSetup){
                if (setupsFixed) {
                    fixSetups(d,false,pl,ppp);
                    setupsFixed = false;
                }
            }else 
                setupsFixed = fixSetups(d,true,pl,ppp);
                            
            iterILS++;
        }
                        
        tFinal = steady_clock::now();
        diffTime = duration_cast<duration<double> > (tFinal-tStarIni);
    }while(diffTime <= vinte && iter <= maxIter);    
        
    //saida = instancia >> fo >> tempo até o otimo >> tempo total de execucao >> tempo total no PPP >> tempo total no VRP
    std::cout << "\n" << d.name << " " << fos[1] << " " << fos[2] << " " << fos[3] << " >> " << timeStar << " " << diffTime.count() << " " << totPPPtime << " " << totVRPtime << " >>";
    //quantidade de utilizacoes de cada solucao inicial para o roteamento
    for (int ls = 0; ls <= 4; ls++) std::cout << " " << sini[ls];
    std::cout << " >>";
    //quantidade de melhorias encontrada por cada busca inter-rota
    for (int ls = 0; ls <= 8; ls++) std::cout << " " << iers[ls];
    std::cout << " >>";
    //quantidade de melhorias encontrada por cada busca intra-rota
    for (int ls = 0; ls <= 7; ls++) std::cout << " " << iars[ls];
    std::cout << " >> " << totalOfImprovements;
    std::cout << " >> " << totalOfPerturbations;
    
    solprp star;
    star.p = cpProd(ps);
    star.r = cpRout(rs);
    star.vrp = std::vector<vrs> (vrps);
    star.fo = std::vector<double> (fos);
    star.ppp = ppp;    
    star.bestTime = timeStar;
    export2Sol(d,star);
    
    return double(fos[0]);
}

int main(int argc, char *argv[]){        
    data d;
    d.name = argv[1];    

    d.pathIn = "/home/varatan/Documentos/Dropbox/02-ufmg/02-doutorado/01-pesquisa/02-tese/03-implementacoes/00-instances/04-fortes/01-gerador-de-spmp/RPRP/";    
    d.pathOut = "/home/varatan/Documentos/Dropbox/02-ufmg/02-doutorado/01-pesquisa/02-tese/03-implementacoes/02-prp/04-prp-v-indx-h/00-outputs/02-td/01/";
    
    read(d);                        
    double fo = tils(d);
    //std::cout << fo;

    return 0;    
}
