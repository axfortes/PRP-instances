#include <iostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <random>
#include <algorithm>
#include "boost/dynamic_bitset.hpp"

using namespace std;

typedef struct
{
    char type; //tipo, (u)niforme
    int n; //numero de clientes
    int n1; //numero total de nos, considerando o deposito 0 e a copia n+1
    int m; //numero de linhas produtivas
    int P;  //numero de produtos
    int T;  //numero de periodos;

    //informacoes sobre a planta produtiva
    std::vector<int>  C;                //capacidade de produção de cada linha de cada item
    std::vector<double> l;              //custo de setup/ativacao de cada item em cada linha
    std::vector<double> u;              //custo unitário de produção de cada item em cada lkinha;


    //informacoes dos clientes
    std::vector<std::vector<std::vector<int> > > d;   //demanda periodica de cada item por cada cliente
    std::vector<std::vector<int> > a;                 //tempo de viagem de i para j
    std::vector<int> s;                              //tempo de servico de i
    std::vector<std::vector<double> > B;              //custo de atrasar o produto p para o cliente i

    //informacoes da frota
    int V;                              //numero de veiculos
    std::vector<int> Q;                 //capacidade maxima de transporte de determinado veiculo
    std::vector<int> e;                 //custo de setup/ativação de cada veiculo
    int H;                              //horizonte de tempo para entrega

    //d-00214510%
    //informacoes gerais
    std::vector<std::vector<double> > c; //custo de transporte de i para j
    std::vector<double> x;               //coordenada x do no i
    std::vector<double> y;               //coordenada y do no i
    std::vector<std::vector<double> > h; //custo de estocagem do produto p no no i
    std::vector<std::vector<int> > U;   //limite superior de estocagem de cada item na planta produtiva
    std::vector<std::vector<int> > I0;   //inventario inicial em cada no, deve satisfazer a demanda do periodo

    std::vector<std::vector<int> > N;           //closest neighbors
    std::vector<boost::dynamic_bitset<> >Nb;    //binary set of closest neighbors
    std::vector<double> Delta;                  //estimated visitation cost

} DATA;

void gerador(DATA &d);
void imprime(DATA d);
void vSort(std::vector<int> &V);
void help();
void createNeighborhood(std::vector<std::vector<double> > c, int n, std::vector<std::vector<int > > &N, int nSize, std::vector<boost::dynamic_bitset<> >&Nb);
void selection_sort(std::vector<double> &cost,  std::vector<int> &Ni, int tam);

int main()
{
    std::vector<int> nN;
    std::vector<int> nP;
    std::vector<int> nT;
    std::vector<int> nV;

    //primeiras instancias
    /*
    nN.push_back(5);
    nN.push_back(10);
    nN.push_back(15);
    //*/
    /*
    //segundas instancias
    nN.push_back(25);
    nN.push_back(50);
    nN.push_back(75);
    nN.push_back(100);

    nP.push_back(3);
    nP.push_back(4);
    nP.push_back(5);

    nT.push_back(3);
    nT.push_back(7);
    nT.push_back(10);
    nT.push_back(14);
    nT.push_back(21); //para segunda leva de instancias

    /*
    //primeiras intancias
    nV.push_back(3);
    nV.push_back(5);
    nV.push_back(7);
    //*/
    /*
    //intancias maiores para rotear
    nV.push_back(10);
    nV.push_back(11);
    nV.push_back(12);
    nV.push_back(13);
    nV.push_back(14);
    nV.push_back(15);

    DATA d;
    int inst = 0;
    //for (int nn = 0; nn < 3; nn++){
    for (int nn = 0; nn <= 3; nn++){
        d.n = nN[nn];
        for (int pp = 0; pp < 3; pp++){
            d.P = nP[pp];
            for (int tt = 0; tt < 4; tt++){
                d.T = nT[tt];
                //for (int vv = 0; vv < 3; vv++){
                for (int vv = 0; vv <= 5; vv++){
                    d.V = nV[vv];
                    gerador(d);
                    inst++;
                }
            }
        }
    }
    std::cout << "\n\tTotal de instancias geradas foi " << inst;
    //*/

    //criando instancia minima para validar formulações do AMPL e C++
    DATA d;
    d.n = 4;
    d.V = 2;
    d.T = 2;
    d.P = 1;
    gerador(d);
    //apagar apos validacao dos codigos
    //*/
    return 0;
}

void gerador(DATA &d){
    // Esse método é do C++11 pra poder gerar um número aleatório considerando uma distribuição UNIFORME
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937_64 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    //gerando as coordenadas dos pontos
    //std::uniform_real_distribution<> pos(0,100); //primeira leva de intancias
    std::uniform_real_distribution<> pos(0,200);
    //std::normal_distribution<double> pos(0.0,100.0);
    d.x = std::vector<double> (d.n+2,0.0);  //coordenada x do no i
    d.y = std::vector<double> (d.n+2,0.0);
    //reservar as posições 0 e n+1 para os nós artificiais
    for (int i = 0; i <= d.n; i++){
        d.x[i] = pos(gen);
        d.y[i] = pos(gen);
        if (i == 0){
            d.x[d.n+1] = d.x[0];
            d.y[d.n+1] = d.y[0];
        }
        //std::cout << "\ni: " << i << " -> x: " << d.x[i] << " -> y: " << d.y[i];
    }

    double dMed = 0.0;
    int cnt = 0;
    int maxC = 0;
    int minC = 999999;
    d.c = std::vector<std::vector<double> > (d.n+2, std::vector<double> (d.n+2, 0.0));
    //std::cout << "\n";
    for (int i = 0; i <= d.n+1; i++){
        for (int j = 0; j <= d.n+1; j++){
            d.c[i][j] = round(sqrt((d.x[i] - d.x[j])*(d.x[i] - d.x[j]) + (d.y[i] - d.y[j])*(d.y[i] - d.y[j])));
            //if (i != j && i < j) std::cout << "c[" << i << "][" << j << "] = " << d.c[i][j] << " ";
            dMed += d.c[i][j];
            cnt++;
            if (d.c[i][j] > maxC) maxC = d.c[i][j];
            if (d.c[i][j] < maxC && d.c[i][j] > 0) minC = d.c[i][j];
        }
        //std::cout << "\n";
    }
    //std::cout << "\n" << "\nc_med: " << dMed/cnt;

    d.N = std::vector<std::vector<int > > (d.n);
    d.Nb = std::vector<boost::dynamic_bitset<> > (d.n, boost::dynamic_bitset<> (d.n,0));
    createNeighborhood(d.c, d.n, d.N, 2, d.Nb);
    //Neighborhoods created

    d.Delta = std::vector<double> (d.n + 1, 0);
    for (int i = 1; i <= d.n; i++){
        int j = d.N[i-1][0];
        int k = d.N[i-1][1];
        d.Delta[i] = std::min(2*d.c[0][i],0.5*(d.c[i][j] + d.c[k][i]));
        //std::cout << "\nDelta[" << i << "] = " << d.Delta[i];
    }

    std::uniform_int_distribution<> dem(5,50);
    d.d = std::vector<std::vector<std::vector<int> > > (d.T + 1,std::vector<std::vector <int> > (d.P + 1, std::vector<int>(d.n + 1, 0)));
    int dt = 0;
    std::vector<std::vector<int> > DemTot (d.T + 1, std::vector<int> (d.P + 1,0));
    std::vector<int> dTot (d.P,0);
    std::vector<std::vector<int> > demMedia (d.n + 1, std::vector<int> (d.P + 1, 0));
    for (int t = 1; t <= d.T; t++){
        //std::cout << "\n>t: " << t;
        for (int p = 1; p <= d.P; p++){
            //std::cout << "\n>>p: " << p;
            for (int i = 1; i <= d.n; i++){
                d.d[t][p][i] = dem(gen);
                //std::cout << "\n>>>i: " << i << " -> d: " << d.d[t][p][i];
                DemTot[t][p] += d.d[t][p][i];
                demMedia[i][p] += d.d[t][p][i];
                dt += d.d[t][p][i];
            }
            //std::cout << "\nDemTot - t: " << t << " - p: " << p << " -> " << DemTot[t][p];
            dTot[p] += DemTot[t][p];
        }
    }

    d.U = std::vector<std::vector<int> > (d.P + 1, std::vector<int> (d.n + 1, 0));
    for (int i = 0; i <= d.n; i++){
        //std::cout << "\>ni: " << i;
        for (int p = 1; p <= d.P; p++){
            if (i == 0){
                std::uniform_int_distribution<> holdCapacity(0.05*dTot[p],0.3*dTot[p]);
                d.U[p][i] = holdCapacity(gen); //posso estocar toda a producao periodica mais 50%
            }else if (i > 0){
                demMedia[i][p] = std::ceil(demMedia[i][p]/d.T);
                d.U[p][i] = std::ceil(demMedia[i][p]*2);
            }
            //std::cout << "\n>>p: " << p << " - > U: " << d.U[p][i];
        }
    }

    int demandTotal = 0;
    //as capacidades de produção por produto
    d.C = std::vector<int> (d.P + 1, 0);
    d.l = std::vector<double> (d.P + 1,0); //custo de preparação de produção (setup)
    d.u = std::vector<double> (d.P + 1,0.0); //custo de producao unitario
    std::uniform_int_distribution<> prodCost(2,8);
    std::uniform_int_distribution<> setupCost(1000,5000);
    //std::cout << "\ndemanda total: ";
    for(int p = 1; p <= d.P; p++) {
        std::uniform_int_distribution<> prodCapacity(round(0.5*dTot[p]),round(0.8*dTot[p]));
        demandTotal += dTot[p];
        d.C[p] = prodCapacity(gen);
        d.l[p] =  setupCost(gen);
        d.u[p] =  prodCost(gen);
        //std::cout << "\np: " << p << " -> C = " << d.C[p] << " -> l = " << d.l[p] << " -> u = " << d.u[p];
    }

    //definindo o I0 igual a demanda do primeiro período para cada cliente
    d.I0 = std::vector<std::vector<int> > (d.P + 1, std::vector<int> (d.n + 1,0));
    for (int p = 1; p <= d.P; p++){
        for (int i = 1; i <= d.n; i++){
            d.I0[p][i] = d.d[1][p][i];
            //std::cout << "\nI0[" << p << "][" << i << "] = " << d.I0[p][i];
        }
    }

    int medQ = ceil(demandTotal/d.T); //carga média por periodo
    //std::cout << "\nV:" << d.V <<"\nmedQ: " << medQ;
    //getchar();
    d.Q = std::vector<int> (d.V + 1, 0);
    d.e = std::vector<int> (d.V + 1, 0);
    std::uniform_int_distribution<> maxVehicleLoad(round(0.8*medQ/d.V),round(medQ/d.V));
    std::uniform_int_distribution<> setupVehicle(500,2000);
    for (int v = 1; v <= d.V; v++){
        d.Q[v] = maxVehicleLoad(gen);
        d.e[v] = setupVehicle(gen);
        //std::cout << "\nv: " << v << " - Q = " << d.Q[v] << " - e = " << d.e[v];
    }
    //ordenando os custos e capacidades
    vSort(d.Q);
    vSort(d.e);

    //custos de estocagem
    std::uniform_int_distribution<> holding(3,5);
    d.h = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n+1,0.0));
    double maxHoldCost = 0;
    for (int p = 1; p <= d.P; p++){
        //std::cout << "\np: " << p;
        for (int i = 0; i <= d.n; i++){
            d.h[p][i] = holding(gen);
            //std::cout << "\ni: " << i << " -> h = " << d.h[p][i];
            //if (d.h[p][i] > maxHoldCost) maxHoldCost = d.h[p][i];
        }
    }

    //https://www.agatetepe.com.br/estimando-custo-de-atraso/, pode ser até 10 vezes os custo real
    std::uniform_int_distribution<> backorder(5,10);
    d.B = std::vector<std::vector<double> > (d.P + 1, std::vector<double> (d.n+1,0.0));
    for (int p = 1; p <= d.P; p++){
        std::cout << "\np: " << p;
        for (int i = 1; i <= d.n; i++){
            //d.B[p][i] = 100*backorder(gen);
            //d.B[p][i] = d.h[p][i]*backorder(gen);
            d.B[p][i] = (d.h[p][i] + d.u[p] + d.Delta[i])*backorder(gen);
            //std::cout << "\ni: " << i << " -> B = " << d.B[p][i];
        }
    }


    //d.H = 7*(minC + 0.5*(maxC - minC));
    d.H = 10*60;
    //std::cout << "\nmaxC: " << maxC << " " << " - minC: " << minC << " - H: " << d.H;

    //d.st = std::vector<int> (d.nC + 1, floor(maxC/15));
    d.s = std::vector<int> (d.n + 1, 10);

    imprime(d);
}

void imprime(DATA d)
{
    ofstream arq;

    string n = std::to_string(d.n);
    //string m = std::to_string(d.m);
    string T = std::to_string(d.T);
    string P = std::to_string(d.P);
    string V = std::to_string(d.V);
    string name = "X" + n + "C" + P + "P" + T + "T" + V + "V" + ".dat";
    string outputfile = "/home/varatan/Documentos/Dropbox/02-ufmg/02-doutorado/01-pesquisa/02-tese/03-implementacoes/00-instances/04-fortes/01-gerador-de-spmp/x2/" + name;
    //string outputfile = "/home/varatan/Documentos/Dropbox/02-ufmg/02-doutorado/01-pesquisa/02-tese/03-implementacoes/00-instances/04-fortes/01-gerador-de-spmp/x3/" + name;

    //convertendo string para char
    char *cName = new char[outputfile.length()+1];
    memcpy(cName, outputfile.c_str(), outputfile.length() + 1);
    arq.open(cName);

    if (!arq.is_open()) help();

    arq << d.n;
    //arq << "\n" << d.m;
    arq << "\n" << d.T;
    arq << "\n" << d.P;
    arq << "\n" << d.V;
    arq << "\n" << d.H;

    for (int v = 1; v <= d.V; v++)
        arq << "\n" << d.Q[v] << "\t" << d.e[v];

    for (int i = 0; i <= d.n; i++){
        arq << "\n" << fixed << setprecision(1) << d.x[i] << "\t" << d.y[i];
        if (i > 0){
            arq << "\t" << d.s[i];
            for (int p = 1; p <= d.P; p++)
                arq << "\t" << d.I0[p][i] << "\t" << d.h[p][i] << "\t" << d.U[p][i] << "\t" << d.B[p][i];
        }else if (i == 0){
            for (int p = 1; p <= d.P; p++)
                arq << "\t" << d.I0[p][i] << "\t" << d.h[p][i] << "\t" << d.U[p][i] << "\t" << d.C[p] << "\t" << d.l[p]<< "\t" << d.u[p];
        }
    }

    for (int p = 1; p <= d.P; p++){
        for (int i = 1; i <= d.n; i++){
            arq << "\n";
            for (int t = 1; t <= d.T; t++){
                arq << d.d[t][p][i] << "\t";
            }
        }
    }
    arq.close();
}

void vSort(std::vector<int> &V) {
	for (int fixo = 0; fixo < V.size() - 1; fixo++) {
		int menor = fixo;

		for (int i = menor + 1; i < V.size(); i++)
			if (V[i] < V[menor]) menor = i;

		if (menor != fixo) {
			int sV = V[fixo];
			V[fixo] = V[menor];
			V[menor] = sV;
		}
	}
}

void help()
{
    std::cout << std::endl << std::endl << "exec [data file] \n " << std::endl;
    exit(1);
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
