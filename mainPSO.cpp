
#include <iostream>

#include <OptFrame/Heuristics/EvolutionaryAlgorithms/PSO.hpp>
#include <algorithm>
#include <functional>
#include <iostream>

#include <OptFCore/FCore.hpp>
#include <OptFrame/Core.hpp>
#include <OptFrame/Heuristics/Heuristics.hpp> // many metaheuristics here...
#include <OptFrame/Scanner++/Scanner.hpp>
#include <OptFrame/Util/Matrix.hpp>
#include <OptFrame/Util/printable.h>

using namespace std;
using namespace optframe;

// next definitions come here within namespace
// this also works when defining in global scope (same as 'class')
namespace example_pso {

class ProblemContext // 100% copiado do brkga
{
public:
   int n; //number of days
   vector<double> dadosreais;

// load data from Scanner
void load(Scanner& scanner)
{
   n = *scanner.nextInt();      // reads number of days
   int acu = 0;
   //
   vector<int> dadosreais(n);
   //
   for (int i = 0; i < n; i++) 
   {
      scanner.next();
      acu += *scanner.nextInt(); // does cumulated sum of cases
      dadosreais[i] = *scanner.nextDouble(); 
   }
} //load  
}; //class

ProblemContext pPSO; // onde eu guardo as variaveis todas

using ESolutionPSO = std::pair<
  std::vector<double>, // first part of search space element: solution (representation)
  Evaluation<double>   // second part of search space element: evaluation (objective value)
  >;

double poptot = 513584;
double v = 1/7;
double eta = 1/10;

Evaluation<double>
fevaluatePSO(const std::vector<double>& s)
{
   double fit = 0;
   double S, I, R, U;
   
   // to be estimated
   I = s[0];
   U = s[1];
   double tau = s[2];
   double f = s[3];

   S = poptot - I - U;

   for (unsigned i = 0; i < pPSO.n; i++)
   {
      S += - (tau / poptot) * S * (I + U);
      I += (tau / poptot) * S * (I + U) - v * I;
      R += f * v * I - eta * R;
      U += (1 - f) * v * I - eta * U;

      fit += pow(R - pPSO.dadosreais[i], 2);
   }
   return Evaluation<double>{ fit/pPSO.n };
}
  
// Evaluate
FEvaluator<ESolutionPSO, MinOrMax::MINIMIZE> evPSO{
   fevaluatePSO
};

//
} // example_pso

// import everything on main()
using namespace std;
using namespace optframe;
using namespace scannerpp;
using namespace example_pso;

RandGen myRandGen; // good random generator

int
main()
{
   Scanner scanner{ File{ "berlin52.tsp" } };
   pPSO.load(scanner);

   int seed = 1000012345;


   // Particle Swarm Optimization
   // PSO(Evaluator<S, XEv, XES>& evaluator, unsigned pop_size, unsigned iter_max, vector<double> cI, vector<double> cS, RandGen& _rg)

   myRandGen.setSeed(seed);
   int nParam = 1; // param size

   //quotas for I0, U0, tau, f
   vector<double> cI = {1, 1, 0, 0};
   vector<double> cS = {100, 100, 0.1, 1.0};

   PSO<> myPSO{
      evPSO,
      10, // pop_size
      50, // iter_max
      cI, // number of parameters is extracted from cI.size()
      cS,
      myRandGen
   };
   myPSO.setVerboseR();

   auto statusPSO = myPSO.search(5.0); // 5.0 seconds max

   std::cout << "BEST = " << myPSO.best->second.evaluation()
             << " sol: " << myPSO.best->first << std::endl;

   return 0;
}
