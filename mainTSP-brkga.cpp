//g++ -I../../src -std=c++17 -fconcepts mainTSP-brkga.cpp

#include <iostream>

#include <OptFrame/Heuristics/EvolutionaryAlgorithms/BRKGA.hpp>
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
using namespace scannerpp;

// next definitions come here within namespace
// this also works when defining in global scope (same as 'class')
namespace TSP_brkga {

// define TSP solution type as 'vector<int>', using 'double' as evaluation type
using ESolutionTSP = std::pair<
  std::vector<int>,  // first part of search space element: solution (representation)
  Evaluation<double> // second part of search space element: evaluation (objective value)
  >;

// TSP problem context and data reads
class ProblemContext
{
public:
   int n;               // number of clients
   Matrix<double> dist; // distance matrix (Euclidean)
   // load data from Scanner
   void load(Scanner& scanner)
   {
      n = *scanner.nextInt();      // reads number of clients
      dist = Matrix<double>(n, n); // initializes n x n matrix
      //
      vector<double> xvalues(n);
      vector<double> yvalues(n);
      //
      for (int i = 0; i < n; i++) {
         scanner.next();
         xvalues[i] = *scanner.nextDouble(); // reads x
         yvalues[i] = *scanner.nextDouble(); // reads y
      }
      // calculate distance values, for every client pair (i,j)
      for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++)
            dist(i, j) = distance(xvalues.at(i), yvalues.at(i), xvalues.at(j), yvalues.at(j));
   }
   // euclidean distance (double as return)
   double distance(double x1, double y1, double x2, double y2)
   {
      return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
   }
};
// Create TSP Problem Context
ProblemContext pTSP;

Evaluation<double>
fevaluate(const std::vector<int>& s)
{
   double f = 0;
   for (int i = 0; i < int(pTSP.n) - 1; i++)
      f += pTSP.dist(s[i], s[i + 1]);
   f += pTSP.dist(s[int(pTSP.n) - 1], s[0]);
   return Evaluation<double>{ f };
}

// Evaluate
FEvaluator<ESolutionTSP, MinOrMax::MINIMIZE>
  ev{
     fevaluate
  };

// ===========================
// decoder function for BRKGA (and alike)

pair<Evaluation<double>, vector<int>>
fDecode(const vector<double>& rk)
{
   vector<pair<double, int>> v(rk.size());
   int k = 0;
   for (unsigned i = 0; i < v.size(); i++)
      v[k] = pair<double, int>(rk[i], i);

   sort(v.begin(), v.end(), [](const pair<double, int>& i, const pair<double, int>& j) -> bool {
      return i.first < j.first;
   });

   // R = vector<int>
   vector<int> p(v.size());
   for (unsigned i = 0; i < v.size(); i++)
      p[i] = v[i].second;

   Evaluation<double> e = ev.evaluate(p);
   return make_pair(e, p);
}

// evaluator random keys (for TSP)
FDecoderRK<std::vector<int>, Evaluation<>, double, MinOrMax::MINIMIZE> decoder{
   fDecode
};

//
} // TSP_brkga

// import everything on main()
using namespace std;
using namespace optframe;
using namespace scannerpp;
using namespace TSP_brkga;

class MyRandomKeysInitPop : public InitialPopulation<std::vector<double>, Evaluation<double>>
{
   using RSK = std::vector<double>;

private:
   int sz;
   double prc;

public:
   MyRandomKeysInitPop(int size, double precision)
     : sz{ size }
     , prc{ precision }
   {
   }

   Population<RSK, Evaluation<double>> generatePopulation(unsigned populationSize, double timelimit) override
   {
      Population<RSK, Evaluation<double>> pop;

      for (unsigned i = 0; i < populationSize; i++) 
      {
         vector<double>* d = new vector<double>(sz);
         for (int j = 0; j < sz; j++)
            d->at(j) = ( rand() % static_cast<int>(prc) ) / static_cast<int>(prc);
         pop.push_back(d);
      }

      return pop;
   }
};

int
main()
{
   long seed = 12345678; // setting initial random seed

   // load data into problem context 'pTSP'
   Scanner scanner{ File {"berlin52.tsp"} };
   
   pTSP.load(scanner);
   std::cout << "file loaded" << std::endl;

   for (int i = 32; i > 1; i--){ // loop to test brkga seed boundaries

      ::srand(seed + i*i); // pseudo randomizes sytem random

      InitialPopulation<vector<double>, ESolutionTSP::second_type>* initPop =
      new MyRandomKeysInitPop(pTSP.n, pow(2, i)); // passing key_size

      // Parameters BRKGA
      // (C1): Evaluator<S, XEv>& _evaluator, int key_size, unsigned numGen, unsigned _popSize, double fracTOP, double fracBOT, double _probElitism) :

      //eprk, pTSP.n, 1000, 30, 0.4, 0.3, 0.6
      BRKGA<ESolutionTSP::first_type, ESolutionTSP::second_type, double, ESolutionTSP> brkga(
      decoder,
      *initPop,
      1000,
      30,
      0.4,
      0.3,
      0.6);

      //delete MyRandomKeysInitPop(pTSP.n, pow(2, i)); // killing sample size to avoid memory loss
      // unable to do the thing above because it is not a pointer
   }

   std::cout << "FINISHED" << std::endl;
   return 0;
}

//g++ -I../../src -std=c++17 -fconcepts mainTSP-brkga.cpp