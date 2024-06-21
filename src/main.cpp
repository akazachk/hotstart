#include <cstdio>
#include <iostream> // cerr
#include <string>
#include <vector>
#include <limits> // numeric_limits<int>::max()

#ifdef USE_CLP_SOLVER
  #include <OsiClpSolverInterface.hpp>
  using SolverInterface = OsiClpSolverInterface; ///< SolverInterface is set to OsiClpSolverInterface
#else
  #include <OsiSolverInterface.hpp>
  using SolverInterface = OsiSolverInterface; ///< SolverInterface is set to OsiSolverInterface (this will not work)
#endif

/// @brief Overload solve from hot start because of issues
bool solveFromHotStart(OsiSolverInterface* const solver, const int col,
    const bool isChangedUB, const double origBound, const double newBound,
    int& numhotStartViolations, const int MAX_NUM_HOT_START_VIOLS = 2);

/// @brief Overload #solveFromHotStart for when hot start is not disabled
inline bool solveFromHotStart(OsiSolverInterface* const solver, const int col,
    const bool isChangedUB, const double origBound, const double newBound,
    const int MAX_NUM_HOT_START_VIOLS = 2) {
  int numHotStartViolations = 0;
  return solveFromHotStart(solver, col, isChangedUB, origBound, newBound, numHotStartViolations, MAX_NUM_HOT_START_VIOLS);
} // solveFromHotStart

/// @brief Generate SplitDisjunction for all variables that are fractional, using hot starts
int generateSplitDisjunctions(
    std::vector<int>& fracCoreSelected,
    const OsiSolverInterface* const si,
    const double EPS = 1e-7,
    const double AWAY = 1e-3,
    const int MAX_NUM_HOT_START_VIOLS = 2);

/// @brief Use hot starting to check split disjunctions
int main(int argc, char** argv) {
  if (argc == 1) {
    fprintf(stderr, "*** ERROR: Need to specify filename as first argument.\n");
    return 1;
  }
  const int MAX_NUM_HOT_START_VIOLS = (argc >= 3) ? std::stoi(argv[2]) : 2;
  const double EPS = 1e-7;

  const std::string filename = argv[1];
  printf("Instance: %s\n", filename.c_str());

  OsiSolverInterface* solver = NULL;
  solver = new SolverInterface;
  solver->readMps(filename.c_str());
  solver->initialSolve();
  if (!solver->isProvenOptimal()) {
    fprintf(stderr, "*** ERROR: Solver LP relaxation not optimal.\n");
  } else {
    // Solver is proven optimal; try generating splits
    std::vector<int> fracCore;
    generateSplitDisjunctions(fracCore, solver, EPS, 1e-3, MAX_NUM_HOT_START_VIOLS);
  }

  if (solver) { delete solver; }
  return 0;
} /* main */

/// @details Some problems with normal OsiSolverInterface::solverFromHotStart() can happen,
/// such as with arki001, in which we can have OsiSolverInterface::isIterationLimitReached().
/// The plan then is to try one more time with a hot start, and otherwise switch to resolving
/// (depending on the parameter \p MAX_NUM_HOT_START_VIOLS).
///
/// @return OsiSolverInterface::isProvenOptimal()
bool solveFromHotStart(
    /// [in,out] Original solver
    OsiSolverInterface* const solver,
    /// [in] Which variable is being modified
    const int col,
    /// [in] Is the upper bound changing?
    const bool isChangedUB,
    /// [in] Old value of the variable bound
    const double origBound,
    /// [in] New value of the variable bound
    const double newBound,
    /// [in,out] If in a prior call, hot start seemed to not work, we might disable it, in which case we need to use normal resolving
    int& numHotStartViolations,
    /// [in] Maximum number of \p numHotStartViolations before switching to resolve
    const int MAX_NUM_HOT_START_VIOLS) {
  if (MAX_NUM_HOT_START_VIOLS  < 0) {
    solver->solveFromHotStart();
    return (solver->isProvenOptimal());
  }

  if (numHotStartViolations < MAX_NUM_HOT_START_VIOLS) {
    solver->solveFromHotStart();

    if (solver->isIterationLimitReached()) {
      numHotStartViolations++;
      // This sometimes happens, e.g., with arki001
      solver->unmarkHotStart();
      solver->resolve();
      if (isChangedUB) {
        solver->setColUpper(col, origBound);
      } else {
        solver->setColLower(col, origBound);
      }
      solver->resolve();
      solver->markHotStart();
      if (isChangedUB) {
        solver->setColUpper(col, newBound);
      } else {
        solver->setColLower(col, newBound);
      }
      solver->solveFromHotStart();
      if (solver->isIterationLimitReached()) {
        numHotStartViolations++;
      }
    }
  } // check if hot start is not disabled

  if (numHotStartViolations >= MAX_NUM_HOT_START_VIOLS) {
    // Else, hot start is disabled, revert to resolve
    solver->resolve();
  }
  return (solver->isProvenOptimal());
} /* solveFromHotStart overload */

/// @return number of split disjunctions generated
int generateSplitDisjunctions(
    std::vector<int>& fracCoreSelected,
    const OsiSolverInterface* const si,
    const double EPS,
    const double AWAY,
    const int MAX_NUM_HOT_START_VIOLS) {
  std::vector<int> fracCore = si->getFractionalIndices(AWAY);
  if (fracCore.size() == 0) {
    return 0;
  }

  int num_splits = 0;
  fracCoreSelected.clear();
  fracCoreSelected.reserve(fracCore.size());

  // Set up solver for hot start
  SolverInterface* solver;
  try {
    solver = dynamic_cast<SolverInterface*>(si->clone());
  } catch (std::exception& e) {
    fprintf(stderr,
        "Unable to clone solver into desired SolverInterface.\n");
    exit(1);
  }
#ifdef USE_CLP
  try {
    const int hot_start_iter_limit = std::numeric_limits<int>::max();
    solver->setIntParam(OsiMaxNumIterationHotStart, hot_start_iter_limit);
    solver->setSpecialOptions(16); // use standard strong branching rather than clp's
  } catch (std::exception& e) {
    // It's okay, we can continue
  }
#endif

  std::vector<double> sortCriterion;
  sortCriterion.reserve(fracCore.size());
  std::vector<double> fracCoreVal(solver->getColSolution(), solver->getColSolution() + solver->getNumCols());

  // START OF MAIN HOT START CODE
  solver->enableFactorization();
  solver->markHotStart();

  // Loop over each fractional variable and check each side of the split disjunction
  for (int var : fracCore) {
    const double val = fracCoreVal[var];
    const double floorxk = std::floor(val);
    const double ceilxk = std::ceil(val);

    if (!si->isInteger(var)) {
      fprintf(stderr,
          "Chosen variable %d is not an integer variable.\n", var);
      exit(1);
    }

    const bool isFloor = (std::abs((val) - (floorxk)) <= EPS);
    const bool isCeil = (std::abs((val) - (ceilxk)) <= EPS);
    if (isFloor || isCeil) {
      fprintf(stderr,
          "Chosen variable %d is not fractional (value: %1.6e).\n", var, val);
      exit(1);
    }

    const double origLB = solver->getColLower()[var];
    const double origUB = solver->getColUpper()[var];
    bool downBranchFeasible = true, upBranchFeasible = true;
    double downBound = std::numeric_limits<double>::max();
    double upBound = std::numeric_limits<double>::max();

    // Check down branch
    solver->setColUpper(var, floorxk);
    solveFromHotStart(solver, var, true, origUB, floorxk, MAX_NUM_HOT_START_VIOLS);
    if (solver->isProvenOptimal()) {
      downBound = solver->getObjValue();
    } else if (solver->isProvenPrimalInfeasible()) {
      downBranchFeasible = false;
    } else {
      // Something strange happened
      fprintf(stderr,
          "Down branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      exit(1);
    }
    solver->setColUpper(var, origUB);

    // Return to previous state
    // solver->solveFromHotStart();
    solveFromHotStart(solver, var, true, origUB, origUB, MAX_NUM_HOT_START_VIOLS);

    // Check up branch
    solver->setColLower(var, ceilxk);
    solveFromHotStart(solver, var, false, origLB, ceilxk, MAX_NUM_HOT_START_VIOLS);
    if (solver->isProvenOptimal()) {
      upBound = solver->getObjValue();
    } else if (solver->isProvenPrimalInfeasible()) {
      upBranchFeasible = false;
    } else {
      // Something strange happened
      fprintf(stderr,
          "Up branch is neither optimal nor primal infeasible on variable %d (value %e).\n",
          var, val);
      exit(1);
    }
    solver->setColLower(var, origLB);

    // Return to original state
    // solver->solveFromHotStart();
    solveFromHotStart(solver, var, false, origLB, origLB, MAX_NUM_HOT_START_VIOLS);

    // Check if some side of the split is infeasible
    if (!downBranchFeasible || !upBranchFeasible) {
      if (!downBranchFeasible && !upBranchFeasible) {
        // Infeasible problem
        fprintf(stderr,
            "Infeasible problem due to integer variable %d (value %e).\n",
            var, val);
        exit(1);
      }
    }
    else {
      num_splits++;
      fracCoreSelected.push_back(var);
      sortCriterion.push_back(CoinMin(downBound, upBound));
    }
  } // loop through fractional core
  solver->unmarkHotStart();
  solver->disableFactorization();
  if (solver) {
    delete solver;
  }

  // Sort by decreasing strong branching lb
  std::vector<unsigned> sortIndex(fracCoreSelected.size());
  for (unsigned i = 0; i < sortIndex.size(); i++) {
    sortIndex[i] = i;
  }
  std::sort(sortIndex.begin(), sortIndex.end(),
      [&](const unsigned i, const unsigned j)
      { return sortCriterion[i] > sortCriterion[j]; } );

  #ifdef TRACE
  { // DEBUG
    for (int i = 0; i < num_splits; ++i) {
      printf("Var: %d", fracCore[i]);
      printf("\tSort criterion: %f\n", sortCriterion[sortIndex[i]]);
    }
  }
  #endif
  return num_splits;
} /* generateSplitDisjunctions */
