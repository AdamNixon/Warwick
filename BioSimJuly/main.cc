/**************************************************************************
 
  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library 
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2015 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger 
  Copyright (C) 2004 - 2015 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2015 Mirko Kraenkel
  Copyright (C) 2006 - 2015 Christoph Gersbacher
  Copyright (C) 2006 - 2015 Martin Nolte
  Copyright (C) 2011 - 2015 Tobias Malkmus
  Copyright (C) 2012 - 2015 Stefan Girke
  Copyright (C) 2013 - 2015 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner
  Copyright (C) 2015        Marco Agnese
  Copyright (C) 2015        Martin Alkaemper


  The dune-fem module is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
**************************************************************************/
#include <config.h>

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#warning "Experimental grid extensions required for GeoGridPart. Reconfigure with --enable-experimental-grid-extensions to enable GeoGridPart."
int main() { return 1; }
#else

// iostream includes
#include <iostream>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

#warning DEFORMATION
// include geometrty grid part
#include <dune/fem/gridpart/geogridpart.hh>
// include description of surface deformation
#include "deformation.hh"

// include header for heat model
#include "heat.hh"

#include "heatmodel.hh"
#include "heatscheme.hh"

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
double algorithm ( HGridType &grid, int step )
{
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 6 > FunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 3 > FunctionSpaceTypeG;
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 1 > SurFunctionSpaceType;
  // create time provider
  Dune::Fem::GridTimeProvider< HGridType > timeProvider( grid );

  // we want to solve the problem on the leaf elements of the grid
  //! [Setup the grid part for a deforming domain]
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > HostGridPartType;
  HostGridPartType hostGridPart( grid );
  
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceTypeG, HostGridPartType, POLORDER > DeformationSpaceType;
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DeformationSpaceType > DeformationDiscreteFunctionType;
   // assuming dune-istl used in femscheme
  DeformationSpaceType deformationSpace(hostGridPart);
 
  typedef DeformationCoordFunction< HGridType::dimensionworld > DeformationType;
  DeformationType analyticDeformation;
  typedef Dune::Fem::GridFunctionAdapter< DeformationType, HostGridPartType > DiscreteDeformationType;
  DiscreteDeformationType discreteDeformation( "deformation", analyticDeformation, hostGridPart, 3 );

  DeformationDiscreteFunctionType deformation("deformation",deformationSpace);
  interpolate( discreteDeformation, deformation );
  DeformationDiscreteFunctionType deformationG("deformationG",deformationSpace);
  interpolate( discreteDeformation, deformationG );

  typedef Dune::Fem::GeoGridPart< DeformationDiscreteFunctionType > GridPartType;
  GridPartType gridPart( deformation );
  GridPartType gridPartG( deformationG );
  //! [Setup the grid part for a deforming domain]

// Velocity test
 typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceTypeG, GridPartType, POLORDER > DiscreteFunctionSpaceTypeG;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceTypeG > DiscreteFunctionTypeG;
DiscreteFunctionSpaceTypeG cortexSpace(gridPart);
DiscreteFunctionTypeG cortex("ball", cortexSpace );
DiscreteFunctionSpaceTypeG IntialSpace(gridPart);
DiscreteFunctionTypeG cortex1("ball", IntialSpace );
DiscreteFunctionSpaceTypeG velocitySpace(gridPart);
DiscreteFunctionTypeG velocity("ball", cortexSpace );

     for (int i=0; i<deformation.blocks(); ++i)
    {
	*cortex.block(i)=*deformation.block(i);
	*cortex1.block(i)=*deformation.block(i);
	*velocity.block(i)=*deformation.block(i);
    }

 typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SurFunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceTypeS;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceTypeS > DiscreteFunctionTypeS;
DiscreteFunctionSpaceTypeS phaseSpace(gridPart);
DiscreteFunctionTypeS phasefield("phase field", phaseSpace );

   for (int i=0; i<phasefield.blocks(); ++i)
    {
	*phasefield.block(i)=0; 
    }

  // type of the mathematical model used
  typedef TimeDependentCosinusProduct< FunctionSpaceType > GeoProblemType;
  typedef HeatModelG< FunctionSpaceType, FunctionSpaceTypeG, GridPartType, DiscreteFunctionTypeG, DiscreteFunctionTypeS  > GeoModelType;
  typedef TimeDependentSurfaceProblem< SurFunctionSpaceType > SurProblemType;
  typedef HeatModelS< SurFunctionSpaceType, GridPartType,  DiscreteFunctionTypeG, DiscreteFunctionTypeS > SurModelType; 

  GeoProblemType problemG( timeProvider );
  SurProblemType problemS( timeProvider );

  double area=0., volume=0.;
  // implicit model for left hand side
  GeoModelType implicitModelG( problemG, gridPart, cortex, cortex1, phasefield, true, false, area, volume );
  GeoModelType setupimplicitModelG( problemG, gridPart, cortex, cortex1, phasefield, true, true, area, volume );
  // explicit model for right hand side
  GeoModelType explicitModelG( problemG, gridPart, cortex, cortex1, phasefield, false, false, area, volume );
  GeoModelType setupexplicitModelG( problemG, gridPart, cortex, cortex1,  phasefield, false, true, area, volume );

 SurModelType implicitModelS( problemS, gridPart, velocity, phasefield, true );
 SurModelType explicitModelS( problemS, gridPart, velocity, phasefield, false );

  // create heat scheme
  typedef HeatScheme< GeoModelType, GeoModelType, SurModelType, SurModelType > SchemeType;
  SchemeType scheme( gridPart, gridPartG, implicitModelG, explicitModelG, setupimplicitModelG, setupexplicitModelG, implicitModelS, explicitModelS, area, volume );

  typedef Dune::Fem::GridFunctionAdapter< SurProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problemS, gridPartG, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionTypeS *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solutionSS()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );
  

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "heat.endtime", 2.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("heat.reducetimestepfactor", 1 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "heat.timestep", 0.125 );

  timeStep *= pow(dtreducefactor,step);
 
  //! [time loop]
  // initialize with fixed time step
  timeProvider.init( timeStep ) ;
 
  // initialize scheme and output initial data
  scheme.initialize();

  scheme.Setup_prepare();
  scheme.Setup_solve();

  std::cout << "setup complete" << std::endl;
    auto SIt = scheme.solutionS().dbegin();
    auto SEnd = scheme.solutionS().dend();
    auto PIt = phasefield.dbegin();
    for ( ; SIt != SEnd; ++SIt, ++PIt)
    {
    *PIt = *SIt;
//std::cout << "solution_i: " << *PIt << std::endl;
    } 

  // time loop, increment with fixed time step
  for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
  //! [time loop]
  {
    dataOutput.write( timeProvider );
    scheme.prepare();
    //! [Set the new time to move to new surface]
    //deformation.setTime( timeProvider.time() + timeProvider.deltaT() ); 
   
     // solve once - but now we need to reassmble
    scheme.solve(true);
 
    //! [Set the new time to move to new surface]
    double count=0;
    double radius=0;
    auto deformIt = deformationG.dbegin();
    auto deformEnd = deformationG.dend();
    auto velIt = velocity.dbegin();
    auto coxIt = cortex.dbegin();
    auto solIt = scheme.solutionGrid().dbegin();
    for ( ; deformIt != deformEnd; ++deformIt, ++solIt, ++velIt )
    {
      *deformIt = *solIt; 
      *velIt=*solIt; 
     // *coxIt=*solIt;
      count+=1;
      radius += *solIt**solIt;
     // std::cout << "solution_i: " << *solIt << std::endl;
    }
   
  scheme.solveS(true); 
count=0.;
radius=0.;
    auto SIt = scheme.solutionS().dbegin();
    auto SEnd = scheme.solutionS().dend();
    auto PIt = phasefield.dbegin();
    for ( ; SIt != SEnd; ++SIt, ++PIt)
    {
    *PIt = *SIt;
count+=1;
radius+=*SIt;
    } 
std::cout << "phase: " << radius/count << std::endl;

  }
 

  // output final solution
  dataOutput.write( timeProvider );
 
  // select norm for error computation
  typedef Dune::Fem::L2Norm< GridPartType > NormType;
  NormType norm( gridPartG ); 
  return norm.distance( gridExactSolution, scheme.solutionS() );

}


// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "../data/parameter" );

  // type of hierarchical grid
  typedef Dune::GridSelector::GridType  HGridType ;

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "heat.repeats", 0 );

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "heat.level" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

	// refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // calculate first step
  double oldError = algorithm( grid, (repeats > 0) ? 0 : -1 );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    const double newError = algorithm( grid, step );
    const double eoc = log( oldError / newError ) / M_LN2;
    if( Dune::Fem::MPIManager::rank() == 0 )
    {
      std::cout << "Error: " << newError << std::endl;
      std::cout << "EOC( " << step << " ) = " << eoc << std::endl;
    }
    oldError = newError;
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
#endif
