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
#ifndef ELLIPT_FEMSCHEME_HH
#define ELLIPT_FEMSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/solver/cginverseoperator.hh>

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// local includes
#include "probleminterface.hh"

#include "model.hh"

#include "normal.hh"
#include "rhs.hh"
#include "elliptic.hh"

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

// FemScheme
//----------

/*******************************************************************************
 * template arguments are:
 * - GridPsrt: the part of the grid used to tesselate the
 *             computational domain
 * - Model: description of the data functions and methods required for the
 *          elliptic operator (massFlux, diffusionFlux)
 *     Model::ProblemType boundary data, exact solution,
 *                        and the type of the function space
 *******************************************************************************/
template < class ModelG, class ModelS >
class FemScheme
{
public:
  //! type of the mathematical model
  typedef ModelG ModelTypeG ;
  typedef ModelS ModelTypeS ;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelTypeG::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename ModelTypeG :: FunctionSpaceType   FunctionSpaceType;
  typedef typename ModelTypeG :: FunctionSpaceTypeG   FunctionSpaceTypeG;
  typedef typename ModelTypeS :: FunctionSpaceType   FunctionSpaceTypeS;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceTypeG, GridPartType, POLORDER > DiscreteFunctionSpaceTypeG;
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceTypeS, GridPartType, POLORDER > DiscreteFunctionSpaceTypeS;

  // choose type of discrete function, Matrix implementation and solver implementation
#if HAVE_DUNE_ISTL && WANT_ISTL
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
 // typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
//  typedef Dune::Fem::ISTLGMResOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
  typedef Dune::Fem::ISTLBICGSTABOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;

  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceTypeG > DiscreteFunctionTypeG;

  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceTypeS > DiscreteFunctionTypeS;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionTypeS, DiscreteFunctionTypeS > LinearOperatorTypeS;
  typedef Dune::Fem::ISTLCGOp< DiscreteFunctionTypeS, LinearOperatorTypeS > LinearInverseOperatorTypeS;
#else
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceTypeG > DiscreteFunctionTypeG;

  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceTypeS > DiscreteFunctionTypeS;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionTypeS, DiscreteFunctionTypeS > LinearOperatorTypeS;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionTypeS > LinearInverseOperatorTypeS;
#endif

  /*********************************************************/

  //! define Laplace operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelTypeG > EllipticOperatorType;
  typedef DifferentiableEllipticOperator< LinearOperatorTypeS, ModelTypeS > EllipticOperatorTypeS;

  FemScheme( GridPartType &gridPart,
	     GridPartType &gridPartG,
             const ModelTypeG& implicitModel,
	     const ModelTypeG& setupimplicitModel, const ModelTypeS& implicitModelS )
    : implicitModel_( implicitModel ),
      setupimplicitModel_( setupimplicitModel ),
      implicitModelS_( implicitModelS ),
      gridPart_( gridPart ),
      gridPartG_( gridPartG ),
      discreteSpace_( gridPart_ ),
      discreteSpaceG_( gridPartG_ ),
      discreteSpaceS_( gridPart_ ),
      discreteSpaceSS_( gridPartG_ ),
      solution_( "solution", discreteSpace_ ),
      solutionGrid_( "solutionGrid", discreteSpaceG_ ),
      rhs_( "rhs", discreteSpace_ ),
      solutionS_( "solution", discreteSpaceS_ ),
      solutionSS_( "solution", discreteSpaceSS_ ),
      rhsS_( "rhs", discreteSpaceS_ ),
      // the elliptic operator (implicit)
      implicitOperator_( implicitModel_, discreteSpace_ ),
      setupimplicitOperator_( setupimplicitModel_, discreteSpace_ ),
      implicitOperatorS_( implicitModelS_, discreteSpaceS_ ),
      // create linear operator (domainSpace,rangeSpace)
      linearOperator_( "assembled elliptic operator", discreteSpace_, discreteSpace_ ),
      linearOperatorS_( "assembled elliptic operator", discreteSpaceS_, discreteSpaceS_ ),
      // tolerance for iterative solver
      solverEps_( Dune::Fem::Parameter::getValue< double >( "poisson.solvereps", 1e-8 ) )
  {
    // set all DoF to zero
    solution_.clear();
    solutionS_.clear();
  }

  DiscreteFunctionType &solution()
  {
    return solution_;
  }
  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  DiscreteFunctionTypeG &solutionGrid()
  {
    return solutionGrid_;
  }
  const DiscreteFunctionTypeG &solutionGrid() const
  {
    return solutionGrid_;
  }
  
    DiscreteFunctionTypeS &solutionS()
  {
    return solutionS_;
  }
  const DiscreteFunctionTypeS &solutionS() const
  {
    return solutionS_;
  }

    DiscreteFunctionTypeS &solutionSS()
  {
    return solutionSS_;
  }
  const DiscreteFunctionTypeS &solutionSS() const
  {
    return solutionSS_;
  }


  //! setup the right hand side
  void prepare()
  {
    // assemble rhs
    assembleRHS ( implicitModel_, implicitModel_.rightHandSide(), implicitModel_.neumanBoundary(), rhs_ );
    // set boundary values to the rhs
    implicitOperator_.prepare( implicitModel_.dirichletBoundary(), rhs_ );
  }

void Setup_solve()
  {
   setupimplicitOperator_.jacobian( solution_ , linearOperator_); 
   // set up initial condition to satify dirichlet b.c. which is needed for iterative solvers
    setupimplicitOperator_.prepare( rhs_, solution_ ); 
    // inverse operator using linear operator
     LinearInverseOperatorType invOp( linearOperator_, solverEps_, solverEps_ );
    // solve system
     invOp( rhs_, solution_ );
  }

  void solve ( bool assemble )
  {
 
    //! [Solve the system]
    if( assemble )
    {  
      // assemble linear operator (i.e. setup matrix)
      implicitOperator_.jacobian( solution_ , linearOperator_); // adam: solution_ is u_old!
    }
    // set up initial condition to satify dirichlet b.c. which is needed for iterative solvers
    implicitOperator_.prepare( rhs_, solution_ ); // adam: one found in elliptic not above
    // inverse operator using linear operator
     LinearInverseOperatorType invOp( linearOperator_, solverEps_, solverEps_ );
    // solve system
     invOp( rhs_, solution_ );
    //! [Solve the system]
  updater ( implicitModel_, implicitModel_.rightHandSide(), solution_, solutionGrid_);

  }

  
    void solveS ( bool assemble )
  {
    //! [Solve the system]
    if( assemble )
    {  
      // assemble linear operator (i.e. setup matrix)
      implicitOperatorS_.jacobian( solutionS_ , linearOperatorS_); // adam: solution_ is u_old!
    }
    // set up initial condition to satify dirichlet b.c. which is needed for iterative solvers
   
   implicitOperatorS_.prepare( rhsS_, solutionS_ ); // adam: one found in elliptic not above
    // inverse operator using linear operator
     LinearInverseOperatorTypeS invOp( linearOperatorS_, solverEps_, solverEps_ );
    // solve system
     
     invOp( rhsS_, solutionS_ );
    //! [Solve the system]

 updaterS ( solutionS_, solutionSS_);

  }

protected:
  const ModelTypeG& implicitModel_;   // the mathematical model
  const ModelTypeG& setupimplicitModel_; 
  const ModelTypeS& implicitModelS_; 

  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with
  GridPartType  &gridPartG_;

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionSpaceTypeG discreteSpaceG_; 
  DiscreteFunctionSpaceTypeS discreteSpaceS_;
  DiscreteFunctionSpaceTypeS discreteSpaceSS_;
  DiscreteFunctionType solution_;   // the unknown
  DiscreteFunctionTypeG solutionGrid_;
  DiscreteFunctionType rhs_;        // the right hand side
  DiscreteFunctionTypeS solutionS_;   // the unknown
  DiscreteFunctionTypeS solutionSS_;   // phase field projected to gridG
  DiscreteFunctionTypeS rhsS_;        // the right hand side

  EllipticOperatorType implicitOperator_; // the implicit operator
  EllipticOperatorType setupimplicitOperator_; // the implicit operator
  EllipticOperatorTypeS implicitOperatorS_; 


  LinearOperatorType linearOperator_;  // the linear operator (i.e. jacobian of the implicit)
  LinearOperatorTypeS linearOperatorS_; 

  const double solverEps_ ; // eps for linear solver
};

#endif // end #if ELLIPT_FEMSCHEME_HH
