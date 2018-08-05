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
#ifndef HEAT_FEMSCHEME_HH
#define HEAT_FEMSCHEME_HH

#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/space/common/interpolate.hh>

// local includes
#include "normal.hh"
#include "femscheme.hh"
#include "model.hh"

// HeatScheme
//-----------

template < class ImplicitModelG, class ExplicitModel, class ImplicitModelS, class ExplicitModelS >
struct HeatScheme : public FemScheme<ImplicitModelG, ImplicitModelS>
{
  typedef FemScheme<ImplicitModelG, ImplicitModelS> BaseType;
  typedef typename BaseType::GridType GridType;
  typedef typename BaseType::GridPartType GridPartType;
  typedef typename BaseType::ModelTypeG ImplicitModelType;	// is the model used in femscheme
  typedef ExplicitModel ExplicitModelType;
  typedef typename BaseType::ModelTypeS ImplicitModelTypeS;
  typedef ExplicitModelS ExplicitModelTypeS;
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  HeatScheme( GridPartType &gridPart,
	      GridPartType &gridPartG,
              const ImplicitModelType& implicitModel,
              const ExplicitModelType& explicitModel,
	      const ImplicitModelType& setupimplicitModel,
              const ExplicitModelType& setupexplicitModel,
	      const ImplicitModelTypeS& implicitModelS,
              const ExplicitModelTypeS& explicitModelS,
	      double &area,  double &vol )
  : BaseType(gridPart, gridPartG, implicitModel, setupimplicitModel, implicitModelS),
    explicitModel_(explicitModel),
    setupexplicitModel_(setupexplicitModel),
    explicitModelS_(explicitModelS),
    explicitOperatorS_( explicitModelS_, discreteSpaceS_ ),
    explicitOperator_( explicitModel_, discreteSpace_ ),
    setupexplicitOperator_( setupexplicitModel_, discreteSpace_ ),
    area_( area ),
    vol_( vol )
  {
  }

  void prepare()
  {

   // assembleArea( solutionGrid_, area_); 
    assembleArea( solution_, area_); 
double areaG_;
assembleArea( solutionGrid_, areaG_); 
double volT_;
assembleVol (  implicitModel_.rightHandSide(), solution_,  volT_ ); // 3D
   // assembleVolume( solutionGrid_, vol_); // 2D
   // assemble2DVolume( implicitModel_.rightHandSide(), solution_, vol_);
     // apply constraints, e.g. Dirichlet contraints, to the solution
    explicitOperator_.prepare( explicitModel_.dirichletBoundary(), solution_); // in elliptic.hh  
    // apply explicit operator and also setup right hand side
    explicitOperator_( solution_, rhs_);	
    // apply constraints, e.g. Dirichlet contraints, to the result
    explicitOperator_.prepare( solution_, rhs_ ); // again elliptic.hh function

    explicitOperatorS_.prepare( explicitModelS_.dirichletBoundary(), solutionS_); 
    // apply explicit operator and also setup right hand side
     explicitOperatorS_( solutionS_, rhsS_);		
    // apply constraints, e.g. Dirichlet contraints, to the result
    explicitOperatorS_.prepare( solutionS_, rhsS_ ); 
    
 updaterS ( solutionS_, solutionSS_);
  }

 void Setup_prepare()
  {
  //  assembleArea( solutionGrid_, area_);
    assembleArea( solution_, area_); 
    assembleVol (  implicitModel_.rightHandSide(), solution_,  vol_ );
  //  assembleVolume( solutionGrid_, vol_);
 // assemble2DVolume( implicitModel_.rightHandSide(), solution_, vol_);
     // apply constraints, e.g. Dirichlet contraints, to the solution
    setupexplicitOperator_.prepare( setupexplicitModel_.dirichletBoundary(), solution_); // in elliptic.hh  
    // apply explicit operator and also setup right hand side
    setupexplicitOperator_( solution_, rhs_);	
    // apply constraints, e.g. Dirichlet contraints, to the result
    setupexplicitOperator_.prepare( solution_, rhs_ ); // again elliptic.hh function
    
  }

  void initialize ()
  {
    
    typedef Dune::Fem::GridFunctionAdapter< typename ExplicitModelType::InitialFunctionType, GridPartType > GridInitialFunction;
    GridInitialFunction gridInitialFunction( "initial data", explicitModel_.initialFunction(), solution_.gridPart() );
    interpolate( gridInitialFunction, solution_ );  
 
    typedef Dune::Fem::GridFunctionAdapter< typename ExplicitModelTypeS::InitialFunctionType, GridPartType > GridInitialFunctionS;
    GridInitialFunctionS gridInitialFunctionS( "initial data", explicitModelS_.initialFunction(), solutionS_.gridPart() );
    interpolate( gridInitialFunctionS, solutionS_ );

 updaterS ( solutionS_, solutionSS_);

  }

private:
  using BaseType::gridPart_;
  using BaseType::gridPartG_;
  using BaseType::discreteSpace_;
  using BaseType::discreteSpaceS_;
  using BaseType::solution_;
  using BaseType::solutionGrid_;
  using BaseType::solutionS_;
  using BaseType::solutionSS_;
  using BaseType::implicitModel_;
  using BaseType::setupimplicitModel_;
  using BaseType::implicitModelS_;
  using BaseType::rhs_;
  using BaseType::rhsS_;
  const ExplicitModelType &explicitModel_;
  const ExplicitModelType &setupexplicitModel_;
  const ExplicitModelTypeS &explicitModelS_;
  typename BaseType::EllipticOperatorType explicitOperator_; // the operator for the rhs
  typename BaseType::EllipticOperatorType setupexplicitOperator_;
  typename BaseType::EllipticOperatorTypeS explicitOperatorS_; 
  double &area_;
  double &vol_;
};

#endif // end #if HEAT_FEMSCHEME_HH
