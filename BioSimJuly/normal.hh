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
#ifndef RHS_HH
#define RHS_HH

#include <iostream>

#include <dune/fem/quadrature/cachingquadrature.hh>


// assembleRHS
// -----------

template< class Model, class Function, class DiscreteFunction, class DiscreteFunctionG >
void updater ( const Model &model, const Function &function,  DiscreteFunction &sol, DiscreteFunctionG &solG )
{ 
  solG.clear();
double r=0., s=0.;
double check=0.;
  for (int i=0; i<sol.blocks(); ++i)
    {
	check+=1;
	for (int j=0; j<3; ++j)
      {
      (*solG.block(i))[j]=(*sol.block(i))[j];		// sol is (u,v) \Delta u=w
      }
	r+=std::sqrt((*sol.block(i))[0]*(*sol.block(i))[0]+(*sol.block(i))[1]*(*sol.block(i))[1]+(*sol.block(i))[2]*(*sol.block(i))[2]);
    }   
    std::cout << r/check << "	number of points:  " << check << std::endl;
  sol.communicate();
}

template< class DiscreteFunction, class DiscreteFunctionS >
void updaterS ( DiscreteFunction &sol, DiscreteFunctionS &solS )
{ 
  solS.clear();

double r=0., s=0.;
  for (int i=0; i<sol.blocks(); ++i)
    {
      (*solS.block(i))[0]=(*sol.block(i))[0];		// sol is (u,v) \Delta u=w
r+=(*sol.block(i))[0];
s+=1.;
    }   
//     std::cout << "phase normal:  " << r/s << std::endl;
  solS.communicate();
}

template< class Function, class DiscreteFunction >
void assemble2DVolume ( const Function &function,  DiscreteFunction &sol, double &volume )
{ 
 volume=0.;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  const DiscreteFunctionSpaceType &dfSpace = sol.space();

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    typedef typename Function::RangeType RangeType;
    
RangeType x0, x1, x2;
 sol.evaluate(entity.geometry().corner(0),x0);	
sol.evaluate(entity.geometry().corner(1),x1);	
sol.evaluate(entity.geometry().corner(2),x2);	

   double det=0.; 
  det+=x0[0]*x1[1]*x2[2];
  det+=x1[0]*x2[1]*x0[2];
  det+=x2[0]*x0[1]*x1[2];
  det-=x0[2]*x1[1]*x2[0];
  det-=x1[2]*x2[1]*x0[0];
  det-=x2[2]*x0[1]*x1[0];
det/=2.;
volume+=std::abs(det);


  }
   std::cout << volume << std::endl;
}

template< class DiscreteFunction >
void assembleArea ( const DiscreteFunction solution, double &area )
{ 
  area=0.;
  double volume=0.;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;

  const DiscreteFunctionSpaceType &dfSpace = solution.space();


  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    area+=entity.geometry().volume();

  }

std::cout << "area:	" << area << std::endl;
}

template< class DiscreteFunction >
void assembleVolume ( DiscreteFunction &rhs, double &volume )
{ 
volume=0.;
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;


  const DiscreteFunctionSpaceType &dfSpace = rhs.space();

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
   double det=0.; 
  det+=entity.geometry().corner(0)[0]*entity.geometry().corner(1)[1]*entity.geometry().corner(2)[2];
  det+=entity.geometry().corner(1)[0]*entity.geometry().corner(2)[1]*entity.geometry().corner(0)[2];
  det+=entity.geometry().corner(2)[0]*entity.geometry().corner(0)[1]*entity.geometry().corner(1)[2];
  det-=entity.geometry().corner(0)[2]*entity.geometry().corner(1)[1]*entity.geometry().corner(2)[0];
  det-=entity.geometry().corner(1)[2]*entity.geometry().corner(2)[1]*entity.geometry().corner(0)[0];
  det-=entity.geometry().corner(2)[2]*entity.geometry().corner(0)[1]*entity.geometry().corner(1)[0];
det/=6.;
volume+=std::abs(det);
  }
  //  normalise 
 std::cout << "volume	" << volume << std::endl;
}



template< class Function, class DiscreteFunction >
void assembleVol ( const Function &function, DiscreteFunction &sol, double &vol )
{
vol=0.;
  
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  const DiscreteFunctionSpaceType &dfSpace = sol.space();

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;
  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  const int quadOrder = 2*dfSpace.order()+1;

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    const typename Function::LocalFunctionType localFunction =
             function.localFunction( entity);
    LocalFunctionType solLocal = sol.localFunction( entity );
    typedef typename Function::RangeType RangeType;
    
    RangeType xc0, xc1, xc2;
    solLocal.evaluate( entity.geometry().local( entity.geometry().corner(0) ) ,xc0);	
    solLocal.evaluate( entity.geometry().local( entity.geometry().corner(1) ) ,xc1);	
    solLocal.evaluate( entity.geometry().local( entity.geometry().corner(2) ) ,xc2);		

    const RangeType w1 = xc0-xc1;
    const RangeType w2 = xc0-xc2; 
    const RangeType w3 = xc0+xc1+xc2;
    RangeType nu2;
    nu2[0]=w1[1]*w2[2]-w1[2]*w2[1];
    nu2[1]=-(w1[0]*w2[2]-w1[2]*w2[0]);
    nu2[2]=w1[0]*w2[1]-w1[1]*w2[0];
    double norm=0., localVol=0.;
		for( int localCol = 0; localCol < 3; ++localCol )
		{
		norm += nu2[localCol]*nu2[localCol];
		localVol +=nu2[localCol]*w3[localCol]; //nu2[localCol]*w3[localCol];
		}
    localVol /=3.; // avg of xc0 xc1 xc2
    //localVol /=std::sqrt(norm);
    localVol /=2.;  						// area on \Gamma(t) and normalising the norm cancel but for factor 2
   // localVol *=entity.geometry().volume();	// or int over gamma^0
    localVol /=3.; // Dimension of surface
    vol +=std::abs(localVol);
  }
    std::cout << "volume	" << vol << std::endl;
  sol.communicate();
}


#endif // #ifndef RHS_HH
