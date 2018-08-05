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
#ifndef DEFORMATION_HH
#define DEFORMATION_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/geometrygrid/coordfunction.hh>
#include <dune/fem/space/common/functionspace.hh>

// DeformationCoordFunction
// ------------------------

template< int dimWorld >
struct DeformationCoordFunction
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunction ( const double time = 0.0 )
  : time_( time )
  {}

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    //const double newtime = (time_ < 0.5 ? 0.0 : (time_ < 0.75 ? 4*(time_ - 0.5) : 1.0));
 /*    const double newtime = std::min( time_+0.75, 1.0 );

   
    const double r1 = std::abs( x[ 0 ] );
    const double target = (1.0 - (r1*r1))*((r1*r1) + 0.05) + (r1*r1)*sqrt(1.0 - (r1*r1));

    const double r2 = std::sqrt( x[1]*x[1] + x[2]*x[2] );
    const double factor = std::exp( -2*newtime )*r2 + (1.0 - std::exp( -2*newtime ))*target;

    y[ 0 ] = 2 * x[ 0 ] + newtime*(x[ 0 ] > 0 ? 2.0 : -1.0 )*x[ 0 ];
    y[ 1 ] = factor * x[ 1 ] / (r2 + 0.000001);
    y[ 2 ] = factor * x[ 2 ] / (r2 + 0.000001); */ 
  /* */
   const double r1 = std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    const double r0 = 4.;
 const double r2 = std::sqrt(x[0]*x[0]+x[1]*x[1]);
   
 y[0] = r0*x[ 0 ]/r1;
   y[1] = r0*x[ 1 ]/r1;
   y[2] = r0*x[ 2 ]/r1;

  y[0] = r0*x[0]/r1;
  y[1] = r0*x[1]/r1;


const double r3 = std::sqrt(y[0]*y[0]+y[1]*y[1]);

//y[2]=r3;
/* */
y[2] = 1.5-0.5*std::cos(M_PI*r3/2.);

if(r3>2.)
{y[2]=std::sqrt(4.-(r3-2.)*(r3-2.));}

if(x[2]<0.)
{y[2]=-y[2];} 

//y *= 7./4;

/*
y[2]= r0*x[2]/r1;
if(y[2]>0){	
	if(r3<2.5)
	{
	//y[2]= y[2]+std::sqrt(4.-r3*r3);
	y[2]= std::sqrt(r0*r0-2.5*2.5)+1.5+1.5*std::cos(M_PI*pow(r3*r3/(2.5*2.5),0.75));	
	}
}
  */
  /*
  y[0] = r0*x[0]/r1;
  y[1] = (0.3*y[0]*y[0]+0.6*r0)*x[1]/r1;
  y[2] = (0.3*y[0]*y[0]+0.6*r0)*x[2]/r1; */
 /*  
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];  
 */
  }

  void setTime ( const double time ) { time_ = time; }

private:
  double time_;
};

//! deformation depending on a discrete function
template <class DiscreteFunctionType>
class DeformationDiscreteFunction
: public Dune::DiscreteCoordFunction< double, 3, DeformationDiscreteFunction< DiscreteFunctionType > >
{
  typedef Dune::DiscreteCoordFunction< double, 3, DeformationDiscreteFunction< DiscreteFunctionType > > BaseType;

  typedef typename DiscreteFunctionType :: GridType GridType ;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType ;
  typedef typename DiscreteFunctionType :: RangeType  RangeType ;
public:
  DeformationDiscreteFunction ( const DiscreteFunctionType& vertices )
  : vertices_( vertices )
  {}

  template< class HostEntity , class RangeVector >
  void evaluate ( const HostEntity &hostEntity, unsigned int corner,
                  RangeVector &y ) const
  {
    DUNE_THROW(Dune::NotImplemented,"evaluate not implemented for codim > 0");
  }

  template <class RangeVector>
  void evaluate ( const typename GridType :: template Codim<0>::Entity &entity,
		  unsigned int corner,
                  RangeVector &y ) const
  {
    y = entity.geometry()[corner];

    return;
    typedef typename GridType::ctype  ctype;
    enum { dim = GridType::dimension };

    const Dune::ReferenceElement< ctype, dim > &refElement
      = Dune::ReferenceElements< ctype, dim >::general( entity.type() );

    LocalFunctionType localVertices = vertices_.localFunction( entity );

    localVertices.evaluate( refElement.position( corner, dim ), y );
  }

  void setTime ( const double time )
  {
  }

protected:
  const DiscreteFunctionType& vertices_;
};

#endif // #ifndef DEFORMATION_HH

