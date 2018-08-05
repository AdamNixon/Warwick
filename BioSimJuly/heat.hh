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
#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include <iostream>

#include "temporalprobleminterface.hh"

template <class FunctionSpace>
class TimeDependentCosinusProduct : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;

  TimeDependentCosinusProduct( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0); 
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    // no exact solution - method needed for initial data
    phi = RangeType(0);
    const double r1 = std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	
    phi[0] = x[0];
    phi[1] = x[1];
    phi[2] = x[2];
   // phi[3] = 2.*x[0]/r1;
   // phi[4] = 2.*x[1]/r1;
   // phi[5] = 2.*x[2]/r1;
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // no exact solution
    ret = JacobianRangeType(0);
  }
    //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }
};

template <class FunctionSpace>
class TimeDependentSurfaceProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;

  TimeDependentSurfaceProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
                 RangeType& phi) const
  {
    phi = RangeType(0); 
  }

  //! the exact solution
  virtual void u(const DomainType& x,
                 RangeType& phi) const
  {
    // no exact solution - method needed for initial data
double n0, n1, n2;
n0 = std::sqrt(x[0]*x[0]);
n1 = std::sqrt(x[1]*x[1]);
n2 = std::sqrt(x[2]*x[2]);

phi = RangeType(0); 
/*
double a=3., b=0.125;
if (x[1]>a){
phi= std::cos(b*M_PI*(n1-a))*std::cos(b*M_PI*(n1-a));
}
if (x[1]>11.){
phi = RangeType(1);
} */
//phi = sin(M_PI*x[0]/4.)*sin(M_PI*x[0]/4.);

  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
                         JacobianRangeType& ret) const
  {
    // no exact solution
    ret = JacobianRangeType(0);
  }
    //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }
};

#endif // #ifndef POISSON_HH
