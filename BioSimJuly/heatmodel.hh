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
#ifndef HEAT_MODEL_HH
#define HEAT_MODEL_HH

#include <cassert>
#include <cmath>

#include <iostream>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

#include "temporalprobleminterface.hh"
#include "model.hh"



template< class FunctionSpace, class GridPart, class VeloDF, class PhaseF >
struct HeatModelS : protected DiffusionModel<FunctionSpace,GridPart>
{
  typedef DiffusionModel<FunctionSpace,GridPart> BaseType;
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

// the types in heat.hh
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef TemporalProblemInterface< FunctionSpaceType > ProblemType ;

  typedef typename BaseType::ProblemType InitialFunctionType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  typedef VeloDF DeformationDiscreteFunctionType;
  typedef typename DeformationDiscreteFunctionType :: LocalFunctionType LocalVelocityFunctionType;
  typedef typename VeloDF::DiscreteFunctionSpaceType DeformationDiscreteFunctionSpaceType;
  typedef typename LocalVelocityFunctionType::RangeType                    RangeRangeType;
  typedef typename LocalVelocityFunctionType::JacobianRangeType            RangeJacobianRangeType;


typedef PhaseF DeformationDiscreteFunctionTypeS;
  typedef typename DeformationDiscreteFunctionTypeS :: LocalFunctionType LocalReferenceFunctionTypeS;
  typedef typename PhaseF::DiscreteFunctionSpaceType ReferenceDiscreteFunctionSpaceTypeS;
  typedef typename LocalReferenceFunctionTypeS::RangeType                    RangeRangeTypeS;
  typedef typename LocalReferenceFunctionTypeS::JacobianRangeType            RangeJacobianRangeTypeS;

  static const int dimRange = FunctionSpaceType::dimRange;

  //! constructor taking problem reference, time provider,
  //! time step factor( either theta or -(1-theta) ),
  //! flag for the right hand side
  HeatModelS( const ProblemType& problem,
             const GridPart &gridPart,
	     DeformationDiscreteFunctionType &velocity,
	     DeformationDiscreteFunctionTypeS &phasef,
             const bool implicit )
    : BaseType(problem,gridPart),
      timeProvider_(problem.timeProvider()),
      implicit_( implicit ),
      phasef_( phasef ),
      vel_(velocity),
      L_0(Dune::Fem::Parameter::getValue< double >("blebbing.l_0", 1. ) ),
      timeStepFactor_( 0 )
  {
    // get theta for theta scheme
    const double theta = Dune::Fem::Parameter::getValue< double >("heat.theta", 0.5 );
    if (implicit)
      timeStepFactor_ = theta ;
    else
      timeStepFactor_ = -( 1.0 - theta ) ;
  }

  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &flux ) const
  {
//    linSource( value, gradient, entity, x, value, gradient, flux );
    // the explicit model should also evaluate the RHS
     const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
      // evaluate right hand side
 	RangeType Dwell;
	double l_B=0.025, eps_=0.5;
/**/double dist=0., norm=0.;

    const LocalVelocityFunctionType startLocal = vel_.localFunction(entity);
    RangeRangeType xi;
    startLocal.evaluate(x,xi);

 RangeRangeType y0, y1, y2;
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(0) ) ,y0);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(1) ) ,y1);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(2) ) ,y2);	
const DomainType v1 = y0-y1;
const DomainType v2 = y0-y2; 
    DomainType u_c;
    u_c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    u_c[1]=-(v1[0]*v2[2]-v1[2]*v2[0]);
    u_c[2]=v1[0]*v2[1]-v1[1]*v2[0];
double u_c_norm = std::sqrt(u_c[0]*u_c[0]+u_c[1]*u_c[1]+u_c[2]*u_c[2]);
u_c /= u_c_norm;

    u_c[0]= xGlobal[0]-L_0*u_c[0];
    u_c[1]= xGlobal[1]-L_0*u_c[1];
    u_c[2]= xGlobal[2]-L_0*u_c[2];


   for (unsigned int i=0;i<3;++i)
	{
    dist = xi[i]-u_c[i];
	norm+= dist*dist;
	}
   norm= std::sqrt(norm);
flux =norm;


  }

  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const JacobianRangeType &gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   const JacobianRangeType &gradient,
                   RangeType &flux ) const
  {
 /*   const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType m, dist;
  if( implicit_ )
    {
    const LocalVelocityFunctionType velLocal = vel_.localFunction(entity);
    RangeRangeType xi;
    velLocal.evaluate(x,xi);
    flux *= timeProvider_.deltaT();
    }
*/
    // add term from time derivative
    flux = value;
  }

   //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
   // linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  //! return the diffusive flux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar,
                          const JacobianRangeType &gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
 /* */  if( implicit_ )
    {
    double eps_=0.5;
 //   flux = gradient;
    flux *=eps_;
    }
    flux *= timeProvider_.deltaT(); 
  }
  template< class Entity, class Point >
  void alpha(const Entity &entity, const Point &x, 
             const RangeType &value,
             RangeType &val) const
  {
    BaseType::alpha(entity,x,value,val);
  }
  template< class Entity, class Point >
  void linAlpha(const RangeType &uBar, 
                const Entity &entity, const Point &x, 
                const RangeType &value,
                RangeType &val) const
  {
    BaseType::linAlpha(uBar,entity,x,value,val);
  }
  //! exact some methods from the problem class
  bool hasDirichletBoundary () const
  {
    return BaseType::hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return BaseType::hasNeumanBoundary() ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return BaseType::isDirichletIntersection(inter,dirichletComponent) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar,
          const Entity &entity,
          const Point &x,
          RangeType &u ) const
  {
    BaseType::g(uBar,entity,x,u);
  }

  // return Fem::Function for Dirichlet boundary values
  typename BaseType::DirichletBoundaryType dirichletBoundary( ) const
  {
    return BaseType::dirichletBoundary();
  }

  //! return reference to Problem's time provider
  const TimeProviderType & timeProvider() const
  {
    return timeProvider_;
  }

  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }

  DeformationDiscreteFunctionType &velocity()
  {
    return vel_;
  }
 /* const DeformationDiscreteFunctionType &velocity() const
  {
    return vel_;
  }
*/
protected:
  using BaseType::problem_;
  const TimeProviderType &timeProvider_;
  bool implicit_;
  double timeStepFactor_;
  const double L_0;

private:
 DeformationDiscreteFunctionType &vel_;
 DeformationDiscreteFunctionTypeS &phasef_;

};

template< class FunctionSpace, class FunctionSpaceG, class GridPart, class surfaceRef, class PhaseF  >
struct HeatModelG : protected DiffusionModel<FunctionSpace,GridPart>
{
  typedef DiffusionModel<FunctionSpace,GridPart> BaseType;
  typedef FunctionSpace FunctionSpaceType;
  typedef FunctionSpaceG FunctionSpaceTypeG;		// add
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef TemporalProblemInterface< FunctionSpaceType > ProblemType ;

  typedef typename BaseType::ProblemType InitialFunctionType;

  typedef Dune::Fem::TimeProviderBase TimeProviderType;

  typedef surfaceRef DeformationDiscreteFunctionType;
  typedef typename DeformationDiscreteFunctionType :: LocalFunctionType LocalReferenceFunctionType;
  typedef typename surfaceRef::DiscreteFunctionSpaceType ReferenceDiscreteFunctionSpaceType;
  typedef typename LocalReferenceFunctionType::RangeType                    RangeRangeType;
  typedef typename LocalReferenceFunctionType::JacobianRangeType            RangeJacobianRangeType;

typedef PhaseF DeformationDiscreteFunctionTypeS;
  typedef typename DeformationDiscreteFunctionTypeS :: LocalFunctionType LocalReferenceFunctionTypeS;
  typedef typename PhaseF::DiscreteFunctionSpaceType ReferenceDiscreteFunctionSpaceTypeS;
  typedef typename LocalReferenceFunctionTypeS::RangeType                    RangeRangeTypeS;
  typedef typename LocalReferenceFunctionTypeS::JacobianRangeType            RangeJacobianRangeTypeS;

  static const int dimRange = FunctionSpaceType::dimRange;

  //! constructor taking problem reference, time provider,
  //! time step factor( either theta or -(1-theta) ),
  //! flag for the right hand side
  HeatModelG( const ProblemType& problem,
             const GridPart &gridPart,
	     DeformationDiscreteFunctionType &cortex,
	     DeformationDiscreteFunctionType &start,
	     DeformationDiscreteFunctionTypeS &phasef,
             const bool implicit,
	     const bool setup,
	     double &area,
	     double &vol )
    : BaseType(problem,gridPart),
      timeProvider_(problem.timeProvider()),
      implicit_( implicit ),
      setup_( setup),
      area_( area ),
      vol_( vol ),
      start_( start ),
      cortex_( cortex ),
      phasef_( phasef ),
      M_KA(Dune::Fem::Parameter::getValue< double >("blebbing.m_ka", 1. ) ),
      M_KB(Dune::Fem::Parameter::getValue< double >("blebbing.m_kb", 1. ) ),
      X_0(Dune::Fem::Parameter::getValue< double >("blebbing.x_0", 1. ) ),
      P_0(Dune::Fem::Parameter::getValue< double >("blebbing.p_0", 1. ) ),
      K_0(Dune::Fem::Parameter::getValue< double >("blebbing.k_0", 1. ) ),
      U_B(Dune::Fem::Parameter::getValue< double >("blebbing.u_b", 1. ) ),
      L_0(Dune::Fem::Parameter::getValue< double >("blebbing.l_0", 1. ) ),
      Omega(Dune::Fem::Parameter::getValue< double >("blebbing.omega", 2.1677e-08 ) ),
      timeStepFactor_( 0 )
  {
    // get theta for theta scheme
    const double theta = Dune::Fem::Parameter::getValue< double >("heat.theta", 0.5 );
    if (implicit)
      timeStepFactor_ = theta ;
    else
      timeStepFactor_ = -( 1.0 - theta ) ;
  }
  
 
  template< class Entity, class Point >
  void source ( const Entity &entity,
                const Point &x,
                const RangeType &value,
                const JacobianRangeType &gradient,
                RangeType &flux ) const
  {
	// time mass/setupidentity
	for( int localCol = 0; localCol < 3; ++localCol )
	{
	flux[localCol] += value[localCol];	// u +Mu
	}

if (!setup_)
	{
double norm=0.;
int N=3;
RangeType nu;
const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
 const LocalReferenceFunctionType cortexLocal = cortex_.localFunction(entity);
 const LocalReferenceFunctionType startLocal = start_.localFunction(entity);
 RangeRangeType xc0, xc1, xc2;
 
cortexLocal.evaluate( entity.geometry().local( entity.geometry().corner(0) ) ,xc0);	
cortexLocal.evaluate( entity.geometry().local( entity.geometry().corner(1) ) ,xc1);	
cortexLocal.evaluate( entity.geometry().local( entity.geometry().corner(2) ) ,xc2);	

const DomainType w1 = xc0-xc1;
    const DomainType w2 = xc0-xc2; 
    DomainType nu2;
    nu2[0]=w1[1]*w2[2]-w1[2]*w2[1];
    nu2[1]=-(w1[0]*w2[2]-w1[2]*w2[0]);
    nu2[2]=w1[0]*w2[1]-w1[1]*w2[0];
  double dirc=0.;  
		for( int localCol = 0; localCol < N; ++localCol )
		{
		norm += nu2[localCol]*nu2[localCol];
		nu[localCol]=nu2[localCol];
		dirc += nu2[localCol]*xGlobal[localCol];
		}
	
	if (dirc <0.)
	{
	nu *=-1.;
	}
		
nu/=std::sqrt(norm);

double dist=0., norm1=0., norm2=0.;
  
    RangeRangeType xi;
    cortexLocal.evaluate(x,xi);
RangeType binders(0);

 RangeRangeType y0, y1, y2;
 
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(0) ) ,y0);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(1) ) ,y1);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(2) ) ,y2);	
const DomainType v1 = y0-y1;
const DomainType v2 = y0-y2; 
    DomainType u_c;
    u_c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    u_c[1]=-(v1[0]*v2[2]-v1[2]*v2[0]);
    u_c[2]=v1[0]*v2[1]-v1[1]*v2[0];
double u_c_norm = std::sqrt(u_c[0]*u_c[0]+u_c[1]*u_c[1]+u_c[2]*u_c[2]);
u_c /= u_c_norm;

    u_c[0]= xGlobal[0]-L_0*u_c[0];
    u_c[1]= xGlobal[1]-L_0*u_c[1];
    u_c[2]= xGlobal[2]-L_0*u_c[2];

   for (unsigned int i=0;i<3;++i)
	{
    binders[i] = value[i]-u_c[i];
	norm2 += binders[i]*binders[i];
	}
norm2=std::sqrt(norm2);

   for (unsigned int i=0;i<N;++i)
	{
     	u_c[i] +=L_0*binders[i]/(norm2+1e-10);
	}

if (norm2 < 0.02 )
{
//std::cout << "error" << std::endl;
}
 
u_c *=K_0/(1.+std::exp(2.*(norm2-U_B)/1e-10) );
u_c *= 1.+500./(1+std::exp(2.*(norm2-0.03)/1e-10)); // fail safe
u_c /= Omega;
u_c *= timeProvider_.deltaT();
flux+=u_c;

//nu /=3*vol_;	// 2D
nu /=vol_+1e-10;
nu *=  P_0; // 2.4572e-05; //1e-05;	// WARNING
nu /= Omega;

nu *= timeProvider_.deltaT();


flux+=nu;
}
  }

  template< class Entity, class Point >
  void linSource ( const RangeType& uBar,
                   const JacobianRangeType &gradientBar,
                   const Entity &entity,
                   const Point &x,
                   const RangeType &value,
                   const JacobianRangeType &gradient,
                   RangeType &flux ) const
  {
	flux=value;
	if (!setup_){
		double norm2=0.;
		const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
 		const LocalReferenceFunctionType cortexLocal = cortex_.localFunction(entity);
 		const LocalReferenceFunctionType startLocal = start_.localFunction(entity);
	    RangeRangeType xi;
    	cortexLocal.evaluate(x,xi);
		RangeType binders(0), Linker(0);

 		RangeRangeType y0, y1, y2;
 
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(0) ) ,y0);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(1) ) ,y1);	
startLocal.evaluate( entity.geometry().local( entity.geometry().corner(2) ) ,y2);	
const DomainType v1 = y0-y1;
const DomainType v2 = y0-y2; 
    DomainType u_c;
    u_c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    u_c[1]=-(v1[0]*v2[2]-v1[2]*v2[0]);
    u_c[2]=v1[0]*v2[1]-v1[1]*v2[0];
double u_c_norm = std::sqrt(u_c[0]*u_c[0]+u_c[1]*u_c[1]+u_c[2]*u_c[2]);
u_c /= u_c_norm;

    u_c[0]= xGlobal[0]-L_0*u_c[0];
    u_c[1]= xGlobal[1]-L_0*u_c[1];
    u_c[2]= xGlobal[2]-L_0*u_c[2];

   for (unsigned int i=0;i<3;++i)
	{
    binders[i] = uBar[i]-u_c[i];
	norm2 += binders[i]*binders[i];
	Linker[i]=value[i];
	}
norm2=std::sqrt(norm2);

Linker *=K_0/(1.+std::exp(2.*(norm2-U_B)/1e-10) );
Linker *= 1.+500./(1+std::exp(2.*(norm2-0.03)/1e-10)); // fail safe
Linker /= Omega;
Linker *= timeProvider_.deltaT();
flux+=Linker;
	
	}
	
  }
  

   //! return the diffusive flux
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity,
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
flux=0.;
double para=2.*area_;
//para /=161.; //455.;
para = 0.3;
int N=3;
if (!setup_)
	{
		
    const LocalReferenceFunctionType cortexLocal = cortex_.localFunction(entity);
    RangeJacobianRangeType dxc;
    RangeRangeType xc;
    cortexLocal.jacobian(x,dxc);
    cortexLocal.evaluate(x,xc);	
	
	const LocalReferenceFunctionType StartLocal = start_.localFunction(entity);
    RangeJacobianRangeType dxt;
    RangeRangeType xt;
    StartLocal.jacobian(x,dxt);
    StartLocal.evaluate(x,xt);	
	flux=0.;
	
	double norm =0., normc=0., norm1=0.;
	double m_ka= 17e-6;
		for( int localCol = 0; localCol < N; ++localCol )
		{
		 for( int localRow = 0; localRow < 3; ++localRow )
		 {
		 normc += dxt[localCol][localRow]*dxt[localCol][localRow];
		 norm  += gradient[localCol][localRow]*gradient[localCol][localRow];
		 //normc  += dxc[localCol][localRow]*dxc[localCol][localRow];
		 }
		}
//std::cout << normc << "			"	<<  norm << std::endl;
   const LocalReferenceFunctionTypeS phasefLocal = phasef_.localFunction(entity);
    RangeRangeTypeS c;
    phasefLocal.evaluate(x,c);	
   
	for( int localCol = 0; localCol < N; ++localCol )
		{
	for( int localRow = 0; localRow < 3; ++localRow )
		 {
		flux[localCol][localRow] += gradient[localCol][localRow];
		flux[localCol][localRow] *= std::sqrt(normc)/(std::sqrt(norm)+1e-10);
		flux[localCol][localRow] *= X_0;
		//flux[localCol][localRow] *= std::sqrt(norm);
		//flux[localCol][localRow] *= para;
		flux[localCol][localRow] *= M_KA;// m_ka; // /2.;	
		flux[localCol][localRow] *= timeProvider_.deltaT();	// hit both u terms
		flux[localCol][localRow] /= Omega;
		}
		}
 // std::cout << para << std::endl;
	}

 if (setup_)
 {
	for( int localCol = 0; localCol < N; ++localCol )
	{
	flux[localCol+3]=gradient[localCol];
	//flux[localCol+3]*=para*para;	
	}
 }

  }
  
  //! return the diffusive flux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar,
                          const JacobianRangeType &gradientBar,
                          const Entity &entity,
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
flux=0;
	if (!setup_)
	{
	double m_ka= 17e-6, norm=0.;
	JacobianRangeType bend(0);

    const LocalReferenceFunctionType cortexLocal = cortex_.localFunction(entity);
    RangeJacobianRangeType dxc;
    RangeRangeType xc;
    cortexLocal.jacobian(x,dxc);
    cortexLocal.evaluate(x,xc);	

   const LocalReferenceFunctionTypeS phasefLocal = phasef_.localFunction(entity);
    RangeRangeTypeS c;
    phasefLocal.evaluate(x,c);	
//std::cout << c << std::endl;
int N=3;
double para=2.*area_;
//para /=161.; //455.;
para=0.3;
	for( int localCol = 0; localCol < N; ++localCol )
		{
	for( int localRow = 0; localRow < 3; ++localRow )
		 {
		bend[localCol][localRow] += gradient[localCol][localRow];
		bend[localCol][localRow] *= M_KA; //m_ka; // /2.;
		//bend[localCol][localRow] *= para*para;
		}
		}

	
	for( int localCol = 0; localCol < N; ++localCol )
		{
	for( int localRow = 0; localRow < 3; ++localRow )
		 {
		flux[localCol][localRow] -= gradient[localCol+3][localRow];	
		flux[localCol][localRow] *= M_KB; // 0.00003690547; //8e-13;	//0.14e-6/2.;
		//flux[localCol][localRow] *= para*para;
		flux[localCol+3][localRow] += gradient[localCol][localRow];
		//flux[localCol+3][localRow] *= para*para;
		flux[localCol][localRow] += bend[localCol][localRow];
		flux[localCol][localRow] *= timeProvider_.deltaT();	// hit both u terms
		flux[localCol][localRow] /= Omega;
		}
		}
	}

  }
  
    //! return the diffusive flux

  
  template< class Entity, class Point >
  void alpha(const Entity &entity, const Point &x, 
             const RangeType &value,
             RangeType &val) const
  {
    BaseType::alpha(entity,x,value,val);
  }
  template< class Entity, class Point >
  void linAlpha(const RangeType &uBar, 
                const Entity &entity, const Point &x, 
                const RangeType &value,
                RangeType &val) const
  {
    BaseType::linAlpha(uBar,entity,x,value,val);
  }
  //! exact some methods from the problem class
  bool hasDirichletBoundary () const
  {
    return BaseType::hasDirichletBoundary() ;
  }
  bool hasNeumanBoundary () const
  {
    return BaseType::hasNeumanBoundary() ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter, Dune::FieldVector<bool,dimRange> &dirichletComponent ) const
  {
    return BaseType::isDirichletIntersection(inter,dirichletComponent) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar,
          const Entity &entity,
          const Point &x,
          RangeType &u ) const
  {
    BaseType::g(uBar,entity,x,u);
  }

  // return Fem::Function for Dirichlet boundary values
  typename BaseType::DirichletBoundaryType dirichletBoundary( ) const
  {
    return BaseType::dirichletBoundary();
  }

  //! return reference to Problem's time provider
  const TimeProviderType & timeProvider() const
  {
    return timeProvider_;
  }

  const InitialFunctionType &initialFunction() const
  {
    return problem_;
  }

  typename BaseType::RightHandSideType rightHandSide(  ) const
  {
    return BaseType::rightHandSide();
  }

protected:
  using BaseType::problem_;
  const TimeProviderType &timeProvider_;
  bool implicit_;
  bool setup_;
  double timeStepFactor_;
  const double M_KA, M_KB, X_0, K_0, P_0, Omega, U_B, L_0;

private:
  double &area_;
  double &vol_;
  DeformationDiscreteFunctionType &cortex_;
  DeformationDiscreteFunctionType &start_;
 DeformationDiscreteFunctionTypeS &phasef_;
};

#endif // #ifndef HEAT_MODEL_HH
