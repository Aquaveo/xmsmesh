//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing
/// \copyright (C) Copyright Aquaveo 2018. Distributed under FreeBSD License
/// (See accompanying file LICENSE or https://aqaveo.com/bsd/license.txt)
//------------------------------------------------------------------------------
#pragma once

//----- Included files ---------------------------------------------------------
#include <vector>
#include <xmscore/stl/vector.h>
#include <xmscore/misc/base_macros.h>
#include <xmscore/misc/boost_defines.h>

//----- Forward declarations ---------------------------------------------------

//----- Namespace declaration --------------------------------------------------
namespace xms
{
class MePolyOffsetterOutput;
class InterpBase;
//----- Constants / Enumerations -----------------------------------------------

//----- Structs / Classes ------------------------------------------------------
/// \brief Redistributes the point locations on a polygon based on a size
/// function
/// \see MePolyRedistributePtsImpl
class MePolyRedistributePts
{
public:
  static BSHP<MePolyRedistributePts> New();
  virtual ~MePolyRedistributePts();

  /// \cond
  virtual void SetSizeFunc(BSHP<InterpBase> a_interp) = 0;
  virtual void SetSizeFuncFromPoly(const VecPt3d& a_outPoly,
                                   const VecPt3d2d& a_inPolys,
                                   double a_sizeBias) = 0;
  virtual void SetConstantSizeFunc(double a_size) = 0;
  virtual void SetConstantSizeBias(double a_sizeBias) = 0;
  virtual void SetUseCurvatureRedistribution(double a_featureSize,
                                             double a_meanSpacing,
                                             double a_minimumCurvature,
                                             bool a_smooth) = 0;
  virtual VecPt3d Redistribute(const VecPt3d& a_polyLine) = 0;
  virtual double SizeFromLocation(const Pt3d& a_location) = 0;

  virtual std::string ToPyRepr() const = 0;

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MePolyRedistributePts);
  /// \endcond

protected:
  MePolyRedistributePts();
};
//----- Function prototypes ----------------------------------------------------

} // namespace xms
