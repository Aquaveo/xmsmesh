#pragma once
//------------------------------------------------------------------------------
/// \file
/// \ingroup meshing_detail
/// \copyright (C) Copyright Aquaveo 2018. Distributed under the xmsng
///  Software License, Version 1.0. (See accompanying file
///  LICENSE_1_0.txt or copy at http://www.aquaveo.com/xmsng/LICENSE_1_0.txt)
//------------------------------------------------------------------------------

//----- Included files ---------------------------------------------------------

// 3. Standard library headers
#include <vector>

// 4. External library headers

// 5. Shared code headers
#include <xmscore/stl/vector.h>
#include <xmscore/misc/base_macros.h> // for XM_DISALLOW_COPY_AND_ASSIGN
#include <xmscore/misc/boost_defines.h>

//----- Forward declarations ---------------------------------------------------

namespace xms
{
//----- Forward declarations ---------------------------------------------------

class TrTin;
class MePolyRedistributePts;

////////////////////////////////////////////////////////////////////////////////
class MeRelaxer
{
public:
  static BSHP<MeRelaxer> New();

  /// \cond
  virtual void Relax(const VecInt& a_fixedPoints, boost::shared_ptr<TrTin> a_tin) = 0;
  virtual bool SetRelaxationMethod(const std::string& a_relaxType) = 0;
  virtual void SetPointSizer(BSHP<MePolyRedistributePts> a_sizer) = 0;
  /// \endcond

protected:
  MeRelaxer();
  virtual ~MeRelaxer();

private:
  XM_DISALLOW_COPY_AND_ASSIGN(MeRelaxer)
}; // MeRelaxer

} // namespace xms
