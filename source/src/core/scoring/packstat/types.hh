// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/types.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_types_hh
#define INCLUDED_core_scoring_packstat_types_hh

#include <core/scoring/packstat/types.fwd.hh>

#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <utility/VirtualBase.hh>


#include <utility/vector1.hh> // AUTO IWYU For vector1


namespace core {
namespace scoring {
namespace packstat {

typedef core::Real PackstatReal;

typedef utility::vector1<PackstatReal> Floats;
typedef Floats::iterator FloatIter;
typedef Floats::const_iterator FloatCIter;

typedef utility::vector1<core::Real> Reals;
typedef Reals::iterator RealIter;
typedef Reals::const_iterator RealCIter;

typedef numeric::xyzVector<PackstatReal> XYZ;
typedef utility::vector1<XYZ> XYZs;
typedef XYZs::iterator XYZIter;
typedef XYZs::const_iterator XYZCIter;


struct Sphere {

	Sphere( XYZ xyz_in, PackstatReal rad_in, core::Size _id = 0 ) :
		xyz(xyz_in),
		radius(rad_in),
		sasa(0),
		id(_id),
		aid(1,_id)
	{}

	Sphere( XYZ xyz_in, PackstatReal rad_in, core::id::AtomID _aid ) :
		xyz(xyz_in),
		radius(rad_in),
		sasa(0),
		id(0),
		aid(_aid)
	{}

	XYZ xyz;
	PackstatReal radius,sasa;
	core::Size id;
	core::id::AtomID aid;
};

typedef utility::vector1<Sphere> Spheres;
typedef Spheres::iterator SphereIter;
typedef Spheres::const_iterator SphereCIter;

struct PosePackData : public utility::VirtualBase {
	Spheres spheres;
	XYZs centers;
	utility::vector1<std::string> labels;
};


} // namespace packstat
} // namespace scoring
} // namespace core


#endif
