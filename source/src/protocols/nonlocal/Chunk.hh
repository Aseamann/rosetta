// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/Chunk.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_CHUNK_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_CHUNK_HH

// Unit headers
#include <protocols/nonlocal/Chunk.fwd.hh>

// External headers
#include <boost/math/distributions/normal.hpp>

// Package headers
#include <protocols/nonlocal/Region.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <numeric/random/DistributionSampler.fwd.hh>
#include <utility/VirtualBase.hh>

#include <memory>

namespace protocols {
namespace nonlocal {

/// @class An object that contains the start/stop points of a continuous
/// region of sequence and knowledge of how to sample insertion positions
/// from it according to an underlying distribution. Certain residues within
/// the region may not be movable. Current behavior mimics but improves upon
/// existing end-biased selection.
class Chunk : public utility::VirtualBase {
	typedef boost::math::normal Normal;
	typedef core::Size Size;
	typedef core::kinematics::MoveMapOP MoveMapOP;

public:
	// -- Construction and Assignment -- //

	/// @brief Default constructor
	/// <region> describes the contiguous stretch of residues.
	/// <movable> describes the modifiable degrees of freedom in the system.
	Chunk( RegionOP region, MoveMapOP movable);

	/// @brief Copy constructor
	Chunk(const Chunk& other);

	/// @brief Destructor
	~Chunk() override ; // auto-removing definition from header{}

	/// @brief Assignment operator
	Chunk& operator=(const Chunk& other);

	// -- Accessors -- //

	/// @brief Chooses an allowable insertion position on [start, stop] according
	/// to the probability distribution
	core::Size choose() const;

	/// @brief Lower boundary of this chunk
	core::Size start() const;

	/// @brief Upper boundary of this chunk
	core::Size stop() const;

	/// @brief Returns the length of this region
	core::Size length() const;

	/// @brief Returns true if at least one position on [start(), stop()] is movable
	bool is_movable() const;

	/// @brief Returns true if there is at least one valid insertion position in
	/// the closed region [start(), stop()], false otherwise.
	bool valid() const;

private:
	RegionOP region_;
	MoveMapOP movable_;
	std::unique_ptr<numeric::random::DistributionSampler<Normal> > sampler_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_CHUNK_HH_
