>// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   energy_methods/SSElementMotifContactEnergy.hh
/// @brief  Will's motif score to determine how well packed the protein core is
/// @author TJ Brunette


#ifndef INCLUDED_core_energy_methods_SSElementMotifEnergy_hh
#define INCLUDED_core_energy_methods_SSElementMotifEnergy_hh


// Package headers
#include <basic/datacache/CacheableData.hh>

#include <energy_methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>



// Utility headers


namespace core {
namespace energy_methods {


class CenPairMotifEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent;

public:
	CenPairMotifEnergy();

	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const {
		return utility::pointer::make_shared< CenPairMotifEnergy >();
	}

	/// @brief Called at the end of the energy evaluation.
	virtual void finalize_total_energy( pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const;


	virtual void indicate_required_context_graphs( utility::vector1< bool > & ) const {};


	virtual core::Size version() const;

private:
	core::scoring::motif::MotifHashManager *mman_;
	utility::vector1<Size> aalist_;
};

} // scoring
} // core

#endif
