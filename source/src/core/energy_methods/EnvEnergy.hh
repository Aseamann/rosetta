// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/EnvEnergy.hh
/// @brief  EnergyMethod that evaluates env and cbeta energies.
/// @author Phil Bradley


#ifndef INCLUDED_core_energy_methods_EnvEnergy_hh
#define INCLUDED_core_energy_methods_EnvEnergy_hh

// Unit Headers

// Package headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/EnvPairPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace energy_methods {



class EnvEnergy : public core::scoring::methods::ContextDependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextDependentOneBodyEnergy  parent;
public:


	EnvEnergy();


	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;


	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// totals
	) const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}


	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	// const-ref to scoring database
	core::scoring::EnvPairPotential const & potential_;
	core::Size version() const override;

};


}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
