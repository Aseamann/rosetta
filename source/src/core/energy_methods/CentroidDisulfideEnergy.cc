// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/CentroidDisulfideEnergy.cc
/// @brief  Centroid Disulfide Energy class implementation
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   2/4/09


// Unit headers
#include <core/energy_methods/CentroidDisulfideEnergy.hh>
#include <core/energy_methods/CentroidDisulfideEnergyCreator.hh>

// Package headers
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>
#include <core/scoring/disulfides/CentroidDisulfideEnergyContainer.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/Methods.hh>
#include <basic/Tracer.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the CentroidDisulfideEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
CentroidDisulfideEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< CentroidDisulfideEnergy >(
		core::scoring::ScoringManager::get_instance()->get_CentroidDisulfidePotential()
	);
}

core::scoring::ScoreTypes
CentroidDisulfideEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( dslfc_cen_dst );
	sts.push_back( dslfc_cb_dst );
	sts.push_back( dslfc_ang );
	sts.push_back( dslfc_cb_dih );
	sts.push_back( dslfc_bb_dih );
	return sts;
}


static basic::Tracer TR( "core.energy_methods.CentroidDisulfideEnergy" );

CentroidDisulfideEnergy::CentroidDisulfideEnergy(
	core::scoring::disulfides::CentroidDisulfidePotential const & potential
) :
	parent( utility::pointer::make_shared< CentroidDisulfideEnergyCreator >() ),
	potential_( potential )
{}

CentroidDisulfideEnergy::~CentroidDisulfideEnergy() = default;

// EnergyMethod Methods:

core::scoring::methods::EnergyMethodOP
CentroidDisulfideEnergy::clone() const
{
	return utility::pointer::make_shared< CentroidDisulfideEnergy >( potential_ );
}


void
CentroidDisulfideEnergy::setup_for_scoring(
	pose::Pose & pose,
	core::scoring::ScoreFunction const & ) const
{
	using namespace core::scoring::methods;
	using namespace core::scoring::disulfides;

	if ( pose.energies().long_range_container( centroid_disulfide_energy ) == nullptr ) {
		CentroidDisulfideEnergyContainerOP dec( new CentroidDisulfideEnergyContainer( pose ) );
		pose.energies().set_long_range_container( centroid_disulfide_energy, dec );
	} else {
		CentroidDisulfideEnergyContainerOP dec = CentroidDisulfideEnergyContainerOP (
			utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer > ( pose.energies().nonconst_long_range_container( centroid_disulfide_energy ) ));
		dec->update( pose );
	}
}

void
CentroidDisulfideEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & ) const
{}


// TwoBodyEnergy Methods:

void
CentroidDisulfideEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd1.has_variant_type( core::chemical::REPLONLY ) || rsd2.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}


	Energy cbcb_distance_score;
	Energy centroid_distance_score;
	Energy cacbcb_angle_1_score;
	Energy cacbcb_angle_2_score;
	Energy cacbcbca_dihedral_score;
	Energy backbone_dihedral_score;

	//Require cysteines
	if ( rsd1.aa() != chemical::aa_cys || rsd2.aa() != chemical::aa_cys ) return;
	//Require Centroid
	if ( rsd1.type().mode() != chemical::CENTROID_t ||
			rsd2.type().mode() != chemical::CENTROID_t ) {
		return;
	}

	core::scoring::disulfides::CentroidDisulfideEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer const > (
		pose.energies().long_range_container( core::scoring::methods::centroid_disulfide_energy )
		)
	);
	//Require they're bonded
	if ( ! dec->residue_forms_disulfide( rsd1.seqpos() ) ||
			dec->other_neighbor_id( rsd1.seqpos() ) != (Size) rsd2.seqpos() ) {
		return;
	}

	potential_.score_disulfide(
		rsd1, rsd2,
		cbcb_distance_score,
		centroid_distance_score,
		cacbcb_angle_1_score,
		cacbcb_angle_2_score,
		cacbcbca_dihedral_score,
		backbone_dihedral_score
	);

	emap[ core::scoring::dslfc_cen_dst ] += centroid_distance_score;
	emap[ core::scoring::dslfc_cb_dst ]  += cbcb_distance_score;
	emap[ core::scoring::dslfc_ang ]     += (cacbcb_angle_1_score + cacbcb_angle_2_score)*.5;
	emap[ core::scoring::dslfc_cb_dih ]  += cacbcbca_dihedral_score;
	emap[ core::scoring::dslfc_bb_dih ]  += backbone_dihedral_score;
}


bool CentroidDisulfideEnergy::defines_intrares_energy( core::scoring::EnergyMap const & ) const
{
	return false;
}


void CentroidDisulfideEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{}

// LongRangeTwoBodyEnergy methods
core::scoring::methods::LongRangeEnergyType
CentroidDisulfideEnergy::long_range_type() const
{
	return core::scoring::methods::centroid_disulfide_energy;
}


bool
CentroidDisulfideEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const
{
	using namespace core::scoring::methods;
	if ( ! pose.energies().long_range_container( centroid_disulfide_energy ) ) return false;

#ifdef NDEBUG
	core::scoring::disulfides::CentroidDisulfideEnergyContainerCOP dec(
		utility::pointer::static_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer const > (
		pose.energies().long_range_container( centroid_disulfide_energy)
		)
	);
#else
	core::scoring::disulfides::CentroidDisulfideEnergyContainerCOP dec(
		utility::pointer::dynamic_pointer_cast< core::scoring::disulfides::CentroidDisulfideEnergyContainer const > (
		pose.energies().long_range_container( centroid_disulfide_energy )
		)
	);
	debug_assert( dec != nullptr );
#endif
	return dec->disulfide_bonded( res1, res2 );
}
core::Size
CentroidDisulfideEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace energy_methods
} // namespace core

