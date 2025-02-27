// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/ClearCompositionConstraintsMover.cc
/// @brief Clears all AACompositionConstraints from a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/aa_composition/ClearCompositionConstraintsMover.hh>
#include <protocols/aa_composition/ClearCompositionConstraintsMoverCreator.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/datacache/DataMap.fwd.hh>

//Auto Headers

#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.aa_composition.ClearCompositionConstraintsMover" );

namespace protocols {
namespace aa_composition {

using namespace core;
using namespace core::scoring;
using namespace constraints;
using namespace utility::tag;




/// @brief Default Constructor
///
ClearCompositionConstraintsMover::ClearCompositionConstraintsMover():
	protocols::moves::Mover( ClearCompositionConstraintsMover::mover_name() )
	//TODO initialize variables here
{
}

/// @brief Copy Constructor
///
ClearCompositionConstraintsMover::ClearCompositionConstraintsMover( ClearCompositionConstraintsMover const &/*src*/ ):
	protocols::moves::Mover( ClearCompositionConstraintsMover::mover_name() )
	//TODO initialize variables here
{
}


/// @brief Destructor
///
ClearCompositionConstraintsMover::~ClearCompositionConstraintsMover()= default;

/// @brief Copy this object and return a pointer to the copy.
///
protocols::moves::MoverOP ClearCompositionConstraintsMover::clone() const { return utility::pointer::make_shared< protocols::aa_composition::ClearCompositionConstraintsMover >( *this ); }

/// @brief Create a new object of this type and return a pointer to it.
///
protocols::moves::MoverOP ClearCompositionConstraintsMover::fresh_instance() const { return utility::pointer::make_shared< ClearCompositionConstraintsMover >(); }

/// @brief Returns the name of this mover ("ClearCompositionConstraintsMover").
///

/// @brief Actually apply the mover to a pose.
///
void
ClearCompositionConstraintsMover::apply( Pose &pose )
{
	pose.clear_sequence_constraints();
	if ( TR.visible() ) {
		TR << "Removed all sequence constraints from the pose." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Parse RosettaScripts XML tag to set up the mover.
///
void
ClearCompositionConstraintsMover::parse_my_tag(
	TagCOP const /*tag*/,
	basic::datacache::DataMap &/*datamap*/
)
{
	//This mover takes no options.
	return;
}

std::string ClearCompositionConstraintsMover::get_name() const {
	return mover_name();
}

std::string ClearCompositionConstraintsMover::mover_name() {
	return "ClearCompositionConstraintsMover";
}

void ClearCompositionConstraintsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Clears all composition constraints from the pose.  Leaves geometric constraints intact.", attlist );
}

std::string ClearCompositionConstraintsMoverCreator::keyname() const {
	return ClearCompositionConstraintsMover::mover_name();
}

protocols::moves::MoverOP
ClearCompositionConstraintsMoverCreator::create_mover() const {
	return utility::pointer::make_shared< ClearCompositionConstraintsMover >();
}

void ClearCompositionConstraintsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ClearCompositionConstraintsMover::provide_xml_schema( xsd );
}


} // aa_composition
} // protocols
