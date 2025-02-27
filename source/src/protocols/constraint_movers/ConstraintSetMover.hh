// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Assigns a ConstraintSet to a pose. Reads and creats ConstraintSet from file via command line option -constraints::cst_file,
/// unless a ConstraintSet is supplied via the constructor or the constraint_set() method.
/// @author ashworth
/// @author Setters/getters added by Vikram K. Mulligan (vmulligan@flatironinstitute.org) to allow configuration
/// from C++ or Python code.

#ifndef INCLUDED_protocols_constraint_movers_ConstraintSetMover_hh
#define INCLUDED_protocols_constraint_movers_ConstraintSetMover_hh

#include <protocols/constraint_movers/ConstraintSetMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>



#include <unordered_map>


namespace protocols {
namespace constraint_movers {

class ConstraintSetMover : public protocols::moves::Mover {

public:
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSetCOP ConstraintSetCOP;

public:
	ConstraintSetMover();
	~ConstraintSetMover() override;

	/// @brief Options constructor that sets the mover name.
	/// @note This does NOT set the consraint file!  To do that, call
	/// ConstraintSetMover::constraint_file() after construction.
	ConstraintSetMover( std::string const & );

	void read_options();
	void constraint_file( std::string const & );

	void constraint_set( ConstraintSetCOP );
	ConstraintSetOP constraint_set();
	ConstraintSetCOP constraint_set() const;

	void apply( Pose & ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Associates relevant options with the ConstraintSetMover class
	static void register_options();

	void
	parse_my_tag( TagCOP, basic::datacache::DataMap & ) override;

	void add_constraints( bool const a ){ add_constraints_ = a; }
	bool add_constraints() const { return add_constraints_; }

	/// @brief Set the (centroid) constraint file name.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void set_cst_file( std::string const & file_name );

	/// @brief Set the (full-atom) constraint file name.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void set_cst_fa_file( std::string const & file_name );

	/// @brief Set the constraint map file name.
	void set_cst_map_file( std::string const & file_name );

	/// @brief Set the constraint map file from previously loaded file contents in memory.
	void set_cst_map_file_contents( std::string const & file_contents );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	ConstraintSetOP constraint_set_low_res_;
	ConstraintSetOP constraint_set_high_res_;
	std::string cst_file_;
	std::string cst_fa_file_;
	bool add_constraints_; // dflt false; if true add the constraints, rather than replacing

	std::unordered_map< std::string, std::string > file_to_input_map_;
	std::string current_input_name_; // from the DataMap

};

} // moves
} // protocols

#endif
