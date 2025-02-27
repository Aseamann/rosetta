// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/EnzymaticMover.fwd.hh
/// @brief  Declarations and simple accessor/mutator definitions for the base class EnzymaticMover
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_EnzymaticMover_HH
#define INCLUDED_protocols_enzymatic_movers_EnzymaticMover_HH

// Unit headers
#include <protocols/enzymatic_movers/EnzymaticMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <utility/vector1.hh> // AUTO IWYU For vector1


namespace protocols {
namespace enzymatic_movers {

typedef std::pair< core::uint, std::string > ReactionSite;

/// @details  This is a base class for any Mover that modifies the sequence of a Pose in a way that simulates the
/// action of an enzyme.
/// Any Mover inheriting from this base class must provide an enzyme family corresponding to a directory of enzyme data
/// and must implement the perform_reaction() method, which modifies, adds, or removes (a) Residue(s).
/// The core machinery of this base class uses the enzymatic data to search for potential reaction sites.
class EnzymaticMover : public moves::Mover {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	EnzymaticMover();

	/// @brief  Constructor with enzyme family provided
	EnzymaticMover( std::string const & enzyme_family );

	/// @brief  Copy constructor
	EnzymaticMover( EnzymaticMover const & ) = default;

	// Destructor
	~EnzymaticMover() override = default;


public: // Standard Rosetta methods ///////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Generate string representation of EnzymaticMover for debugging purposes.
	void show( std::ostream & output=std::cout ) const override;


	// Mover methods
	/// @brief  Return the name of the Mover.
	std::string get_name() const override = 0;  // This method must be overridden.

	moves::MoverOP clone() const override = 0;

	moves::MoverOP fresh_instance() const override = 0;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	static utility::tag::XMLSchemaComplexTypeGeneratorOP xml_schema_complex_type_generator();


	/// @brief  Apply the corresponding move to <input_pose>.
	void apply( core::pose::Pose & input_pose ) override;


public: // Citation Management ////////////////////////////////////////////////

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

protected: // Accessors/Mutators //////////////////////////////////////////////
	/// @brief  Set the family name of this simulated enzyme.
	void set_enzyme_family( std::string const & family_name );

	/// @brief   Actually perform the virtual reaction for this specific EnzymaticMover.
	/// @params  <input_pose>:  The Pose to be acted on.
	/// @params  <site>:        An index for the reactive site to be acted on.
	/// @params  <cosubstrate>: A string providing information to the Enzymatic Mover to add or subtract the proper
	/// atoms or Residues from the Pose.  How this parameter is handled is specific to each distinct EnzymaticMover.
	virtual void perform_reaction(
		core::pose::Pose & input_pose,
		core::uint const site,
		std::string const & cosubstrate="" ) = 0;  // This method must be overridden.

public:
	/// @brief    Get the family name of this simulated enzyme.
	/// @details  This EnzymaticMover is a member of this enyzme family.
	std::string
	get_enzyme_family() const
	{
		return enzyme_family_;
	}

	/// @brief  Set the species name of this simulated enzyme.
	void set_species( std::string const & species_name );

	/// @brief    Get the species name of this simulated enzyme.
	/// @details  This EnzymaticMover is limited to reactions known to occur in the returned species.
	std::string
	get_species() const
	{
		return species_name_;
	}


	/// @brief  Set the specific name of this simulated enzyme.
	void set_enzyme( std::string const & enzyme_name );

	/// @brief  Return the specific name of this simulated enzyme.
	std::string
	get_enzyme() const
	{
		return enzyme_name_;
	}


	/// @brief    Directly set the efficiency of this enzyme, ignoring whatever is in the database.
	/// @details  An efficiency of 0.45 means that this enzyme will perform its reaction at any recognized site 45% of
	/// the time.
	void
	set_efficiency( core::Real setting )
	{
		efficiency_ = setting;
	}

	/// @brief    Get the efficiency of this enzyme.
	/// @details  An efficiency of 0.45 means that this enzyme will perform its reaction at any recognized site 45% of
	/// the time.
	core::Real
	get_efficiency() const
	{
		return efficiency_;
	}


	/// @brief  Do not perform a reaction at this site, even if it is within a consensus sequence match for this enzyme.
	void exclude_site( core::uint seqpos );

	/// @brief  Pass a list of sequence positions that are forbidden from being modified, even if they are within a
	/// consensus sequence match for this enzyme.
	void
	set_excluded_sites( utility::vector1< core::uint > const & excluded_sites )
	{
		excluded_sites_ = excluded_sites;
	}

	/// @brief  Return a list of sequence positions that are forbidden from being modified, even if they are within a
	/// consensus sequence match for this enzyme.
	utility::vector1< core::uint >
	get_excluded_sites() const
	{
		return excluded_sites_;
	}


	/// @brief  Definitely modify this site, if it is within a consensus sequence match for this enzyme.
	void ensure_site( core::uint seqpos );

	/// @brief  Pass a list of sequence positions that are guaranteed to be modified, if they are within a consensus
	/// sequence match for this enzyme.
	void
	set_ensured_sites( utility::vector1< core::uint > const & ensured_sites )
	{
		ensured_sites_ = ensured_sites;
	}

	/// @brief  Return a list of sequence positions that are guaranteed to be modified, if they are within a consensus
	/// sequence match for this enzyme.
	utility::vector1< core::uint >
	get_ensured_sites() const
	{
		return ensured_sites_;
	}


	/// @brief  Return the current number of reactive sites.
	core::Size
	get_n_reactive_sites() const
	{
		return reaction_sites_.size();
	}

	/// @brief  Return the sequence position of the requested reactive site.
	core::uint
	get_reactive_site_sequence_position( core::uint const index ) const
	{
		return reaction_sites_[ index ].first;
	}

	/// @brief  Return the atom name of the requested reactive site.
	std::string const &
	get_reactive_site_atom_name( core::uint const index ) const
	{
		return reaction_sites_[ index ].second;
	}


	/// @brief  Return the current number of reactive sites.
	core::Size
	get_n_co_substrates() const
	{
		return co_substrates_.size();
	}

	/// @brief  Return the requested cosubstrate of this enzymatic reaction.
	std::string const &
	get_co_substrate( core::uint const index ) const
	{
		return co_substrates_[ index ];
	}


	/// @brief  Set this EnzymaticMover to perform only its major reaction.
	void
	perform_major_reaction_only()
	{
		performs_major_reaction_only_ = true;
	}

	/// @brief  Allow this EnzymaticMover to be promiscuous, performing a random transfer from among its possible
	/// co-substrates.
	void
	perform_all_reactions()
	{
		performs_major_reaction_only_ = false;
	}

	/// @brief  Does this enzyme only perform its major reaction?
	bool
	performs_major_reaction_only() const
	{
		return performs_major_reaction_only_;
	}


public:  // Other Methods /////////////////////////////////////////////////////
	// Access the EnzymeManager to determine which sites on a given Pose are able to be modified.
	void set_pose_reactive_sites( core::pose::Pose const & pose );


private:  // Private methods //////////////////////////////////////////////////
	// Set command-line options.  (Called by init().)
	void set_commandline_options();

	// Initialize data members from arguments.
	void init( std::string const & enzyme_family );


	// Access the EnzymeManager to set efficiency.
	void set_efficiency();

	// Access the EnzymeManager to determine which co-substrates are able to take part in the reaction of this
	// particular enzyme.
	void set_available_co_substrates();


private:  // Private data /////////////////////////////////////////////////////
	std::string enzyme_family_ = "";  // e.g., "glycosyltransferases"
	std::string species_name_ = "";  // e.g., "e_coli" or "h_sapiens"
	std::string enzyme_name_ = "";  // for if we want to give this a real-life enzyme name
	core::Real efficiency_ = 1.0;  // ratio of times this enzyme performs its reaction at a recognized site
	utility::vector1< core::uint > excluded_sites_;  // user-specified sequence positions to be excluded from activity
	utility::vector1< core::uint > ensured_sites_;  // user-specified sequence positions to guarantee activity
	utility::vector1< ReactionSite > reaction_sites_;  // sequence positions and atoms to be reacted
	utility::vector1< std::string > co_substrates_;  // descriptors of co-substrates

	// Options
	bool performs_major_reaction_only_ = false;  // Only the most common reaction will occur.
};  // class EnzymaticMover

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_EnzymaticMover_HH
