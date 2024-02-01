// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief 2023_bootcamp_mover_subclass
/// @author Aseamann (aseamann@unomaha.edu)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// Basic includes
#include <iostream>
#include <numeric>

// Rosetta includes
# include <basic/options/option.hh>
# include <basic/options/keys/in.OptionKeys.gen.hh>
# include <devel/init.hh>
# include <core/import_pose/import_pose.hh>  // importing pdbs
# include <utility/pointer/owning_ptr.hh>  // Needed for every time a smart pointer is used
# include <core/scoring/ScoreFunctionFactory.hh>  // Needed for scoring
# include <core/scoring/ScoreFunction.hh> // Needed for scoring 
# include <core/scoring/Energies.hh> // Showing energies from scoring
# include <numeric/random/random.hh>  // For random number
# include <protocols/moves/MonteCarlo.hh>  // For montecarlo
# include <protocols/moves/PyMOLMover.hh>  // PyMolObserver
# include <core/pack/task/PackerTask.hh>  // for setting up PackerTask
# include <core/pack/task/TaskFactory.hh>  // for setting up TaskFactory
# include <core/pack/pack_rotamers.hh>  // for packing
# include <core/kinematics/MoveMap.hh> // for minimizing
# include <core/optimization/AtomTreeMinimizer.hh> // for minimizing
# include <core/optimization/MinimizerOptions.hh> // for minimizing
# include <protocols/bootcamp/fold_tree_from_ss.hh>  // for fold_tree_from_ss
# include <core/pose/variant_util.hh>  // setup cutpoints
# include <core/scoring/ScoreType.hh>  // for setting weight

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& mypose ){
	// Collect length of residues
    int N = mypose.size();

    // Setup fold tree from secondary structure
    core::kinematics::FoldTreeOP ftree = protocols::bootcamp::fold_tree_from_ss( mypose );
    mypose.fold_tree( *ftree );

    // Set score funciton
    // Add cutpoint residues to pose
    // add_variant_type_to_pose_residue( mypose );
	correctly_add_cutpoint_variants( mypose );
    // Set weight for linear chainbreak
    core::scoring::ScoreType linear_chainbreak_score = core::scoring::linear_chainbreak;
    sfxn_->set_weight(linear_chainbreak_score, 1);
    core::Real score = sfxn_->score( mypose );

    // Print score
    // mypose->energies().show( std::cout );
    std::cout << "Score Results is: " << score << std::endl;

    // Setup monte carlo object
    core::Real temp = 1.0;
    protocols::moves::MonteCarloOP mc ( new protocols::moves::MonteCarlo( mypose, *sfxn_, temp ) );

    // Setup pymol observer
    protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( mypose, true, 0 );

    // Setup minimizer
    core::pose::Pose copy_pose;
    core::kinematics::MoveMap mm;
    mm.set_bb( true );
    mm.set_chi( true );
    core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
    core::optimization::AtomTreeMinimizer atm;

    // Data collection for acceptance rate
    core::Size accepted_num = 0;
    utility::vector1< core::Real > score_vec{};

    int loop_count = num_iteractions_;
    for( int i = 0 ; i < loop_count; i++ ) {
        // Determine random residue
        double uniform_random_number = numeric::random::uniform();
        core::Size randres = static_cast< core::Size > (uniform_random_number * N + 1 );
        // Collect random numbers to adjust phi/psi
        core::Real pert1 = numeric::random::gaussian();
        core::Real pert2 = numeric::random::gaussian();
        // Collect original phi/psi
        core::Real orig_phi = mypose.phi( randres );
        core::Real orig_psi = mypose.psi( randres );
        // Set random phi/psi
        mypose.set_phi( randres, orig_phi + pert1 );
        mypose.set_psi( randres, orig_psi + pert2 );

        // Setup PackerTask to run after move phi/psi
        core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( mypose );
        repack_task->restrict_to_repacking();
        core::pack::pack_rotamers( mypose, *sfxn_, repack_task );

        // Run minimizer
        copy_pose = mypose;
        atm.run( copy_pose, mm, *sfxn_, min_opts );
        mypose = copy_pose;

        // apply monte carlo
        bool result_mc = mc->boltzmann( mypose );

        // Collect monte carlo data
        if ( result_mc ) {
            // Update count of accepted poses
            accepted_num++;
            // Save pose score update
            core::Real result = mc->last_accepted_score();
            score_vec.push_back( result );
        }

        // Print results if full
        if ( ( i % 99 == 0 ) && ( i != 0 ) ) {
            core::Real vec_size = score_vec.size();
            core::Real avg_score =  std::accumulate( score_vec.begin(), score_vec.end(), 0.0) / vec_size;
            TR << "Average of score last accapted poses: " << avg_score << std::endl;
            TR << "Out of the last 100 poses " << accepted_num << " were accepted" << std::endl;
            score_vec.clear();
            accepted_num = 0;
        }

        // Call pymol observer
        the_observer->pymol().apply( mypose );
    }

    // Check pose
    core::Real score2 = sfxn_->score( mypose );
    TR << "Updated Score: " << score2 << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Setters and Getters ///
///////////////////////////

/// @brief Set the number of iteractions
void
BootCampMover::set_num_iterations( core::Size num_iteractions ) {
	num_iteractions_ = num_iteractions;
}

/// @brief Set the score function
void
BootCampMover::set_score_function( core::scoring::ScoreFunctionOP sfxn ) {
	sfxn_ = sfxn;
}

/// @brief Get the number of iteractions
core::Size
BootCampMover::get_num_iterations() const {
	return num_iteractions_;
}

/// @brief Get the score function
core::scoring::ScoreFunctionOP
BootCampMover::get_score_function() const {
	return sfxn_;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "2023_bootcamp_mover_subclass", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Aseamann as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"Aseamann",
		"TODO: institution",
		"austin.seamann@rutgers.edu",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
