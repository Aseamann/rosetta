// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Basic includes
# include <iostream>
#include <numeric>

// Rosetta includes
# include <basic/Tracer.hh>  // Output
# include <basic/options/option.hh>
# include <utility/vector1.hh>
# include <basic/options/keys/in.OptionKeys.gen.hh>
# include <devel/init.hh>
# include <core/import_pose/import_pose.hh>  // importing pdbs
# include <utility/pointer/owning_ptr.hh>  // Needed for every time a smart pointer is used
# include <core/pose/Pose.hh>  // Needed for every time a smart pointer is used
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
# include <core/scoring/ScoreTypes.hh>  // for setting weight

// Setup tracer cout
static basic::Tracer TR( "apps.pilot.aseamann" );

int main( int argc, char ** argv ) {
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
    if ( filenames.size() > 0 ) {
    TR << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    } else {
        TR << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
        return 1;
    }

    // Read in pdb
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
    // Collect length of residues
    int N = mypose->size();

    // Setup fold tree from secondary structure
    core::kinematics::FoldTreeOP ftree = protocols::bootcamp::fold_tree_from_ss( *mypose );
    mypose->fold_tree( *ftree );

    // Set score funciton
    core::Size linear_chainbreak = 1;
    // Add cutpoint residues to pose
    mypose->add_variant_type_to_pose_residue( *mypose );
    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
    // Set weight for linear chainbreak
    core::scoring::ScoreType linear_chainbreak_score = core::scoring::linear_chainbreak;
    sfxn.set_weight(linear_chainbreak_score, 1)
    core::Real score = sfxn->score( *mypose );

    // Print score
    // mypose->energies().show( std::cout );
    std::cout << "Score Results is: " << score << std::endl;

    // Setup monte carlo object
    core::Real temp = 1.0;
    protocols::moves::MonteCarloOP mc ( new protocols::moves::MonteCarlo( *mypose, *sfxn, temp ) );

    // Setup pymol observer
    protocols::moves::PyMOLObserverOP the_observer = protocols::moves::AddPyMOLObserver( *mypose, true, 0 );

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

    int loop_count = 100;
    for( int i = 0 ; i < loop_count; i++ ) {
        // Determine random residue
        double uniform_random_number = numeric::random::uniform();
        core::Size randres = static_cast< core::Size > (uniform_random_number * N + 1 );
        // Collect random numbers to adjust phi/psi
        core::Real pert1 = numeric::random::gaussian();
        core::Real pert2 = numeric::random::gaussian();
        // Collect original phi/psi
        core::Real orig_phi = mypose->phi( randres );
        core::Real orig_psi = mypose->psi( randres );
        // Set random phi/psi
        mypose->set_phi( randres, orig_phi + pert1 );
        mypose->set_psi( randres, orig_psi + pert2 );

        // Setup PackerTask to run after move phi/psi
        core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
        repack_task->restrict_to_repacking();
        core::pack::pack_rotamers( *mypose, *sfxn, repack_task );

        // Run minimizer
        copy_pose = *mypose;
        atm.run( copy_pose, mm, *sfxn, min_opts );
        *mypose = copy_pose;

        // apply monte carlo
        bool result_mc = mc->boltzmann( *mypose );

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
        the_observer->pymol().apply( *mypose );
    }

    // Check pose
    core::Real score2 = sfxn->score( *mypose );
    TR << "Updated Score: " << score2 << std::endl;

    return 0;
}