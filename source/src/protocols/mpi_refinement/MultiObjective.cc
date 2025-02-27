// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/mpi_refinement/MultiObjective.cc
/// @brief
/// @author Hahnbeom Park

#include <protocols/mpi_refinement/util.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> // AUTO IWYU For SwitchResidueTypeSetMover
#include <utility/vector0.hh> // AUTO IWYU For vector0

// Utility headers

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("MPI.LHR.O");

// Methods should be REGISTERED here in order to use as multi-objective
MultiObjective::MultiObjective()
{
	set_defaults();
}

MultiObjective::~MultiObjective()= default;

void
MultiObjective::set_defaults()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	objsfxnOPs_.resize( 0 );
	fobjnames_.resize( 0 );

	fobjnames_.push_back( "score" );

	// Add extra multi objective functions
	if ( !option[ OptionKeys::lh::multi_objective_functions ].user() ) {
		TR << "Add default extra objective function." << std::endl;
		fobjnames_.push_back( "goap" );
		fobjnames_.push_back( "similarity" );

	} else {
		TR << "Setup multi-objective function defined by user." << std::endl;
		utility::vector1< std::string > fobjnames_add
			= option[ OptionKeys::lh::multi_objective_functions ]();

		for ( core::Size i_obj = 1; i_obj <= fobjnames_add.size(); ++i_obj ) {
			fobjnames_.push_back( fobjnames_add[i_obj] );
		}
	}

	// use only if specified by user
	//if( option[ constraints::cst_file ].user() ) fobjnames_.push_back( "cst_cen" );
	//if( option[ constraints::cst_fa_file ].user() ) fobjnames_.push_back( "cst_fa" );

	for ( core::Size i_obj = 1; i_obj <= fobjnames_.size(); ++i_obj ) {
		std::string const score_name( fobjnames_[i_obj] );

		if ( score_name == "goap" ) {
			TR << "- " << i_obj << ". Added Goap potential as 'goap'" << std::endl;
			core::scoring::ScoreFunctionOP sfxn =
				core::scoring::ScoreFunctionFactory::create_score_function( "empty" );
			sfxn->set_weight( core::scoring::goap, 0.01 );
			objsfxnOPs_.push_back( sfxn );

		} else if ( score_name == "score" ) {
			TR << "- " << i_obj << ". Added input score as 'score'" << std::endl;
			core::scoring::ScoreFunctionOP sfxn =
				core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );
			// modify cart_bonded term to be less sensitive to small non-ideality
			// note that this is just for "evaluation", not for sampling
			//sfxn->set_weight( core::scoring::cart_bonded, 0.1 );
			objsfxnOPs_.push_back( sfxn );
		} else if ( score_name == "similarity" ) {
			TR << "- " << i_obj << ". Added similarity as 'similarity'" << std::endl;
			// Just to ensure the length, but never been used
			core::scoring::ScoreFunctionOP sfxn =
				core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );
			objsfxnOPs_.push_back( sfxn );

		} else if ( score_name == "cst_fa" ||
				score_name == "cst_cen" ) {
			TR << "- " << i_obj << ". Added " << score_name << " as '" << score_name << "'" << std::endl;
			core::scoring::ScoreFunctionOP sfxn_cst( new core::scoring::ScoreFunction );
			sfxn_cst->set_weight( core::scoring::atom_pair_constraint, 1.0 );
			sfxn_cst->set_weight( core::scoring::coordinate_constraint, 1.0 );
			objsfxnOPs_.push_back( sfxn_cst );

		} else {
			TR.Warning << "Skip unrecognized multi-objective function: " << score_name << std::endl;
		}
	}

	// dominant cut initialize
	utility::vector1< core::Real > objective_cut_from_cmd =
		option[ lh::objective_dominate_cut ]();

	utility::vector1< core::Real > objective_inc_from_cmd =
		option[ lh::objective_cut_increment ]();

	if ( nobjs() != objective_cut_from_cmd.size() ) {
		utility_exit_with_message( "lh::objective_dominate_cut num. args are different from actual objective function definition!");
	}
	if ( nobjs() != objective_inc_from_cmd.size() ) {
		utility_exit_with_message( "lh::objective_dominate_cut num. args are different from actual objective function definition!");
	}

	obj_dominant_cut_ = objective_cut_from_cmd;
	obj_cut_increment_ = objective_inc_from_cmd;
	iha_cut_ = 25.0;
	iha_penalty_slope_ = 10.0;
	iha_penalty_mode_ = "relative";
	nremain_reset_ = 3;
}

// Check if ss1 is better than ss2 by any obj function
bool
MultiObjective::is_dominant( core::io::silent::SilentStructCOP ss1,
	core::io::silent::SilentStructCOP ss2 )
{
	utility::vector1< std::string > fobj_loc( fobjnames_ );

	for ( core::Size iobj = 1; iobj <= fobj_loc.size(); ++iobj ) {
		std::string const score_name( fobj_loc[iobj] );
		core::Real dobj = ss1->get_energy( score_name ) - ss2->get_energy( score_name );

		// when ss1 no more dominates ss2
		if ( dobj > obj_dominant_cut_[iobj] ) return false;
	}

	return true;
}

// Takes care of objective function "similarity"
void
MultiObjective::calculate_pool_diversity(
	protocols::wum::SilentStructStore &structs ) const
{
	core::Size const n( structs.size() );

	for ( core::Size i = 0; i < n; ++i ) {
		core::io::silent::SilentStructOP ss1 = structs.get_struct( i );
		calculate_structure_diversity( ss1, structs );
	}
}

void
MultiObjective::calculate_pool_diversity(
	protocols::wum::SilentStructStore &structs1,
	protocols::wum::SilentStructStore &structs2
) const
{
	core::Size const n( structs1.size() );

	for ( core::Size i = 0; i < n; ++i ) {
		core::io::silent::SilentStructOP ss1 = structs1.get_struct( i );
		calculate_structure_diversity( ss1, structs2 );
	}
}

void
MultiObjective::calculate_structure_diversity(
	core::io::silent::SilentStructOP ss1,
	protocols::wum::SilentStructStore &structs ) const
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//core::Real simlimit = option[ lh::similarity_reference ]();
	//core::Real simlimit = option[ lh::rms_limit ]();
	core::Real simlimit = option[ cm::similarity_limit ]();
	std::string similarity_method = option[ lh::similarity_method ]();
	std::string similarity_measure = option[ lh::similarity_measure ]();
	core::Size const n( structs.size() );
	core::Real const simtol = option[ lh::similarity_tolerance ]();

	//if measured by similarity to Emin
	utility::vector0< core::Real > dvector( n, 0.0 );
	core::Size nrel( dvector.size() );
	for ( core::Size j = 0; j < n; ++j ) {
		core::io::silent::SilentStructOP ss2 = structs.get_struct( j );
		std::string ssname = ss2->decoy_tag();

		// avoid itself
		if ( ssname == ss1->decoy_tag()  ) {
			nrel--;
			dvector[j] = 0.0;
		} else {
			core::Real dist( 0.0 );
			if ( similarity_measure == "Sscore"  ) {
				core::Real dumm; // Dummy variable for return-by-ref rmsd
				dist = CA_Sscore( ss1, ss2, dumm, true, 2.0 );
				//dist = CA_Sscore( ss1, ss2, dumm, 1.0 );

			} else if ( similarity_measure == "rmsd"  ) {
				CA_Sscore( ss1, ss2, dist, true, 2.0 );
				//CA_Sscore( ss1, ss2, dist, 1.0 );

			} else if ( similarity_measure == "looprmsd"  ) {
				std::string loopstr = option[ lh::loop_string ]();
				utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );
				CA_Sscore( ss1, ss2, dist, loopres, true, 2.0 );
				//CA_Sscore( ss1, ss2, dist, loopres, 1.0 );
			} else { }
			dvector[j] = dist;
			//TR << "d: " << j << " " << dist << std::endl;
		}
	}

	core::Real similarity( 0.0 );
	// Diversity based on similarity sum
	if ( similarity_method == "sum"  ) {
		for ( core::Size j = 0; j < n; ++j ) {
			similarity += dvector[j];
		}
		similarity *= 100.0/nrel;

	} else if ( similarity_method == "sigsum"  ) {
		for ( core::Size j = 0; j < n; ++j ) {
			if ( similarity_measure == "Sscore"  ) {
				if ( dvector[j] >= simlimit ) {
					similarity += 1.0;
				} else if ( dvector[j] > simlimit-simtol ) {
					similarity += 0.5 + 0.5*(dvector[j]-simlimit)/simtol;
				}
			} else {
				if ( dvector[j] <= simlimit ) {
					similarity += 1.0;
				} else if ( dvector[j] <= simlimit+simtol ) {
					similarity += 0.5 - 0.5*(dvector[j]-simlimit)/simtol;
				}
			}
		}
		similarity *= 100.0/nrel;

		// Diversity based on closest one
	} else if ( similarity_method == "min"  ||
			similarity_method == "frontiermin"  ) {
		for ( core::Size j = 0; j < n; ++j ) {
			if ( dvector[j] != 1.0 && dvector[j] > similarity ) similarity = dvector[j]*100.0;
		}
	}

	ss1->add_energy( "similarity", similarity );

}

// update library following Conformational Space Annealing (CSA) logic using seed info & distance cut
bool
MultiObjective::update_library_seeds(protocols::wum::SilentStructStore &structs,
	protocols::wum::SilentStructStore &new_structs,
	core::Real const dcut,
	utility::vector1< core::Size > const seeds, // should be 'poolid'
	std::string const prefix,
	std::string const objname,
	core::Size const maxreplace
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;

	//core::Real const simlimit = option[ OptionKeys::lh::rms_limit ]();
	//std::string const measure = option[ lh::similarity_measure ]();
	// let's use our own... rms_limit is too misleading
	core::Real const simlimit = option[ cm::similarity_limit ]();
	std::string const measure( "Sscore" );
	bool const reverse_order( false );

	// CSA will use only the first fobjnames_ registered, whatever number of multiobjs are assigned
	//std::string objname = fobjnames_[1];

	// info that should be conserved
	utility::vector1< std::string > columns_copy;
	columns_copy.push_back( "poolid" ); columns_copy.push_back( "nuse" ); columns_copy.push_back( "score" ); columns_copy.push_back( objname );

	// rescore silents at the beginning
	core::pose::Pose pose_tmp;
	for ( core::Size iss = 0; iss < structs.size(); ++iss ) {
		core::io::silent::SilentStructOP ss = structs.get_struct( iss );
		ss->fill_pose( pose_tmp );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose_tmp );
		objsfxnOPs_[1]->score( pose_tmp );
		ss->energies_from_pose( pose_tmp );
	}
	for ( core::Size iss = 0; iss < new_structs.size(); ++iss ) {
		core::io::silent::SilentStructOP ss = new_structs.get_struct( iss );
		ss->fill_pose( pose_tmp );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose_tmp );
		objsfxnOPs_[1]->score( pose_tmp );
		ss->energies_from_pose( pose_tmp );
	}

	// store starting reference pool
	protocols::wum::SilentStructStore ref_structs( structs );

	// add iHA penalty if defined
	if ( iha_cut_ > 0.0 && objname == "penaltysum" ) {
		TR << "Penalty activated based on structural difference from reference structure." << std::endl;
		add_init_dev_penalty( new_structs, init_pose(), iha_penalty_mode(), iha_cut(), iha_penalty_slope() );
		add_init_dev_penalty( ref_structs, init_pose(), iha_penalty_mode(), iha_cut(), iha_penalty_slope() ); //in case when they don't have it yet
	}

	// superimpose to lowest energy at the begining only
	SilentStructOP ss0 = structs.get_struct( 0 );
	TR << "superimpose all before starting" << std::endl;
	superimpose_all( ss0, new_structs, columns_copy );
	superimpose_all( ss0, ref_structs, columns_copy );

	TR << "simlimit: " << simlimit << " " << measure << std::endl;

	// first start with removing redundancy
	core::Size nbefore( new_structs.size() );
	TR << "nbefore: " << nbefore << std::endl;
	filter_similar( new_structs, measure, simlimit, objname );

	TR << "New structures, initial/filtered_by_similiarity:" << nbefore << " " << new_structs.size() << std::endl;

	// sort by energy
	ref_structs.sort_by( objname );
	new_structs.sort_by( objname );
	core::Size nref( ref_structs.size() );
	//core::Size nnew( new_structs.size() );

	/* TODO
	core::Real const var_new = get_RMSF( new_structs );
	core::Real const var_ref = get_RMSF( ref_structs );
	TR << "Variation in energy, RefStructs: " << F(8,3,ref_structs.get_struct( 0 )->get_energy( objname )) << " to "
	<< F(8,3,ref_structs.get_struct( nref-1 )->get_energy( objname ))
	<< " / NewStructs: " << F(8,3,new_structs.get_struct( 0 )->get_energy( objname )) << " to "
	<< F(8,3,new_structs.get_struct( nnew-1 )->get_energy( objname )) << std::endl;

	TR << "Variation in struct(RMSF), RefStructs: " << var_ref << " / NewStrcuts: " << var_new << std::endl;
	*/

	// add nuse for seeds
	std::map< core::Size, bool > found_new;
	TR << "Available seeds: ";
	for ( core::Size iseed = 1; iseed <= seeds.size(); ++iseed ) TR << " " << seeds[iseed];
	TR << std::endl;

	for ( core::Size iref = 0; iref < ref_structs.size(); ++iref ) {
		SilentStructOP ss2 = ref_structs.get_struct( iref );
		core::Size const ipool = core::Size(ss2->get_energy( "poolid" ));
		if ( seeds.contains( ipool ) ) {
			core::Size nuse = ss2->get_energy( "nuse" ) + 1;
			ss2->add_energy( "nuse", nuse );
			found_new[ ipool ] = false;
		}
	}

	protocols::wum::SilentStructStore selected;
	core::Real emax( 0.0 );
	core::Size nreplace( 0 );

	for ( core::Size inew = 0; inew < new_structs.size(); ++inew ) {

		SilentStructOP ss1 = reverse_order? new_structs.get_struct( new_structs.size() - inew - 1 ) : new_structs.get_struct( inew );

		core::Real const score1 = ss1->get_energy( objname );
		core::Size const parent( ss1->get_energy( "parent" ) );
		//TR << ss1->decoy_tag() << ": parent " << parent << ", count: " << found_new.count( parent ) << std::endl;
		bool found_close( false );

		ref_structs.sort_by( objname ) ;
		emax = ref_structs.get_struct( nref-1 )->get_energy( objname );

		//TR << "current Emax: " << emax << std::endl;

		// search by reverse order
		core::Real distmin( 99.0 );
		core::Real emin_close( 1e10 );
		core::Size closest( 1 );
		for ( core::Size iref = ref_structs.size()-1; int(iref) >=0; --iref ) {
			SilentStructOP ss2 = ref_structs.get_struct( iref );
			core::Real const score2 = ss2->get_energy( objname );
			//core::Size const ipool = core::Size(ss2->get_energy( "poolid" ));

			core::Real dist = distance( ss1, ss2, measure, false );
			// Dscore = 1 - Sscore
			if ( measure == "Sscore" ) dist = 1.0 - dist;

			if ( dist < distmin ) {
				distmin = dist;
				closest = iref;
			}
			if ( dist < dcut ) found_close = true;
			if ( dist < dcut && score2 < emin_close ) emin_close = score2;

			/*
			if( score1 < score2 ){ // lower
			TR << "Replace_Sim: " << I(3,ipool) << " by " << ss1->decoy_tag()
			<< ", Eref/Emax/Enew/dist: "
			<< F(8,3,score2) << "/" << F(8,3,emax) << "/" << F(8,3,score1)
			<< "/" << F(8,3,dist) << std::endl;

			succeed_substitute_info( ss1, ss2, false );
			ref_structs.store()[iref] = ss1;
			break;
			} else {
			//TR << "Keep " << iref << " " << ss2->decoy_tag() << " " << score2 << std::endl;
			}
			*/
		}

		// Decision
		if ( found_close ) {
			if ( score1 < emin_close && nreplace < maxreplace ) { // lower than any close refstruct
				// replace closest
				SilentStructOP ss2 = ref_structs.get_struct( closest );
				core::Real const score2 = ss2->get_energy( objname );
				core::Size const ipool = core::Size(ss2->get_energy( "poolid" ));

				TR << "Replace_Sim: " << I(3,ipool) << " by " << A(30, ss1->decoy_tag())
					<< ", Eref/Emax/Enew/dist: "
					<< F(8,3,score2) << "/" << F(8,3,emax) << "/" << F(8,3,score1)
					<< "/" << F(8,3,distmin) << std::endl;

				succeed_substitute_info( ss1, ss2, false );
				ref_structs.store()[closest] = ss1;
				nreplace++;

			} else {
				TR << "Pass_Sim   :        "   << A(30, ss1->decoy_tag() )
					<< ", Eref/Emax/Enew/mind: "
					<< F(8,3,emin_close) << "/" << F(8,3,emax) << "/" << F(8,3,score1)
					<< "/" << F(8,3,distmin) << std::endl;
			}

		} else {
			if ( score1 < emax && nreplace < maxreplace ) { // found new->reset
				SilentStructOP ss2 = ref_structs.get_struct( nref-1 );
				core::Real const score2 = ss2->get_energy( objname );
				core::Size const ipool = core::Size(ss2->get_energy( "poolid" ));
				TR << "Replace_Max: " << I(3,ipool) << " by "  << A(30, ss1->decoy_tag() )
					<< ", Eref=Emax/Enew/mind:          " << F(8,3,score2)
					<< "/" << F(8,3,score1) << "/" << F(8,3,distmin) << std::endl;

				succeed_substitute_info( ss1, ss2, true ); //reset!
				nreplace++;

				// in case when the seed found a new one, keep nuse for the seed
				//found_new[ ipool ] = true;
				if ( found_new.count( parent ) == 0 ) {
					TR << "No such parent " << parent << " found among seeds, check it out!" << std::endl;
				} else {
					found_new[ parent ] = true;
				}

				// keep size
				ref_structs.store().pop_back();
				ref_structs.add( ss1 );
			} else {
				TR << "Pass_Max   :       "  << A(30,ss1->decoy_tag())
					<< ",     /Emax/Enew/mind: " << "         /" << F(8,3,emax) << "/" << F(8,3,score1)
					<< "/" << F(8,3,distmin) << std::endl;
			}
		}
	}

	// retag and replace into "structs"
	structs.clear();
	ref_structs.sort_by( objname );

	for ( core::Size i = 1; i <= seeds.size(); ++i ) {
		if ( found_new.at( seeds[i] ) ) {
			TR << "Seed " << seeds[i] << ", found new!" << std::endl;
		} else {
			TR << "Seed " << seeds[i] << ", nothing found new..." << std::endl;
		}
	}

	core::Size nused_structs( 0 );
	core::Size ntotal( ref_structs.size() );
	for ( core::Size i = 0; i < ref_structs.size(); ++i ) {
		SilentStructOP ss = ref_structs.get_struct( i );
		std::stringstream sstream( "" );
		sstream << prefix << "." << i;
		core::Size const poolid( ss->get_energy( "poolid" ) );

		// increase nuse if the pool was seed but
		if ( seeds.contains( poolid ) ) {
			if ( found_new.at( poolid ) ) {
				int nuse = std::max( 0, int(ss->get_energy( "nuse" )-1) );
				ss->add_energy( "nuse", nuse );
			}
		}
		if ( ss->get_energy( "nuse" ) > 0 ) nused_structs++;

		ss->set_decoy_tag( sstream.str() );
		structs.add( ss );
	}

	// reset condition
	if ( ( ntotal-nused_structs ) <= nremain_reset() ) {
		for ( core::Size iss = 0; iss < structs.size(); ++iss ) {
			core::io::silent::SilentStructOP ss = structs.get_struct( iss );
			ss->add_energy( "nuse", 0.0 );
		}
	}

	// report
	structs.sort_by( objname );
	TR << "I/Pool/FObj/Score/Nuse/GDTMM/dist_to_emin" << std::endl;

	SilentStructOP ss_emin = structs.get_struct( 0 );
	for ( core::Size i = 0; i < structs.size(); ++i ) {
		SilentStructOP ss = structs.get_struct( i );
		core::Real dist = distance( ss_emin, ss, measure, false );
		if ( measure == "Sscore" ) dist = 1.0 - dist;

		TR << I(3,i)
			<< " " << I(3,ss->get_energy( "poolid" ))
			<< " " << F(8,3,ss->get_energy( objname ))
			<< " " << F(8,3,ss->get_energy( "score" ))
			<< " " << I(4,ss->get_energy( "nuse" ))
			<< " " << F(8,3,ss->get_energy( "GDTMM_final" ))
			<< F(8,3,dist) << std::endl;
	}

	return true;
} // update_library_seeds

// Update library pool based on NSGAII rule, the multiobj GA
bool
MultiObjective::update_library_NSGAII(protocols::wum::SilentStructStore &structs,
	protocols::wum::SilentStructStore &new_structs,
	core::Size const nmax,
	bool const update_obj_cut
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;

	core::Real const simlimit = option[ OptionKeys::cm::similarity_limit ]();
	std::string const measure = option[ lh::similarity_measure ]();

	// Total pool
	protocols::wum::SilentStructStore totalpool;
	totalpool.clear();
	// Parent pool
	protocols::wum::SilentStructStore prvpool( structs );
	std::string const sim_replace_obj = option[ OptionKeys::lh::sim_replace_obj ]();
	std::string const similarity_method = option[ lh::similarity_method ]();

	// increment on obj_dominant cut
	if ( update_obj_cut ) {
		for ( core::Size i_obj = 1; i_obj <= nobjs(); ++i_obj ) {
			obj_dominant_cut_[i_obj] += obj_cut_increment_[i_obj];
		}
	}
	// report current obj_cut
	TR << "Current obj_dominate_cut: ";
	for ( core::Size i_obj = 1; i_obj <= nobjs(); ++i_obj ) {
		TR << "/ " << fobjnames_[i_obj] << " " << F(6,1,obj_dominant_cut_[i_obj]);
	}
	TR << std::endl;

	core::Size nbefore;
	nbefore = new_structs.size();
	filter_similar( new_structs, measure, simlimit, sim_replace_obj );

	//calculate_pool_diversity( new_structs );
	//filter_similar( new_structs, measure, simlimit, "similarity" );
	TR << "New structures, initial/filtered_by_similiarity:" << nbefore << " " << new_structs.size() << std::endl;
	totalpool.add( new_structs );

	// then filter out within the whole pool
	// but only for unexpired parents
	core::Size const expire_cut = option[ OptionKeys::lh::expire_after_rounds ]();
	for ( core::Size i_p = 0; i_p < structs.size(); ++i_p ) {
		core::io::silent::SilentStructOP ss = structs.get_struct( i_p );
		if ( ss->get_energy( "nuse" ) <= expire_cut ) {
			totalpool.add( ss );
		} else {
			TR << "Limit parent " << ss->decoy_tag() << " based on nuse cut." << std::endl;
		}
	}

	nbefore = totalpool.size();

	core::Size nmax_filt( nmax );
	//if( option[ lh::filter_up_to_maxlib ]() ) nmax_filt = nmax;
	filter_similar( totalpool, measure, simlimit, sim_replace_obj, nmax_filt );

	TR << "Total pool, initial/filtered_by_similarity: " << nbefore << " " << totalpool.size() << std::endl;

	if ( totalpool.size() <= nmax_filt ) {
		structs = totalpool;
		new_structs.clear();
		TR << "NSGAII call, but Direct addition due to empty spaces" << std::endl;
		return true;
	}

	// get list of unused pools

	// entering diversity

	// org:
	totalpool.all_add_energy( "similarity", 0.0 ); //just initialize
	if ( similarity_method != "frontiermin"  ) {
		calculate_pool_diversity( totalpool, totalpool );
		//calculate_pool_diversity( totalpool, prvpool );
	}

	TR << "update_library based on NSGAII, nprv/nnew/nfilt/nmax_filt = ";
	TR << structs.size() << "/" << new_structs.size() << "/";
	TR << totalpool.size() << "/" << nmax_filt << std::endl;

	core::Size time1 = time(nullptr);

	// Be careful - all the indices start from 0
	utility::vector0< protocols::wum::SilentStructStore > frontier( 1 );
	utility::vector0< utility::vector0< core::Size > > ifront( 1 );
	core::Size nstruct( totalpool.size() );
	utility::vector0< utility::vector0< core::Size > > Sp( nstruct );
	utility::vector0< core::Size > np( nstruct );
	utility::vector0< bool > is_left( nstruct, true );

	// Fast non-dominated-sort
	for ( core::Size i_p = 0; i_p < nstruct; ++i_p ) {
		SilentStructCOP p( totalpool.get_struct( i_p ) );
		Sp[i_p].resize( 0 );
		np[i_p] = 0;

		for ( core::Size i_q = 0; i_q < nstruct; ++i_q ) {
			if ( i_p == i_q ) continue;
			SilentStructCOP q( totalpool.get_struct( i_q ) );

			/*
			TR << "i_p/i_q/score/goap/sim: " << i_p << " " << i_q;
			for( core::Size iobj = 1; iobj <= fobjnames_.size(); ++iobj )
			TR << " " << p->get_energy( fobjnames_[iobj] ) - q->get_energy( fobjnames_[iobj] );
			TR << " " << is_dominant( p, q) << " " << is_dominant( q, p ) << std::endl;
			*/

			if ( is_dominant( p, q ) ) { // p dominates q
				Sp[i_p].push_back( i_q );
			} else if ( is_dominant( q, p ) ) { // q dominates p
				np[i_p]++; // N non-dominating
			}
		}

		// if not dominiated by any, set as first-frontier
		if ( np[i_p] == 0 ) {
			frontier[0].add( totalpool.get_struct( i_p ) );
			ifront[0].push_back( i_p );
			for ( core::Size j = 0; j < ifront[0].size(); ++j ) is_left[ ifront[0][j] ] = false;
		}
	}
	frontier[0].all_add_energy( "frontier", 0 );

	TR.Debug << "frontier1 (size " << ifront[0].size() << ": ";
	for ( core::Size j = 0; j < ifront[0].size(); ++j ) TR << " " << ifront[0][j];
	TR.Debug << ")" << std::endl;

	// Debug
	for ( core::Size j = 0; j < nstruct; ++j ) {
		TR.Debug << "i_p/np/Sp: " << j << " " << np[j] << " | ";
		for ( core::Size k = 0; k < Sp[j].size(); ++k ) TR.Debug << " " << Sp[j][k];
		TR.Debug << std::endl;
	}

	core::Size time2 = time(nullptr);

	// Assign frontiers
	core::Size i( 0 );
	core::Size n_in_frontier( frontier[0].size() );
	protocols::wum::SilentStructStore sorted_pool;

	while ( true ) {
		// Stop if assigned enough
		if ( n_in_frontier >= nmax ) break;

		sorted_pool.add( frontier[i] );

		// update diversity on currently added so far if sorted_pool is large enough
		// experimental
		if ( similarity_method == "frontiermin" ) {
			calculate_pool_diversity( totalpool, sorted_pool );
		}

		protocols::wum::SilentStructStore Q; // Next frontier
		utility::vector0< core::Size > iQ;
		protocols::wum::SilentStructStore const &P = frontier[i];
		utility::vector0< core::Size > const &iP = ifront[i];

		for ( core::Size istr = 0; istr < P.size(); ++istr ) {
			core::Size i_p = iP[istr];

			// When p pops out, p's next frontier members' dominant number reduces by 1
			for ( core::Size jstr = 0; jstr < Sp[i_p].size(); ++jstr ) {
				core::Size i_q = Sp[i_p][jstr];
				np[i_q]--; // mark that

				if ( np[i_q] == 0 ) { // If all the dominating members are already poped out
					Q.add( totalpool.get_struct( i_q ) );
					iQ.push_back( i_q );
				}
			}
		}

		i++;
		Q.all_add_energy( "frontier", i );
		frontier.push_back( Q );
		ifront.push_back( iQ );
		for ( core::Size j = 0; j < iQ.size(); ++j ) is_left[ iQ[j] ] = false;
		n_in_frontier += Q.size();

		TR.Debug << "frontier" << i+1 << "(size " << ifront[i].size() << ": ";
		for ( core::Size j = 0; j < ifront[i].size(); ++j ) TR.Debug << " " << ifront[i][j];
		TR.Debug << "), left: " << int(nmax) - int(n_in_frontier) << std::endl;
	}

	// Add into sorted_pool
	// up to n-1 frontier, store all
	//for( core::Size i_front = 0; i_front < frontier.size()-1; ++i_front ){
	//  sorted_pool.add( frontier[i_front] );
	//}

	// for the last frontier, sort by certain objective (e.g. the first objective)
	// Make sure your objective is reasonable for doing this;
	protocols::wum::SilentStructStore &last_frontier = frontier[ frontier.size()-1 ];
	//utility::vector0< core::Size > &ilast = ifront[ frontier.size()-1 ];

	core::Size time3 = time(nullptr);

	//last_frontier.sort_by( fobjnames(1) );
	//last_frontier.sort_by( "goap" );
	// try using "similarity to the already selected pool" for picking among the last frontier
	for ( core::Size j = 0; j < last_frontier.size(); ++j ) {
		if ( sorted_pool.size() > 0 ) {
			calculate_structure_diversity( last_frontier.get_struct(j), sorted_pool );
		} else {
			calculate_structure_diversity( last_frontier.get_struct(j), last_frontier );
		}
	}
	last_frontier.sort_by( "similarity" );
	core::Size time4 = time(nullptr);

	core::Size const nleft( nmax - sorted_pool.size() );

	if ( frontier.size() == 1 ) {
		TR << "NSGAII: Adding on pool, " << sorted_pool.size() << " from first frontier" << std::endl;
	} else {
		TR << "NSGAII: Adding on pool, " << sorted_pool.size() << " from first " << frontier.size() - 2 << " frontiers";
		TR << " and " << nleft << " from the last frontier" << std::endl;
	}

	TR << "Time required for NSGA (front 1/n-1/n): " << time2-time1 << " " << time3-time2 << " " << time4-time2 << std::endl;

	// new_structs will return the removed ones
	new_structs.clear();
	for ( core::Size i_left = 0; i_left < last_frontier.size(); ++i_left ) {
		if ( i_left < nleft ) {
			sorted_pool.add( last_frontier.get_struct( i_left ) );
		} else {
			new_structs.add( last_frontier.get_struct( i_left ) );
		}
	}

	// remaining frontiers
	core::Size const final_frontier( frontier.size() );
	for ( core::Size j = 0; j < nstruct; ++j ) {
		if ( is_left[j] ) {
			totalpool.get_struct( j )->add_energy( "frontier", final_frontier );
			new_structs.add( totalpool.get_struct( j ) );
		}
	}

	// Update library
	structs.clear();
	structs.add( sorted_pool );

	// finally add diversity info : different from number used for frontier decision,
	// but will give estimation to the next iteration
	//calculate_pool_diversity( structs );

	return true;
}

// only for CSA
void
MultiObjective::succeed_substitute_info( core::io::silent::SilentStructOP ss_sub,
	core::io::silent::SilentStructCOP ss_ref,
	bool const reset
) const
{
	core::Size const nuse = reset? 0:ss_ref->get_energy( "nuse" );
	core::Size const poolid( ss_ref->get_energy( "poolid" ) );

	ss_sub->add_energy( "nuse", nuse );
	ss_sub->add_energy( "poolid", poolid );
}

void
MultiObjective::filter_similar( protocols::wum::SilentStructStore &structs,
	std::string const measure,
	core::Real const criteria,
	std::string const score_for_priority,
	core::Size const nmax
)
{

	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size const nstruct( structs.size() );
	utility::vector0< bool > filtered( nstruct, false );
	utility::vector0< core::Size > nuse( nstruct, 0 );

	core::Size nfilt( 0 );
	core::Size nfilt_max = (nstruct >= nmax) ? nstruct - nmax : 0;

	TR << "start filtering from nstruct=" << structs.size() << ", " << "DistMeasure/Dcut: "
		<< measure << "/" << F(8,3,criteria) << std::endl;

	for ( core::Size i = 1; i < nstruct; ++i ) {
		SilentStructOP ss1 = structs.get_struct( i );
		core::Real score1 = ss1->get_energy( score_for_priority );
		std::string const iname = ss1->decoy_tag();

		for ( core::Size j = 0; j < i; ++j ) {
			SilentStructOP ss2 = structs.get_struct( j );
			core::Real score2 = ss2->get_energy( score_for_priority );
			std::string const jname = ss2->decoy_tag();

			if ( nfilt >= nfilt_max ) break;

			if ( filtered[j] ) continue;


			// only Sscore for now
			bool is_similar( false );
			core::Real dist = distance( ss1, ss2, measure, false );
			if ( measure == "Sscore" ) dist = 1.0 - dist;

			if ( dist < criteria ) is_similar = true;
			if ( !is_similar ) continue;

			// get larger nuse if similar
			if ( score1 < score2 ) {
				filtered[j] = true;
				nfilt++;
				nuse[i] = nuse[j] > nuse[i] ? nuse[j] : nuse[i];
				TR << "remove " << jname << " on " << iname << ", dist/dE " << F(8,3,dist) << "/" << F(8,3,score2 - score1) << std::endl;
				break;
			} else {
				filtered[i] = true;
				nfilt++;
				nuse[j] = nuse[j] > nuse[i] ? nuse[j] : nuse[i];
				TR << "remove " << iname << " on " << jname << ", dist/dE " << F(8,3,dist) << "/" << F(8,3,score1 - score2) << std::endl;
			}
		} // j

		if ( nfilt >= nfilt_max ) break;
	} // i

	core::Size i( 0 );
	protocols::wum::SilentStructStore outstructs;
	for ( protocols::wum::SilentStructStore::const_iterator it = structs.begin();
			it != structs.end(); ++it, ++i ) {
		if ( !filtered[i] ) {
			core::io::silent::SilentStructOP ss( *it );
			ss->add_energy( "nuse", nuse[i] ); //updated based on max within similar
			outstructs.add( ss );
		}
	}

	structs = outstructs;
}

void
MultiObjective::add_objective_function_info(
	core::io::silent::SilentStructOP ss,
	protocols::wum::SilentStructStore & ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Get a copy
	core::io::silent::SilentStructCOP ss_copy = ss->clone();

	// First cleanup useless information so as to make it "readable"
	ss->clear_energies();

	// Equally add information from ss no matter how they are generated
	ss->add_energy( "state", ss_copy->get_energy("state") );
	ss->add_energy( "samplemethod", ss_copy->get_energy("samplemethod") );
	ss->add_energy( "round", ss_copy->get_energy("round") );
	ss->add_energy( "nuse", ss_copy->get_energy("nuse") );
	ss->add_energy( "iter", ss_copy->get_energy("iter") );
	ss->add_energy( "NMmode", ss_copy->get_energy("NMmode") );
	ss->add_energy( "NMscale", ss_copy->get_energy("NMscale") );

	// copy loopinfo
	core::Size nloop = (core::Size)( ss->get_energy( "nloop" ) );
	if ( nloop > 0 ) {
		//ss->add_string_value( "loops_nsampled", ss_copy->get_string_value( "loops_nsampled" ) );
		for ( core::Size iloop = 1; iloop <= nloop; ++iloop ) {
			std::stringstream sstream( "" );
			sstream << "nsampled_loop" << iloop;
			ss->add_energy( sstream.str(), ss_copy->get_energy( sstream.str() ) );
		}
	}

	// Get native info
	core::chemical::ResidueTypeSetCOP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueTypeSetCOP rsd_set_cen
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	core::pose::Pose pose0;
	core::import_pose::pose_from_file( pose0, *rsd_set, option[ in::file::s ](1) , core::import_pose::PDB_file);
	add_poseinfo_to_ss( *ss, pose0, "_i" );

	if ( option[ in::file::native ].user() ) {
		core::pose::Pose native_pose;
		core::import_pose::pose_from_file( native_pose, *rsd_set, option[ in::file::native ]() , core::import_pose::PDB_file);
		add_poseinfo_to_ss( *ss, native_pose, "" );
	}

	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	// Iter over multi objective functions
	for ( core::Size i_obj = 1; i_obj <= nobjs(); ++ i_obj ) {
		std::string score_name( fobjnames( i_obj ) );
		core::Real score_i = ss->get_energy( score_name );
		//TR << "Fobj: " << score_name << std::endl;
		// Calculate energy if empty yet
		if ( score_name == "score" ) {
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set );
			objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->energies_from_pose( pose_tmp );

		} else if ( score_name == "similarity" ) {
			continue;

		} else if ( score_name == "cst_cen" ) {
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set_cen );
			core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose_tmp );
			score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->add_energy( "cst_cen", score_i );

		} else if ( score_name == "cst_fa" ) {
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp, *rsd_set );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose_tmp );
			score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->add_energy( "cst_fa", score_i );

		} else if ( score_i == 0.0 && score_name != "score" ) {
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp );
			score_i = objsfxnOPs_[i_obj]->score( pose_tmp );
			ss->add_energy( score_name, score_i );
		}
	}

	// add ediff if both score and goap exists
	if ( fobjnames_.contains("score") && fobjnames_.contains("goap") ) {
		ss->add_energy( "ediff", ss->get_energy("goap") - ss->get_energy("score") );
		ss->add_energy( "esum", ss->get_energy("goap") + ss->get_energy("score") );
	}

	// constraints

	// finally add diversity info
	// skip?
	//calculate_structure_diversity( ss, sstore );
}

void
MultiObjective::add_objective_function_info( protocols::wum::SilentStructStore & sstore ) const
{
	core::Size const nstruct( sstore.size() );
	for ( core::Size i = 0; i < nstruct; ++i ) {
		add_objective_function_info( sstore.get_struct( i ), sstore );
	}
}


std::string
MultiObjective::formatted_objs_values( core::io::silent::SilentStruct const &ss ) const
{
	std::stringstream sstream;

	for ( core::Size iobj = 1; iobj <= nobjs(); ++iobj ) {
		std::string const score_name( fobjnames(iobj) );
		sstream << " | " << F(6,1, ss.get_energy( score_name ) );
	}
	return sstream.str();
}

std::string
MultiObjective::formatted_objs_names() const
{
	std::stringstream sstream;

	for ( core::Size iobj = 1; iobj <= nobjs(); ++iobj ) {
		std::string const scorename( fobjnames(iobj) );
		sstream << " " << std::setw(7) << scorename;
	}
	return sstream.str();
}

}
}
