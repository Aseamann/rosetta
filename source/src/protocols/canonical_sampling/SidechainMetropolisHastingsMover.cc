// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/SidechainMetropolisHastingsMover.cc
/// @brief SidechainMetropolisHastingsMover methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/SidechainMetropolisHastingsMoverCreator.hh>

#include <protocols/simple_moves/sidechain_moves/SidechainMoverBase.hh>

// protocols headers
#include <basic/datacache/DataMap.fwd.hh>

#include <protocols/jd2/util.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>

// core headers

#include <core/pose/Pose.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>


// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

#include <protocols/canonical_sampling/TemperatureController.hh> // AUTO IWYU For TemperatureController

// C++ Headers

static basic::Tracer tr( "protocols.canonical_sampling.SidechainSidechainMetropolisHastingsMover" );

namespace protocols {
namespace canonical_sampling {

using namespace core;
using namespace scoring;

std::string
SidechainMetropolisHastingsMoverCreator::keyname() const {
	return SidechainMetropolisHastingsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SidechainMetropolisHastingsMover >();
}

std::string
SidechainMetropolisHastingsMoverCreator::mover_name() {
	return "SidechainMetropolisHastings";
}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover() :
	stride_( 1000 )
{}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover( core::Size stride ) :
	stride_( stride )
{}

SidechainMetropolisHastingsMover::SidechainMetropolisHastingsMover( SidechainMetropolisHastingsMover const & ) = default;

SidechainMetropolisHastingsMover::~SidechainMetropolisHastingsMover() = default;

bool
SidechainMetropolisHastingsMover::pass_metropolis(core::Real delta_energy, core::Real last_proposal_density_ratio ) const {
	core::Real boltz_factor = delta_energy / monte_carlo()->temperature();
	if ( tr.Trace.visible() ) {
		tr.Trace << " temperature: " << monte_carlo()->temperature()
			<< " deltaE= " << delta_energy
			<< " boltzman=" << boltz_factor
			<< " lpd= " << last_proposal_density_ratio << std::endl;
	}
	core::Real probability = std::exp( std::min( 40.0, std::max( -40.0, boltz_factor ))) *  last_proposal_density_ratio ;
	if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {
		return false;
	} else {
		return true;
	}
}

void
SidechainMetropolisHastingsMover::apply( core::pose::Pose & pose )
{
	prepare_simulation( pose );

	scoring::ScoreFunction const& sfxn = monte_carlo()->score_function();
	pack::interaction_graph::SimpleInteractionGraphOP ig;
	ig = utility::pointer::make_shared< pack::interaction_graph::SimpleInteractionGraph >(); //commented out debug
	ig->set_scorefunction( sfxn );

	utility::vector1< Real > new_chi;
	Real current_energy = sfxn(pose);

	// std::string const traj_file_tag( jd2::current_output_name() );
	//  counters_.reset();

	//ek for fast sidechain sampling and internal mc trials
	utility::vector1< conformation::ResidueOP > current;
	// utility::vector1< conformation::ResidueOP > previous;
	// utility::vector1< pack::dunbrack::ChiVector > chi_vectors;
	// utility::vector1< pack::dunbrack::RotVector > rot_vectors;


	current.resize(pose.size());
	// previous.resize(pose.size());

	// rot_vectors.resize( pose.size() );
	// chi_vectors.resize( pose.size() );

	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		current[ i ] = utility::pointer::make_shared< core::conformation::Residue >( pose.residue( i ) );
	}


	runtime_assert( ig != nullptr );

	ig->initialize( pose );
	Real last_accepted_prop_density( 1.0 );
	Real last_accepted_dE( 0.0 );
	for ( core::Size ct = 1; ct <= ntrials(); ct++ ) {
		protocols::simple_moves::sidechain_moves::SidechainMoverBaseOP move = utility::pointer::dynamic_pointer_cast< protocols::simple_moves::sidechain_moves::SidechainMoverBase > ( random_mover() );
		runtime_assert( move != nullptr ); //fow now only Sidechain Movers...

		core::Size resid = move->suggest_residue_number( pose );
		conformation::ResidueOP new_state( new conformation::Residue( pose.residue( resid ) ) );
		new_state = move->make_move( new_state );
		set_last_move( move );

		Real delta_energy = ig->consider_substitution( resid, new_state, *new_state->nonconst_data_ptr() );
		if ( pass_metropolis( delta_energy, move->last_proposal_density_ratio() ) ) { //ek
			ig->commit_change( resid );
			current_energy -= delta_energy;
			current[ resid ] = new_state;
			set_last_accepted( true );
			last_accepted_prop_density = move->last_proposal_density_ratio();
			last_accepted_dE = delta_energy;
		} else { //rejected metropolis criterion
			ig->reject_change( resid, new_state, *new_state->nonconst_data_ptr() );
			set_last_accepted( false );
		}

		tempering()->temperature_move( current_energy );
		move->observe_after_metropolis( *this );

		core::Size model_count( output_count( ct ) );
		if ( model_count ) {
			for ( core::Size res_i = 1; res_i <= current.size(); res_i++ ) {
				pose.replace_residue( res_i, (*current[ res_i ]), true );
			}
			core::Real const score( sfxn( pose ) );
			if ( std::abs( score-current_energy ) > 1 ) { //threshold 0.1 gives a couple warnings -- but it never drifts apart
				tr.Warning << "Energy mismatch!!! score=" << score << " ig->energy " << current_energy << std::endl;
			}
			nonconst_monte_carlo().set_last_accepted_pose( pose );
			if ( score < nonconst_monte_carlo().lowest_score() ) {
				nonconst_monte_carlo().set_lowest_score_pose( pose );
			}
			protocols::jd2::add_string_real_pair_to_current_job( "prop_density", move->last_proposal_density_ratio() );
			protocols::jd2::add_string_real_pair_to_current_job( "prop_density_accept", last_accepted_prop_density );
			protocols::jd2::add_string_real_pair_to_current_job( "move_dE", delta_energy );
			protocols::jd2::add_string_real_pair_to_current_job( "move_dE_accept", last_accepted_dE );
			protocols::jd2::add_string_string_pair_to_current_job( "move_type", move->type() );
		}

		for ( core::Size i = 1; i <= observers().size(); ++i ) {
			if ( observers()[ i ]->requires_pose() && !model_count ) continue;
			observers()[i]->observe_after_metropolis(*this);
		}

	} //for ntrials

	wind_down_simulation( pose );
}

core::Size
SidechainMetropolisHastingsMover::output_count( core::Size ct ) const {
	if ( ct % stride_ == 0 ) {
		return ct / stride_;
	} else return 0;
}

std::string
SidechainMetropolisHastingsMover::get_name() const
{
	return "SidechainMetropolisHastingsMover";
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMover::clone() const
{
	return utility::pointer::make_shared< protocols::canonical_sampling::SidechainMetropolisHastingsMover >(*this);
}

protocols::moves::MoverOP
SidechainMetropolisHastingsMover::fresh_instance() const
{
	return utility::pointer::make_shared< SidechainMetropolisHastingsMover >();
}

void
SidechainMetropolisHastingsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	stride_ = tag->getOption< core::Size >( "stride", stride_ );
	Parent::parse_my_tag( tag, data );
}


} //moves
} //protocols

