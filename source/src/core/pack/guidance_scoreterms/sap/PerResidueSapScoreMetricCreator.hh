// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/PerResidueSapScoreMetricCreator.hh
/// @brief Creators for PerResidueSapScoreMetric
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetricCreator_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_PerResidueSapScoreMetricCreator_hh

#include <core/simple_metrics/SimpleMetricCreator.hh>

// Protocol headers
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>



namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


class PerResidueSapScoreMetricCreator : public core::simple_metrics::SimpleMetricCreator {
public:


	/// @brief Instantiate a particular SimpleMetric
	core::simple_metrics::SimpleMetricOP
	create_simple_metric() const override;

	/// @brief Return a string that will be used to instantiate the particular SimpleMetric
	std::string
	keyname() const override;

	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition &) const override;
};

} //sap
} //guidance_scoreterms
} //pack
} //core

#endif
