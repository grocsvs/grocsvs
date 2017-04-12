from __future__ import print_function

import argparse
import collections
import json
import logging
import sys

from grocsvs import options as svoptions
from grocsvs import log
from grocsvs import pipeline
from grocsvs import utilities
from grocsvs import stages as svstages


logging.basicConfig(format='%(message)s', level=logging.DEBUG)



def ready_output_dir(options):
    utilities.ensure_dir(options.working_dir)
    utilities.ensure_dir(options.results_dir)
    utilities.ensure_dir(options.log_dir)


def run(options):
    """
    1. create output directories
    2. collect args for each stage
    3. check which stages need to run
    4. iterate through stages and submit jobs
    5. validate that we're done running
    """

    svoptions.validate_options(options)
    ready_output_dir(options)

    stages = get_stages()
    runner = pipeline.Runner(options)

    for stage_name, stage in stages.items():
        runner.run_stage(stage, stage_name)


def clean(options, clean_stage_name):
    stages = get_stages()

    if not clean_stage_name in stages:
        print('*'*20, "ERROR", '*'*20)
        print('Error: unknown stage "{}". Stage must be one of the '.format(clean_stage_name))
        print('following (remember to include surrounding quotation marks):')
        for i, stage_name in enumerate(stages):
            print('{:>3} - "{}"'.format(i+1, stage_name))
        sys.exit(1)

        
    doclean = False
    for stage_name, stage in stages.items():
        if doclean or stage_name == clean_stage_name:
            doclean = True

            stage.clean_all_steps(options)


def get_stages():
    stages = collections.OrderedDict()


    # Pre-processing

    stages["Preflight"] = svstages.preflight.PreflightStep

    stages["Constants"] = svstages.constants.ConstantsStep
    stages["Estimate Read Cloud Parameters"] = svstages.call_readclouds.EstimateReadCloudParamsStep
    stages["Call Read Clouds"] = svstages.call_readclouds.CallReadcloudsStep
    stages["Combine Read Clouds"] = svstages.call_readclouds.CombineReadcloudsStep

    # stages["Filter Fragments"] = svstages.filter_fragments.FilterFragmentsStep
    stages["Sample Info"] = svstages.sample_info.SampleInfoStep
    stages["QC"] = svstages.qc.QCStep


    # Find SV candidates

    stages["Window Barcodes"] = svstages.window_barcodes.WindowBarcodesStep
    stages["Barcode Overlaps"] = svstages.barcode_overlaps.BarcodeOverlapsStep

    stages["SV Candidate Regions"] = \
        svstages.sv_candidate_regions.SVCandidateRegionsStep
    stages["SV Candidates From Regions"] = \
        svstages.sv_candidates.SVCandidatesStep


    ### Initial clustering ###
    
    stages["Refine Breakpoints"] = \
        svstages.refine_grid_search_breakpoints.RefineGridSearchBreakpointsStep
    stages["Combine Refined Breakpoints"] = \
        svstages.refine_grid_search_breakpoints.CombineRefinedBreakpointsStep
    stages["Cluster SVs"] = svstages.cluster_svs.ClusterSVsStep


    ### Assembly ###

    stages["Barcodes From Graphs"] = \
        svstages.barcodes_from_graphs.BarcodesFromGraphsStep

    stages["Collect Reads for Barcodes"] = \
        svstages.collect_reads_for_barcodes.CollectReadsForBarcodesStep

    
    stages["Perform Assembly"] = svstages.assembly.AssemblyStep
    stages["Walk Assemblies"] = svstages.walk_assemblies.WalkAssembliesStep
    stages["Postassembly Merge"] = \
        svstages.postassembly_merge.PostAssemblyMergeStep


    ### Final clustering ###

    stages["Supporting Barcodes"] = \
        svstages.supporting_barcodes.SupportingBarcodesStep
    stages["Pair Evidence"] = svstages.pair_evidence.PairEvidenceStep

    stages["Final Refine"] = \
        svstages.refine_breakpoints.RefineBreakpointsWithAssembliesStep
    stages["Final Cluster SVs"] = svstages.final_clustering.FinalClusterSVsStep


    ### Multi-sample genotyping and postprocessing ###

    stages["Genotyping"] = svstages.genotyping.GenotypingStep
    stages["Merge Genotypes"] = svstages.genotyping.MergeGenotypesStep
    
    #stages["Visualize"] = svstages.visualize.VisualizeStep
    stages["Postprocessing"] = svstages.postprocessing.PostprocessingStep

    return stages


def load_config(config_path):
    try:
        config = json.load(open(config_path))
    except ValueError as err:
        print("Error parsing configuration file '{}': '{}'\n  Check that this is a properly formatted JSON file!".format(config_path, err))
        sys.exit(1)
    options = svoptions.Options.deserialize(config, config_path)
    return options

    
def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Genome-wide Reconstruction of Complex Structural Variants")
    parser.add_argument("config", help="Path to configuration.json file")
    parser.add_argument("--restart", metavar="FROM-STAGE", help="restart from this stage")
    parser.add_argument("--local", action="store_true", help="run locally in single processor mode")
    parser.add_argument("--multiprocessing", action="store_true", help="run locally using multiprocessing")
    parser.add_argument("--debug", action="store_true", help="run in debug mode")

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    
    options = load_config(args.config)
    options.debug = args.debug

    print(options)

    if args.local:
        options.cluster_settings = svoptions.ClusterSettings()
    if args.multiprocessing:
        options.cluster_settings = svoptions.ClusterSettings()
        options.cluster_settings.cluster_type = "multiprocessing"

    if args.restart is not None:
        clean(options, args.restart)

    log.log_command(options, sys.argv)

    return options


def main():
    options = parse_arguments(sys.argv[1:])
    run(options)



if __name__ == '__main__':
    main()
