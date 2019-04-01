from collections import defaultdict
import json
import math
import os

import idseq_dag.util.m8 as m8
import idseq_dag.util.log as log
import idseq_dag.util.command as command
from idseq_dag.engine.pipeline_step import PipelineStep

MAX_NUM_BINS_COVERAGE = 2000
MIN_CONTIG_SIZE = 4

# Transform a range from one scale to another.
# The 'first_range' should be within the interval [first_scale_start, first_scale_end]
# For example, transform_range((3, 5) 0, 10, 100, 200) would return (130, 150)
def transform_range(
    first_range, first_scale_start, first_scale_end, second_scale_start, second_scale_end
):
    def transform_value(value):
        position = (value - first_scale_start) / (first_scale_end - first_scale_start)
        return position * (second_scale_end - second_scale_start) + second_scale_start

    return [transform_value(value) for value in first_range]

def align_range(bound_one, bound_two):
    return (min(bound_one, bound_two), max(bound_one, bound_two))

# Turn (1, 5) into (0, 5) and (5, 1) into (5, 0).
def decrement_lower_bound(bound_one, bound_two):
    if bound_one < bound_two:
        return (bound_one - 1, bound_two)
    else:
        return (bound_one, bound_two - 1)

class PipelineStepGenerateCoverageViz(PipelineStep):
    """Pipeline step to generate JSON file for read alignment visualizations to
    be consumed by the web app.
    """

    def run(self):
        (_blast_m8, reassigned_m8, hit_summary, _refined_counts, _contig_summary, blast_top_m8) = self.input_files_local[0]
        (accession_summary_tab, ) = self.input_files_local[1]
        (contig_coverage_json, ) = self.input_files_local[2]
        (contig_stats_json, ) = self.input_files_local[3]
        (gsnap_deduped_m8, ) = self.input_files_local[4]

        coverage_viz_summary = self.output_files_local()[0]

        valid_contigs = self.process_contig_stats_json(contig_stats_json)

        (accessions_map, taxons_to_accessions, unassigned_reads) = self.process_hitsummary(hit_summary, valid_contigs)

        self.process_accession_summary_tab(accession_summary_tab, accessions_map, taxons_to_accessions)

        contigs_map = self.process_contig_blast_m8(blast_top_m8, valid_contigs)
        reads_map = self.process_deduped_m8(gsnap_deduped_m8, unassigned_reads)

        self.process_contig_coverage_json(contig_coverage_json, contigs_map)

        coverage_viz_obj = self.construct_coverage_viz_obj(accessions_map, contigs_map, reads_map)


        with open(coverage_viz_summary, 'w') as cvs:
            json.dump(taxons_to_accessions, cvs)

        # Create a JSON file for each accession
        coverage_viz_dir = os.path.join(self.output_dir_local, "coverage_viz")
        command.execute(f"mkdir -p {coverage_viz_dir}")
        for accession_id in coverage_viz_obj:
            upload_file = os.path.join(coverage_viz_dir, f"{accession_id}_coverage_viz.json")

            with open(upload_file, 'w') as uf:
                json.dump(coverage_viz_obj[accession_id], uf)

            self.additional_files_to_upload.append(upload_file)

    def process_contig_stats_json(self, contig_stats_json):
        valid_contigs = {}
        with open(contig_stats_json, 'r') as csj:
            contig_read_counts = json.load(csj)

            for contig in contig_read_counts:
                if contig != "*" and contig_read_counts[contig] >= MIN_CONTIG_SIZE:
                    valid_contigs[contig] = True

        return valid_contigs

    def process_hitsummary(self, hit_summary, valid_contigs):
        accessions_map = defaultdict(lambda: {'reads': [], 'contigs': set() })
        # We only allow species taxons for now.
        taxons_to_accessions = defaultdict(set)

        line_count = 0
        with open(hit_summary, 'r') as hs:
            for line in hs:
                line_count += 1
                if line_count % 100000 == 0:
                    log.write(f"{line_count} lines in the hitsummary file processed.")

                values = line.rstrip().split("\t")

                # Only add as contig if the contig is "valid", i.e. it has 4 or more reads.
                if len(values) == 12 and values[7] in valid_contigs:
                    taxons_to_accessions[values[9]].add(values[8])
                    accessions_map[values[8]]["contigs"].add(values[7])
                else:
                    taxons_to_accessions[values[4]].add(values[3])
                    accessions_map[values[3]]["reads"].append(values[0])

        taxons_to_remove = []

        for taxon, accessions in taxons_to_accessions.items():
            total_accessions = len(accessions)
            total_contigs = sum(list(map(lambda accession_id: len(accessions_map[accession_id]["contigs"]), accessions)))
            contigs = list(map(lambda accession_id: accessions_map[accession_id]["contigs"], accessions))
            total_reads = sum(list(map(lambda accession_id: len(accessions_map[accession_id]["reads"]), accessions)))

            # print(f"{taxon} accessions: {total_accessions}, contigs: {total_contigs}, reads: {total_reads}, {contigs}")

            if total_contigs == 0:
                taxons_to_remove.append(taxon)

            taxons_to_accessions[taxon] = list(map(lambda accession_id: {
                "id": accession_id,
                "n_contig": len(accessions_map[accession_id]["contigs"]),
                "n_read": len(accessions_map[accession_id]["reads"])
            }, accessions))

        # Remove all taxons with no mapped contigs.
        for taxon in taxons_to_remove:
            accessions = list(taxons_to_accessions[taxon])
            del taxons_to_accessions[taxon]
            for accession_arr in accessions:
                del accessions_map[accession_arr["id"]]

        unassigned_reads = {}
        for _, accession_obj in accessions_map.items():
            for read in accession_obj["reads"]:
                unassigned_reads[read] = True

        return (accessions_map, taxons_to_accessions, unassigned_reads)

    def process_contig_blast_m8(self, blast_top_m8, valid_contigs):
        contigs = {}
        for _contig_id, _accession_id, _percent_id, _alignment_length, _e_value, _bitscore, line in m8.iterate_m8(blast_top_m8):
            parts = line.split("\t")
            name_parts = parts[0].split("_")

            if parts[0] in valid_contigs:
                contigs[parts[0]] = {
                    "total_length": int(name_parts[3]),
                    "accession": parts[1],
                    "query_start": int(parts[6]),
                    "query_end": int(parts[7]),
                    "subject_start": int(parts[8]),
                    "subject_end": int(parts[9]),
                }

        return contigs

    # We process gsnap.deduped.m8 instead of gsnap.reassigned.m8 because we ignore contigs with read_count < 4.
    # However, these contigs still get reassigned in gsnap.reassigned.m8,
    # and overwrite the original read alignment to the accession, which we want.
    def process_deduped_m8(self, deduped_m8, unassigned_reads_map):
        reads = {}
        for _read_id, _accession_id, _percent_id, _alignment_length, _e_value, _bitscore, line in m8.iterate_m8(deduped_m8):
            parts = line.split("\t")

            if parts[0] in unassigned_reads_map:

                reads[parts[0]] = {
                    # Extract the length from the contig name.
                    "accession": parts[1],
                    "query_start": int(parts[6]),
                    "query_end": int(parts[7]),
                    "subject_start": int(parts[8]),
                    "subject_end": int(parts[9]),
                }

        return reads

    def process_contig_coverage_json(self, contig_coverage_json, contigs_map):
        with open(contig_coverage_json, 'r') as ccj:
            contig_coverage = json.load(ccj)

            for contig_name in contig_coverage:
                if contig_name in contigs_map:
                    contigs_map[contig_name]["coverage"] = contig_coverage[contig_name]["coverage"]

    def process_accession_summary_tab(self, accession_summary_tab, accession_map, taxons_to_accessions):
        with open(accession_summary_tab, 'r') as ast:
            for line in ast:
                values = line.rstrip().split("\t")

                if values[0] in accession_map:
                    accession_map[values[0]]["total_length"] = int(values[2])
                    accession_map[values[0]]["name"] = values[1]

        for taxon, accessions in taxons_to_accessions.items():
            for accession in accessions:
                if accession["id"] in accession_map:
                    accession["name"] = accession_map[accession["id"]]["name"]

    def construct_coverage_viz_obj(self, accessions_map, contigs_map, reads_map):
        coverage_viz_obj = {}

        def get_contig_obj(contig_name, accession_id):
            if contig_name in contigs_map:
                contig_obj = contigs_map[contig_name]
                if contig_obj["accession"] != accession_id:
                    log.write(f"Mismatched contig accessions for {contig_name}: {contig_obj['accession']} (blast) versus {accession_id} (hitsummary)")
                    return None

                # For now, we just return the aligned portion on the accession,
                # even if it's only part of the total contig length.
                return [
                    contig_name,
                    contig_obj["subject_start"],
                    contig_obj["subject_end"]
                ]
            else:
                log.write(f"Could not find contig: {contig_name}")
                # We keep the None so it's easier to see that something went wrong in the output file.
                return None

        def get_read_obj(read_name, accession_id):
            if read_name in reads_map:
                read_obj = reads_map[read_name]
                # hitsummary is more strict than reassigned, since we judge whether a "hit" occurred.
                # Sometimes reassigned will have a value, but hitsummary won't.
                if read_obj["accession"] != accession_id:
                    log.write(f"Mismatched read accessions for {read_name}: {read_obj['accession']} (reassigned) versus {accession_id} (hitsummary)")
                    return None

                # For now, we just return the aligned portion on the accession,
                # even if it's only part of the total read length.
                return [
                    read_name,
                    read_obj["subject_start"],
                    read_obj["subject_end"]
                ]
            else:
                log.write(f"Could not find read: {read_name}")
                # We keep the None so it's easier to see that something went wrong in the output file.
                return None


        for accession_id, accession_data in accessions_map.items():
            total_length = accession_data["total_length"]
            contigs = list(map(lambda contig_name: get_contig_obj(contig_name, accession_id), accession_data["contigs"]))
            reads = list(map(lambda read_name: get_read_obj(read_name, accession_id), accession_data["reads"]))

            (coverage, coverage_bin_size) = PipelineStepGenerateCoverageViz.calculate_accession_coverage(
                accession_id, accession_data, contigs_map, reads_map
            )

            coverage_viz_obj[accession_id] = {
                "total_length": total_length,
                "name": accession_data["name"],
                "contigs": contigs,
                "reads": reads,
                "coverage": coverage,
                "coverage_bin_size": coverage_bin_size
            }

        return coverage_viz_obj

    @staticmethod
    def calculate_accession_coverage(accession_id, accession_data, contigs_map, reads_map, max_num_bins=MAX_NUM_BINS_COVERAGE):
        # Number of bins to calculate coverage for.
        num_bins = min(max_num_bins, accession_data["total_length"])
        bin_size = accession_data["total_length"] / num_bins
        coverage = [0] * num_bins

        for contig_name in accession_data["contigs"]:
            contig_obj = contigs_map[contig_name]
            # Ignore contigs with accession mismatch
            if contig_obj["accession"] != accession_id:
                continue

            # The bins and coverage array are 0-indexed, but subject start/end and coverage start/end are 1-indexed.
            # We convert everything to 0-index here and stay in 0-index for the rest of the code.
            # NOTE: We decrement only the lower bound here basically so that we can treat the discrete integer indices as a continuous interval.
            # We convert back to integer indices when we calculate coverage_arr_start/_end.
            (subject_start, subject_end) = decrement_lower_bound(contig_obj["subject_start"], contig_obj["subject_end"])
            (query_start, query_end) = decrement_lower_bound(contig_obj["query_start"], contig_obj["query_end"])

            # Find all bins that this contig overlaps, and calculate average coverage for each bin separately.
            (bin_start, bin_end) = align_range(subject_start / bin_size, subject_end / bin_size)

            for i in range(math.floor(bin_start), math.ceil(bin_end)):
                # Get the section of the accession that corresponds to the current bin and overlaps with the contig.
                accession_range = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

                # Convert the accession range to a contig range by using the alignment values.
                contig_range = transform_range(accession_range, subject_start, subject_end, query_start, query_end)

                # The contig coverage array should be the same length as the contig length.
                # If not, convert to the appropriate range in the coverage array.
                if contig_obj["total_length"] == len(contig_obj["coverage"]):
                    coverage_range = align_range(contig_range[0], contig_range[1])
                else:
                    coverage_range = transform_range(contig_range, 0, contig_obj["total_length"], 0, len(contig_obj["coverage"]))
                    coverage_range = align_range(coverage_range[0], coverage_range[1])

                (coverage_arr_start, coverage_arr_end) = (math.floor(coverage_range[0]), math.ceil(coverage_range[1]))

                # Get the average coverage for the section of the contig that overlaps with this bin.
                avg_coverage_for_coverage_range = sum(contig_obj["coverage"][coverage_arr_start: coverage_arr_end]) / (coverage_arr_end - coverage_arr_start)

                # Multiply by the proportion of the bin that the contig covers.
                avg_coverage_for_bin = avg_coverage_for_coverage_range * (abs(accession_range[1] - accession_range[0]) / bin_size)

                coverage[i] += avg_coverage_for_bin

        # The processing for reads is similar to contigs above, but the avg coverage on the read is simply 1.
        for read_name in accession_data["reads"]:
            read_obj = reads_map[read_name]
            # Ignore contigs with accession mismatch
            if read_obj["accession"] != accession_id:
                continue

            # The bins and coverage array are 0-indexed, but subject start/end and coverage start/end are 1-indexed.
            # We convert everything to 0-index here and stay in 0-index for the rest of the code.
            # NOTE: We decrement only the lower bound here basically so that we can treat the discrete integer indices as a continuous interval.
            # We convert back to integer indices when we calculate coverage_arr_start/_end.
            (subject_start, subject_end) = decrement_lower_bound(read_obj["subject_start"], read_obj["subject_end"])

            # Find all bins that this contig overlaps, and calculate average coverage for each bin separately.
            (bin_start, bin_end) = align_range(subject_start / bin_size, subject_end / bin_size)

            for i in range(math.floor(bin_start), math.ceil(bin_end)):
                # Get the section of the accession that corresponds to the current bin and overlaps with the read.
                accession_range = [bin_size * max(bin_start, i), bin_size * min(bin_end, i + 1)]

                # The read coverage is 1. Multiply by the proportion of the bin that the read covers.
                avg_coverage_for_bin = (abs(accession_range[1] - accession_range[0]) / bin_size)

                coverage[i] += avg_coverage_for_bin

        coverage = [round(val, 1) for val in coverage]

        return (coverage, bin_size)
