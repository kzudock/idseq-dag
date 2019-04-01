import json
import os
import math
import threading
import time
import traceback
import subprocess
from collections import defaultdict
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.count as count
import idseq_dag.util.s3 as s3
import idseq_dag.util.m8 as m8

from idseq_dag.util.dict import IdSeqDict, IdSeqDictValue, open_file_db_by_extension

MIN_ACCESSIONS_WHOLE_DB_DOWNLOAD = 5000
MAX_ACCESSION_SEQUENCE_LEN = 100000000

FASTA_SEQ_LINE_LENGTH = 80

# TODO(mark): Test that when downloading all of NT, this step works correctly.
class PipelineStepDownloadAccessions(PipelineStep):
    '''
        Download Accessions Based on the hit summary
    '''
    def run(self):
        (_align_m8, _deduped_m8, hit_summary, _orig_counts) = self.input_files_local[0]
        output_reference_fasta = self.output_files_local()[0]
        size_file = self.output_files_local()[1]
        loc_db = s3.fetch_from_s3(
            self.additional_files["loc_db"],
            self.ref_dir_local,
            allow_s3mi=True)

        (read_dict, accession_dict, _selected_genera) = m8.summarize_hits(hit_summary)

        accession_summary = defaultdict(dict)

        db_s3_path = self.additional_attributes["db"]
        db_type = self.additional_attributes["db_type"]
        lineage_db = s3.fetch_from_s3(
            self.additional_files["lineage_db"],
            self.ref_dir_local,
            allow_s3mi=True)

        if len(accession_dict) < MIN_ACCESSIONS_WHOLE_DB_DOWNLOAD:
            self.download_ref_sequences_from_s3(accession_dict, output_reference_fasta, db_type,
                                                loc_db, db_s3_path, accession_summary)
        else:
            # download the whole alignment db
            db_path = s3.fetch_from_s3(db_s3_path, self.ref_dir_local, allow_s3mi=True)
            self.download_ref_sequences_from_file(accession_dict, loc_db, db_path, output_reference_fasta, accession_summary)

        self.generate_size_file(size_file, accession_dict, loc_db, accession_summary)

    def download_ref_sequences_from_file(self, accession_dict, loc_db, db_path,
                                         output_reference_fasta, accession_summary):
        db_file = open(db_path, 'r')
        loc_dict = open_file_db_by_extension(loc_db, IdSeqDictValue.VALUE_TYPE_ARRAY)
        with open(output_reference_fasta, 'w') as orf:
            for accession, taxinfo in accession_dict.items():
                accession_data = self.get_sequence_by_accession_from_file(accession, loc_dict, db_file, accession_summary)
                if accession_data:
                    orf.write(accession_data)
        db_file.close()

    def download_ref_sequences_from_s3(self, accession_dict, output_reference_fasta, db_type,
                               loc_db, db_s3_path, accession_summary):
        ''' Download accessions specified in the selected_genera '''
        threads = []
        error_flags = {}
        semaphore = threading.Semaphore(64)
        mutex = threading.RLock()

        bucket, key = db_s3_path[5:].split("/", 1)
        loc_dict = open_file_db_by_extension(loc_db, IdSeqDictValue.VALUE_TYPE_ARRAY)
        accession_dir = os.path.join(self.output_dir_local, db_type, 'accessions')
        command.execute(f"mkdir -p {accession_dir}")
        for accession, taxinfo in accession_dict.items():
            accession_out_file = os.path.join(accession_dir, accession)
            if os.path.isfile(accession_out_file):
                continue
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepDownloadAccessions.fetch_sequence_for_thread,
                args=[
                    error_flags, accession, accession_out_file, loc_dict,
                    bucket, key, semaphore, mutex
                ])
            t.start()
            threads.append(t)
        for t in threads:
            t.join()
        if error_flags:
            raise RuntimeError("Error in getting sequences by accession list.")

        for accession, taxinfo in accession_dict.items():
            accession_out_file = os.path.join(accession_dir, accession)
            cmd = f"cat {accession_out_file} |head -1"
            seq_name = subprocess.check_output(
                    cmd, executable='/bin/bash', shell=True).decode("utf-8").split(" ", 1)[1].replace("\n", "")
            accession_summary[accession]["name"] = seq_name
        # Combine all the downloaded accessions to a fasta file
        command.execute(f"find {accession_dir}/ -type f | xargs -n 32 -P 1 cat >> {output_reference_fasta}")
    @staticmethod
    def get_sequence_by_accession_from_file(accession_id, loc_dict, db_file, accession_summary):
        entry = loc_dict.get(accession_id)
        if entry:
            range_start = int(entry[0])
            seq_len = int(entry[1]) + int(entry[2])

            db_file.seek(range_start, 0)
            seq_name = db_file.read(int(entry[1]))
            accession_summary[accession]["name"] = seq_name

            if seq_len <= MAX_ACCESSION_SEQUENCE_LEN:
                db_file.seek(range_start, 0)
                return db_file.read(seq_len)

    @staticmethod
    def fetch_sequence_for_thread(error_flags, accession, accession_out_file,
                                  loc_dict, bucket, key,
                                  semaphore, mutex):
        ''' fetch sequence from S3 for the specific accession'''
        try:
            entry = loc_dict.get(accession)
            if entry:
                range_start, name_length, seq_len = [int(e) for e in entry]
                range_end = range_start + name_length + seq_len - 1
                if seq_len <= MAX_ACCESSION_SEQUENCE_LEN:
                    num_retries = 3
                    for attempt in range(num_retries):
                        try:
                            s3.fetch_byterange(range_start, range_end, bucket, key, accession_out_file)
                            break
                        except Exception as e:
                            if attempt + 1 < num_retries:  # Exponential backoff
                                time.sleep(1.0 * (4**attempt))
                            else:
                                msg = f"All retries failed for getting sequence by accession ID {accession}: {e}"
                            raise RuntimeError(msg)

        except:
            with mutex:
                if not error_flags:
                    traceback.print_exc()
                error_flags["error"] = 1
        finally:
            semaphore.release()

    def generate_dummy_fasta_ref(self, output_reference_fasta):
        with open(output_reference_fasta, 'w') as sf:
            sf.write("test size")

    def generate_size_file(self, accession_summary_tab, accession_dict, loc_db, accession_summary):
        loc_dict = open_file_db_by_extension(loc_db, IdSeqDictValue.VALUE_TYPE_ARRAY)

        with open(accession_summary_tab, 'w') as ast:
            for accession_id, _ in accession_dict.items():
                entry = loc_dict.get(accession_id)
                if entry and len(entry) == 3:
                    # The seq_len includes newline characters in the fasta file, so we must correct for this.
                    seq_len = math.floor(FASTA_SEQ_LINE_LENGTH / (FASTA_SEQ_LINE_LENGTH + 1) * int(entry[2]))
                    # If multiple accessions map to the same sequence, they will all be listed together and separated by \x01.
                    # So we need to split here and take the first entry.
                    seq_name = accession_summary[accession_id]["name"].split('\x01')[0]

                    ast.write("%s\t%s\t%s\n" % (accession_id, seq_name, seq_len))
