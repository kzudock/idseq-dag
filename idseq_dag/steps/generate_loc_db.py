''' Generate loc db  '''
import dbm
import shelve
import re
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.command import run_in_subprocess
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.count as count
from idseq_dag.util.dict import IdSeqDict, IdSeqDictValue
BATCH_INSERT_SIZE = 300

class PipelineStepGenerateLocDB(PipelineStep):
    ''' Generate Loc DB for NT/NR '''
    def run(self):
        """
        Generate Loc DB for NT/NR
        """
        db_file = self.input_files_local[0][0]
        loc_db_file = self.output_files_local()[0]
        self.generate_loc_db(db_file, loc_db_file)

    @run_in_subprocess
    def generate_loc_db(self, db_file, loc_db_file):
        loc_db = IdSeqDict(loc_db_file, IdSeqDictValue.VALUE_TYPE_ARRAY)
        batch_list = []
        with open(db_file) as dbf:
            seq_offset = 0
            seq_len = 0
            header_len = 0
            lines = 0
            accession_id = ""
            for line in dbf:
                lines += 1
                if lines % 100000 == 0:
                    log.write(f"{lines/1000000.0}M lines")
                if line[0] == '>':  # header line
                    if seq_len > 0 and len(accession_id) > 0:
                        batch_list.append((accession_id, [seq_offset, header_len, seq_len]))
                        if len(batch_list) >= BATCH_INSERT_SIZE:
                            loc_db.batch_inserts(batch_list)
                            batch_list = []
                    seq_offset = seq_offset + header_len + seq_len
                    header_len = len(line)
                    seq_len = 0
                    s = re.match('^>([^ ]*).*', line)
                    if s:
                        accession_id = s.group(1)
                else:
                    seq_len += len(line)
            if seq_len > 0 and len(accession_id) > 0:
                batch_list.append((accession_id, [seq_offset, header_len, seq_len]))
            loc_db.batch_inserts(batch_list)

    @run_in_subprocess
    def generate_loc_db_old(self, db_file, loc_db_file):
        # TODO: To be deprecated. Using shelve
        loc_dict = shelve.Shelf(dbm.ndbm.open(loc_db_file.replace(".db", ""), 'c'))
        with open(db_file) as dbf:
            seq_offset = 0
            seq_len = 0
            header_len = 0
            lines = 0
            accession_id = ""
            for line in dbf:
                lines += 1
                if lines % 100000 == 0:
                    log.write(f"{lines/1000000.0}M lines")
                if line[0] == '>':  # header line
                    if seq_len > 0 and len(accession_id) > 0:
                        loc_dict[accession_id] = (seq_offset, header_len, seq_len)
                    seq_offset = seq_offset + header_len + seq_len
                    header_len = len(line)
                    seq_len = 0
                    s = re.match('^>([^ ]*).*', line)
                    if s:
                        accession_id = s.group(1)
                else:
                    seq_len += len(line)
            if seq_len > 0 and len(accession_id) > 0:
                loc_dict[accession_id] = (seq_offset, header_len, seq_len)
        loc_dict.close()

    def count_reads(self):
        ''' Count reads '''
        pass

