''' Build Custom Blast Index '''
import os
import urllib.request
from idseq_dag.engine.pipeline_step import PipelineStep
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
import idseq_dag.util.count as count

class PipelineStepBuildCustomBlastIndex(PipelineStep):
    ''' Build Custom Blast Index implementation '''
    def run(self):
        """
          Build custom blast index from an S3 location or a url
        """
        _input_files = self.input_files_local[0] # dummy in this case
        output_tar_file = self.output_files_local()[0]

        db_type = self.additional_attributes['db_type']
        file_source = self.additional_attributes['data_source']
        output_db_name = output_tar_file.replace(".tar", "")

        if file_source.startswith("s3://"):
            db_file = s3.fetch_from_s3(file_source, self.output_dir_local)
        else:
            # Download with wget
            db_file = os.path.join(self.output_dir_local, os.path.basename(file_source))
            urllib.request.urlretrieve(file_source, db_file)
            self.additional_files_to_upload.append(db_file)

        # Build blast index
        if db_file.endswith(".bz2"):
            command.execute(f"bzip2 -dk {db_file}")
            db_file = db_file[:-4]
        elif db_file.endswith(".zip"):
            command.execute(f"unzip {db_file}")
            db_file = db_file[:-4]

        command.execute(f"makeblastdb -in {db_file} -dbtype {db_type} -out {output_db_name}")
        command.execute(f"tar cvf {output_tar_file}  {output_db_name}.*")

    def count_reads(self):
        ''' Count reads '''
        pass


