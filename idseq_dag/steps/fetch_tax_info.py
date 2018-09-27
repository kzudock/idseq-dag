import json
import threading
import re
import time
import wikipedia
from Bio import Entrez
import idseq_dag.util.command as command
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
from idseq_dag.engine.pipeline_step import PipelineStep


class PipelineStepFetchTaxInfo(PipelineStep):
    '''
        fetch tax info based on a list
    '''
    def run(self):
        '''
            1. fetch the taxid -> wikipedia link mappin
            2. fetch wikipedia content
            3. store everything
        '''
        taxid_list = self.input_files_local[0][0]
        (taxid2wiki, taxid2desc) = self.output_files_local()

        taxid2wikidict = {}
        Entrez.email = self.additional_attributes.get("entrez_email", "yunfang@chanzuckerberg.com")

        if s3.check_s3_presence(self.s3_path(taxid2wiki)):
            # generated
            taxid2wiki = s3.fetch_from_s3(self.s3_path(taxid2wiki), taxid2wiki)
            with open(taxid2wiki, "r") as taf:
                for line in taf:
                    (key, val) = line.rstrip("\n").split("\t")
                    taxid2wikidict[key] = val
        else:
            num_threads = self.additional_attributes.get("threads", 16)
            batch_size = self.additional_attributes.get("batch_size", 100)
            self.fetch_ncbi_wiki_map(num_threads, batch_size, taxid_list, taxid2wikidict)
            # output the data
            with open(taxid2wiki, 'w') as taxidoutf:
                for taxid, wikiurl in taxid2wikidict.items():
                    taxidoutf.write(f"{taxid}\t{wikiurl}\n")

        # output dummay for actual wiki content for now
        taxid2wikicontent = {}
        self.fetching_wiki_content(taxid2wikidict, taxid2wikicontent)

        with open(taxid2desc, 'w') as desc_outputf:
            json.dump(taxid2wikicontent, desc_outputf)

    @staticmethod
    def fetching_wiki_content(taxid2wikidict, taxid2wikicontent):
        threads = []
        semaphore = threading.Semaphore(64) # 4x the threads
        mutex = threading.RLock()
        for taxid, url in taxid2wikidict.items():
            m = re.search("curid=(\d+)", url)
            if m:
                pageid = m[1]
                semaphore.acquire()
                t = threading.Thread(
                    target=PipelineStepFetchTaxInfo.
                    get_wiki_content,
                    args=[taxid, pageid, taxid2wikicontent, mutex, semaphore]
                    )
                t.start()
                threads.append(t)
        for t in threads:
            t.join()
    @staticmethod
    def get_wiki_content(taxid, pageid, taxid2wikicontent, mutex, semaphore, max_attempt=3):
        for attempt in range(max_attempt):
            try:
                log.write(f"fetching wiki {pageid} for {taxid}")
                page = wikipedia.page(pageid=pageid)
                output = {"pageid" : page.pageid, "description": page.content[:1000]}
                with mutex:
                    taxid2wikicontent[taxid] = output
                break
            except:
                log.write(f"having trouble fetching {taxid} wiki {pageid} attempt {attempt}")
        semaphore.release()

    @staticmethod
    def fetch_ncbi_wiki_map(num_threads, batch_size, taxid_list, taxid2wikidict):
        threads = []
        semaphore = threading.Semaphore(num_threads)
        mutex = threading.RLock()
        batch = []
        with open(taxid_list, 'r') as taxf:
            for line in taxf:
                taxid = line.rstrip()
                if taxid == 'taxid':
                    continue # header
                batch.append(taxid)
                if len(batch) >= batch_size:
                    semaphore.acquire()
                    t = threading.Thread(
                        target=PipelineStepFetchTaxInfo.
                        get_taxid_mapping_for_batch,
                        args=[batch, taxid2wikidict, mutex, semaphore]
                        )
                    t.start()
                    threads.append(t)
                    batch = []
        if len(batch) > 0:
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepFetchTaxInfo.
                get_taxid_mapping_for_batch,
                args=[batch, taxid2wikidict, mutex, semaphore]
                )
            t.start()
            threads.append(t)
        for t in threads:
            t.join()


    @staticmethod
    def get_taxid_mapping_for_batch(taxids, taxid2wikidict, mutex, semaphore, max_attempt=3):
        taxid_str = ",".join(taxids)
        log.write(f"fetching batch {taxid_str}")
        for attempt in range(max_attempt):
            try:
                handle = Entrez.elink(dbfrom="taxonomy", id=taxid_str, cmd="llinks")
                record = Entrez.read(handle)
                handle.close()

                parsed = {}
                results = record[0]['IdUrlList']['IdUrlSet']
                for result in results:
                    taxid = result['Id']
                    wikiurl = ""
                    for link in result['ObjUrl']:
                        url = str(link['Url'])
                        if re.search('wikipedia.org', url):
                            wikiurl = url
                            break
                    parsed[taxid] = wikiurl
                break
            except:
                log.write(f"failed batch attempt {attempt}")
                time.sleep(5)
        semaphore.release()
        with mutex:
            taxid2wikidict.update(parsed)


    def count_reads(self):
        pass

