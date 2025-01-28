import subprocess
import numpy as np
import sys

_, model, input_file, centrifuge_db, postfix, threads = sys.argv

centrifuge_path = "/global/home/users/jieruixu/jieruixu/sediment_dna/peerj_replication/centrifuge/"

def run_centrifuge(input_file, postfix, centrifuge_db):
    centrifuge_cmd = f"{centrifuge_path}centrifuge -x {centrifuge_path}{centrifuge_db} -S {postfix}.centrifuge --report-file {postfix}.centrifugeLog -U {input_file} -p {threads}"
    print ('Run Centrifuge Classification') 
    subprocess.run(centrifuge_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    print ('Finish running centrifuge')
    kraken_report_cmd = f"{centrifuge_path}centrifuge-kreport -x {centrifuge_path}{centrifuge_db} {postfix}.centrifuge > {postfix}.k2report"
    subprocess.run(kraken_report_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    print ('Finish generating centrifuge report')

if model == 'centrifuge':
    run_centrifuge(input_file, postfix, centrifuge_db)
else:
    print ('Unsupported model')