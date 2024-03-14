# python .\main.py .\MA0003.2
# python .\main.py .\MA0016.1
# python .\main.py .\MA0024.2

from randomized_motif_search import main as rms_main, consensus_from_motifs, consensus_from_pfm
from gibbs_sampler import main as gs_main

import sys
import time

def main():
    print("Running randomized_motif_search.py...")
    start_time = time.time()
    rms_score, rms_motifs = rms_main(sys.argv[1])
    end_time = time.time()
    rms_runtime = end_time - start_time

    print("Running Gibbs Sampler algorithm...")
    start_time = time.time()
    gs_score, gs_motifs = gs_main(sys.argv[1])
    end_time = time.time()
    gs_runtime = end_time - start_time

    true_motif = consensus_from_pfm(sys.argv[1])
    with open('output.txt', 'w', encoding="utf-8") as file:
        file.write("true motif:\t" + true_motif + '\n')
        file.write("RMS motif:\t" + consensus_from_motifs(rms_motifs) + '\n')
        file.write("GS motif:\t" + consensus_from_motifs(gs_motifs) + '\n')
        file.write("RMS score:\t" + str(rms_score) + '\n')
        file.write("GS score:\t" + str(gs_score) + '\n')
        file.write("RMS runtime:\t" + str(rms_runtime) + '\n')
        file.write("GS runtime:\t" + str(gs_runtime) + '\n')

main()


