"""
    This file contains the implementation of the Gibbs Sampler algorithm.
"""
import sys
import random
from randomized_motif_search import random_motifs, profile_from_motifs, motifs_from_profile, calculate_score, read_file, write_file

def gibbs_sampler(dna, k, N = 1000):
    """
    Finds the best motifs using the Gibbs Sampler algorithm.
    :param dna: list of DNA sequences
    :param k: length of the motif
    :param N: number of iterations
    :return: best motifs and score
    """
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    t = len(dna)
    for j in range(N):
        i = random.randint(0, t - 1)
        motifs.pop(i)
        profile = profile_from_motifs(motifs)
        motifs.insert(i, motifs_from_profile([dna[i]], profile, k)[0])
        if calculate_score(motifs) < calculate_score(best_motifs):
            best_motifs = motifs
    return best_motifs, calculate_score(best_motifs)

def main(filename):
    """ 
    Main function for the Gibbs Sampler algorithm.
    """
    dna = read_file(filename)
    best_score = int(1e9)
    best_motifs = []
    for i in range(6, min(20, len(dna[0]))):
        prev_score = best_score
        for j in range(10):
            motifs, score = gibbs_sampler(dna, i)
            if score < best_score or score - prev_score < len(dna) / 1.42:
                best_score = score
                best_motifs = motifs
    # write_file(filename, best_score, best_motifs)
    return best_score, best_motifs

if __name__ == "__main__":
    print("Running Gibbs Sampler algorithm...")
    main(sys.argv[1])
