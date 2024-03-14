"""
    This file contains the code for the Randomized Motif Search algorithm.
"""
import sys
import random

def random_motifs(dna, k):
    """
    Randomly selects k-mers from each DNA sequence.
    :param dna: list of DNA sequences
    :param k: length of the motif
    :return: list of k-mers
    """
    t = len(dna)
    motifs = []
    for i in range(t):
        n = len(dna[i])
        start = random.randint(0, n-k)
        motifs.append(dna[i][start:start+k])
    return motifs

def profile_from_motifs(motifs):
    """
    Creates a profile matrix from the motifs.
    :param motifs: list of k-mers
    :return: profile matrix
    """
    k = len(motifs[0])
    profile = {'a': [1]*k, 'c': [1]*k, 'g': [1]*k, 't': [1]*k}
    for i in range(k):
        for motif in motifs:
            profile[motif[i]][i] += 1
    for key in profile:
        for i in range(k):
            profile[key][i] /= (len(motifs) + 4)
    return profile

def motifs_from_profile(dna, profile, k):
    """
    Creates motifs from the profile matrix.
    :param dna: list of DNA sequences
    :param profile: profile matrix
    :param k: length of the motif
    :return: list of motifs
    """
    motifs = []
    for seq in dna:
        n = len(seq)
        max_prob = -1
        kmer = ''
        for i in range(n-k+1):
            prob = 1
            for j in range(k):
                prob *= profile[seq[i+j]][j]
            if prob > max_prob:
                max_prob = prob
                kmer = seq[i:i+k]
        motifs.append(kmer)
    return motifs

def calculate_score(motifs):
    """
    Calculates the score of the motifs.
    :param motifs: list of k-mers
    :return: score
    """
    k = len(motifs[0])
    t = len(motifs)
    score = 0
    for i in range(k):
        count = {'a': 0, 'c': 0, 'g': 0, 't': 0}
        for motif in motifs:
            count[motif[i]] += 1
        score += t - max(count.values())
    return score

def randomized_motif_search(dna, k):
    """
    Randomized Motif Search algorithm.
    :param dna: list of DNA sequences
    :param k: length of the motif
    :return: list of best motifs
    """
    motifs = random_motifs(dna, k)
    best_motifs = motifs
    while True:
        profile = profile_from_motifs(motifs)
        motifs = motifs_from_profile(dna, profile, k)
        if calculate_score(motifs) < calculate_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs, calculate_score(best_motifs)

def read_file(filename):
    """ Read the input file. """
    filename = filename + '.sites'
    skip = False
    with open(filename, 'r', encoding="utf-8") as file:
        # read line by line
        lines = file.readlines()
        dna = []
        for i, line in enumerate(lines):
            if skip:
                skip = False
                continue
            if line[0] == '>':
                dna.append(lines[i+1].strip().lower())
                skip = True
            else:
                dna[-1] += line.strip().lower()
    return dna

def consensus_from_pfm(filename):
    """ Create a consensus sequence from the pfm file. """
    filename = filename + '.pfm'
    with open(filename, 'r', encoding="utf-8") as file:
        lines = file.readlines()
        # make a 2d matrix from the file
        pfm = []
        for i in range(1, len(lines)):
            pfm.append(list(map(int, list(map(float, lines[i].split())))))
    consensus = ''
    for i in range(len(pfm[0])):
        max_freq = 0
        max_base = ''
        for j in range(4):
            if pfm[j][i] > max_freq:
                max_freq = pfm[j][i]
                max_base = 'acgt'[j]
        consensus += max_base
    return consensus


def consensus_from_motifs(motifs):
    """
    Creates a consensus sequence from the motifs.
    :param motifs: list of k-mers
    :return: consensus sequence
    """
    k = len(motifs[0])
    count = []
    for i in range(k):
        count.append({'a': 0, 'c': 0, 'g': 0, 't': 0})
        for motif in motifs:
            count[i][motif[i]] += 1
    consensus = ''
    for i in range(k):
        consensus += max(count[i], key=count[i].get)
    return consensus
    
def write_file(filename, best_score, best_motifs):
    """ Write the output to the file. """
    true_motif = consensus_from_pfm(filename)
    our_motif = consensus_from_motifs(best_motifs)
    with open('output.txt', 'w', encoding="utf-8") as file:
        file.write("true motif:\t" + true_motif + '\n')
        file.write("our motif:\t" + our_motif + '\n')
        file.write("score: " + str(best_score) + '\n')
        for motif in best_motifs:
            file.write(motif + '\n')

def main(filename):
    """ Main function. """
    dna = read_file(filename)
    best_score = int(1e9)
    best_motifs = []
    for i in range(6, min(20, len(dna[0]))):
        prev_score = best_score
        for j in range(10):
            motifs, score = randomized_motif_search(dna, i)
            if score < best_score or score - prev_score < len(dna) / 1.42:
                best_score = score
                best_motifs = motifs
    # write_file(filename, best_score, best_motifs)
    return best_score, best_motifs

if __name__ == "__main__":
    print("Running randomized_motif_search.py...")
    main(sys.argv[1])
