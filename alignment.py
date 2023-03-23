#!/usr/bin/env python3
"""Author: Janne Rapakko (A08240805)
denovo.py
CSE 190 (Bandeira) - Homework 1"""
import sys # cmd line args

# from functools import reduce


AA_VALUES = {
    "A": 71,
    "R": 156,
    "N": 114,
    "D": 115,
    "C": 103,
    "E": 129,
    "Q": 128,
    "G": 57,
    "H": 137,
    "I": 113,
    "L": 113,
    "K": 128,
    "M": 131,
    "F": 147,
    "P": 97,
    "S": 87,
    "T": 101,
    "W": 186,
    "Y": 163,
    "V": 99
}

# amino acid weight in Da
def getvalue(amino):
    """getvalue(amino) -> int"""
    assert len(amino) == 1
    return AA_VALUES[str(amino)]

# return dict of spectrum
def get_spectrum(spectrum_file):
    """get_spectrum(spectrum_file) -> dict """
    with open(spectrum_file, mode="r", encoding="utf-8") as file:
        spectrum_list = file.readlines()

    spectra = {}
    for line in spectrum_list:
        separator = line.find('\t')
        if separator == -1:
            spectra[int(line[0:])] = 0
        else:
            spectra[int(line[0:separator])] = int(line[separator+1:])

    return spectra

def mod_list(mass):
    mlist = list(map(lambda a: mass + a, range(2,50)))
    return(mlist)


# sequence, SUM(score), modification_maxx, modification_index (0-based)
class MatchScore:
    def __init__(self, seq, mod_mass, mod_idx, spect):
        self.score = 0
        # print(spect)
        # substring from 0 to modifcaiton_index then to end
        if mod_mass != 0:
            self.sequence = seq[:mod_idx + 1] + ((lambda: "+", lambda: "")[mod_mass < 0])() + str(mod_mass) + seq[mod_idx + 1:]
        else:
            self.sequence = seq
        self.modification_mass = mod_mass
        self.modification_index = mod_idx + 1 # 1 indexed for printing

        sequence_values = list(map(lambda aa: getvalue(aa), list(seq)))
        sequence_values[mod_idx] += mod_mass

        b_vals = list()
        for i in range(0, len(sequence_values)):
            if i == 0:
                #base case -- + 1 for H-ion
                b_vals.append(sequence_values[i] + 1)
            else:
                b_vals.append(b_vals[i-1] + sequence_values[i])

        # peaks =  list(map(lambda val: ((lambda: spect[val], lambda: 0)[val in spect.keys()])(), b_vals))
        b_vals_in_spect = list(filter(lambda a: a in spect.keys(), b_vals))

        # print(b_vals_in_spect)

        for peak_idx in b_vals_in_spect:
            self.score += spect[peak_idx]
            # print(self.score)

        # self.score = reduce(lambda a, b: a+b, list(map(lambda a: spect[a], b_vals_in_spect)))

    def __str__(self):
        # print as 1-indexed
        return "{0} {1} {2}".format(self.sequence, self.score, self.modification_index)

    def __lt__(self, other):
        return self.score < other.score

    def get_score(self):
        return self.score

    def get_sequence(self):
        return self.sequence

# prints, void return
def q1a(spectrum_file, sequence, modification_mass):
    """q1a(spectrum_file, sequence) -> print sequence(including modification) score modification_index(1-based)
    str spectrum_file, str sequence, int modification_mass
    q1a Peptide match score using intensity of only main b-ions assuming modification on 1st amino acid"""
    spect = get_spectrum(spectrum_file)
    score = MatchScore(sequence, int(modification_mass), 0, spect)
    print(score)
    sys.exit(0)

# prints, void return
def q1b(spectrum_file, sequence, modification_mass):
    """q1a(spectrum_file, sequence) -> print sequence(including modification) score modification_index(1-based)
    str spectrum_file, str sequence, int modification_mass
    q1a Peptide match score using intensity of only main b-ions with highest value for 1 modification"""
    scores = list()
    spect = get_spectrum(spectrum_file)
    for i in range(0, len(sequence)):
        scores.append(MatchScore(sequence, int(modification_mass), i, spect))
    print(max(scores, key=lambda x: x.get_score()))
    sys.exit(0)

def q2(spectrum_file, sequence):
    # builder -- build each with spectrum matches -- remove if next can't match
    spect = get_spectrum(spectrum_file)
    # print(MatchScore(sequence, 0, 0, spect))
    scores = list()
    scores.append(MatchScore(sequence, 0, 0, spect)) # no modifications

    # filter before creating new MatchScores for each possible modification
    # seq_list = {aa : getvalue(aa) for aa in list(sequence)}
    seq_list = list(map(lambda x: getvalue(x), list(sequence)))

    # seq_list[sequence[0]] += 1 # 1 Da for H-ion
    seq_list[0] += 1

    tmp_score = 0
    # use mod_list once modified
    for possible_mod in seq_list:
        # mlist = list(filter(lambda x: x + tmp_score in spect.keys(), mod_list(seq_list[possible_mod])))
        mlist = list(filter(lambda x: x + tmp_score in spect.keys(), mod_list(possible_mod)))
        # print(list(map(lambda x: x + tmp_score, mlist)))
        # tmp_score += seq_list[possible_mod]
        tmp_score += possible_mod
        for modification in mlist:  # if modification mass is in spectrum
            # TODO: MISSING SuMMATION
            # scores.append(MatchScore(sequence, modification - seq_list[possible_mod], seq_list.index(possible_mod), spect))
            scores.append(MatchScore(sequence, modification - possible_mod, seq_list.index(possible_mod), spect))
            # print(MatchScore(sequence, modification - possible_mod, seq_list.index(possible_mod), spect))

    # print(seq_list)
    # print(mlist)
    # print(spect)
    # print(scores)
    print(max(scores, key=lambda x: x.get_score()))

    # score => sequence (1 mod max) , score collisions
    # check modifications, need to add new sequence for each (prune not in spectrum) need to add flag for mod?

    sys.exit(0)

def ec(spectrum_file, sequence):
    spect = get_spectrum(spectrum_file)
    scores = list()

    all_seqt = list()

    seq_len = len(sequence)
    for i in range(seq_len):
        tmp = ""
        for j in range(i, seq_len):
            tmp += sequence[j]
            all_seqt.append(tmp)

    spect_len = len(spect)

    all_seq = list(filter(lambda x: len(x) <= spect_len, all_seqt))
    # print(all_seq)


    # substrings of sequence -- then

    scores.append(MatchScore(sequence, 0, 0, spect)) # no modifications

    # filter before creating new MatchScores for each possible modification
    # seq_list = {aa : getvalue(aa) for aa in list(sequence)}
    # seq_list = list(map(lambda x: getvalue(x), list(sequence)))
    all_seq_list = list()
    for seq in all_seq:
        all_seq_list.append(list(map(lambda x: getvalue(x), seq)))

    # all_seq_list = list(map(lambda x: getvalue(x), all_seq))

    # print(all_seq_list)

    # print(seq_list)

    # print(spect)

    # seq_list[sequence[0]] += 1 # 1 Da for H-ion
    # seq_list[0] += 1
    # +1 Da for H-ion
    # for seq in all_seq_list
    #     seq[0] += 1

    # tmp_score = 0
    # use mod_list once modified
    for possible_seq in all_seq_list:
        tmp_score = 0
        seq_idx = all_seq_list.index(possible_seq)
        possible_seq[0] += 1 # 1 Da for H-ion for first AA
        for possible_mod in possible_seq:
            # mlist = list(filter(lambda x: x + tmp_score in spect.keys(), mod_list(seq_list[possible_mod])))
            mlist = list(filter(lambda x: x + tmp_score in spect.keys(), mod_list(possible_mod)))
            # print(list(map(lambda x: x + tmp_score, mlist)))
            # tmp_score += seq_list[possible_mod]
            tmp_score += possible_mod
            for modification in mlist:  # if modification mass is in spectrum
                # TODO: MISSING SuMMATION
                # scores.append(MatchScore(sequence, modification - seq_list[possible_mod], seq_list.index(possible_mod), spect))
                scores.append(MatchScore(all_seq[seq_idx], modification - possible_mod, possible_seq.index(possible_mod), spect))
                # print(modification, all_seq[seq_idx], possible_mod, possible_seq)
            # print(MatchScore(sequence, modification - possible_mod, seq_list.index(possible_mod), spect))

    # print(seq_list)
    # print(mlist)
    # print(spect)
    # print(scores)
    # for score in scores:
    #     print(score)
    print(max(scores, key=lambda x: x.get_score()))
    sys.exit(0)


# TODO: argv validations & better err messages
if len(sys.argv) == 1:
    print("No inputs", file=sys.stderr)
    sys.exit(1)

# dynamic fn call
if sys.argv[1] == "q1a":
    q1a(sys.argv[2], sys.argv[3], sys.argv[4])
elif sys.argv[1] == "q1b":
    q1b(sys.argv[2], sys.argv[3], sys.argv[4])
elif sys.argv[1] == "q2a":
    q2(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "q2b":
    q2(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "ec":
    ec(sys.argv[2], sys.argv[3])
else:
    print("Command not found", file=sys.stderr)
    sys.exit(2)