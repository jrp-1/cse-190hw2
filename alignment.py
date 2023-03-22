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
            spectra[int(line[0:-1])] = None
        else:
            spectra[int(line[0:separator])] = int(line[separator+1:-1])

    return spectra

# def peak_summer(val1, val2):


# sequence, SUM(score), modification_maxx, modification_index (0-based)
class MatchScore:
    def __init__(self, seq, mod_mass, mod_idx, spect):
        self.score = 0
        # print(spect)
        # substring from 0 to modifcaiton_index then to end
        self.sequence = seq[:mod_idx + 1] + ((lambda: "+", lambda: "")[mod_mass < 0])() + str(mod_mass) + seq[mod_idx + 1:]
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

        for peak_idx in b_vals_in_spect:
            self.score += spect[peak_idx]
        # self.score = reduce(lambda a, b: a+b, list(map(lambda a: spect[a], b_vals_in_spect)))

    def __str__(self):
        # print as 1-indexed
        return "{0} {1} {2}".format(self.sequence, self.score, self.modification_index)

# prints, void return
def q1a(spectrum_file, sequence, modification_mass):
    """q1a(spectrum_file, sequence) -> print sequence(including modification) score modification_index(1-based)
    str spectrum_file, str sequence, int modification_mass
    q1a Peptide match score using intensity of only main b-ions assuming modification on 1st amino acid"""
    spect = get_spectrum(spectrum_file)
    score = MatchScore(sequence, int(modification_mass), 0, spect)
    print(score)


# TODO: argv validations & better err messages
if len(sys.argv) == 1:
    print("No inputs", file=sys.stderr)
    sys.exit(1)

# dynamic fn call
if sys.argv[1] == "q1a":
    q1a(sys.argv[2], sys.argv[3], sys.argv[4])
elif sys.argv[1] == "q1b":
    q1b(sys.argv[2], sys.argv[3], sys.argv[4])