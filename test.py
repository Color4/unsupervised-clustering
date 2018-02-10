
import pyxDamerauLevenshtein


def calc_distance(str_a, str_b):
    return normalized_damerau_levenshtein_distance(str_a, str_b)

calc_distance("AAAA", "BBBBB")



strings = {"AAAA", "AGATC", "ATTTC", "GGGAAAA"}