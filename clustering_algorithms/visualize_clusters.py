from Bio.SeqUtils.ProtParam import ProteinAnalysis

seq = ProteinAnalysis("KDENDHNKDENDHNKDENDHNKDENDHN")
print(seq.molecular_weight())

