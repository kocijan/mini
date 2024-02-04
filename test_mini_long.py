from minipy import KMer, minimize

TEST_SEQ =  "AACCTTGGACTACGATCGGGGGAACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC" \
            "AACCTTGGACTACGATCGGGGGAACCCCGAACATCTCCTCTCCCATTCTCCCTCCCCTAGAGATTCATTC"

KMER_LEN = 58
WIN_LEN = 9

MINMIZERS = [
    KMer(916718422212990171, 6,  False),
    KMer(4626220760283687045, 11, False),
    KMer(2591599766356106169, 17, False),
    KMer(176679790298658678, 24, False),
    KMer(185623040052900602, 27, False),
    KMer(268886174752216349, 35, False),
    KMer(2750856359964708494, 38, False),
    KMer(1152133110661816871, 47, False),
    KMer(1255617870259713102, 53, False),
    KMer(749452792572590264, 57, False),
    KMer(144082624022435439, 60, False),
    KMer(1325561104408984909, 62, False),
    KMer(2434044891487640860, 69, False),
    KMer(916718422212990171, 76, False)
]

def test_minmize():
  minimizers = minimize(TEST_SEQ, KMER_LEN, WIN_LEN)
  
  assert len(minimizers) == len(MINMIZERS)
  for v, e in zip(minimizers, MINMIZERS):
    assert v.value() == e.value()
    assert v.position() == e.position()
    assert v.strand() == e.strand()

test_minmize()