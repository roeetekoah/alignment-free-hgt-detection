# ---------------- k-mer encoding (collision-free) ----------------
from typing import Set

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"  # 20 canonical AAs
AA_TO_INT = {ch: i for i, ch in enumerate(AA_ALPHABET)}
BASE = len(AA_ALPHABET)
# Fast ASCII lookup table (255 means "invalid / non-canonical")
AA_TABLE = [255] * 256
for aa, i in AA_TO_INT.items():
    AA_TABLE[ord(aa)] = i

# not used in optimized version
'''def encode_kmer_aa(kmer: str) -> int:
    """
    Encode a canonical-AA kmer into a unique integer. Returns -1 if any char not in alphabet.
    """
    v = 0
    for ch in kmer:
        i = AA_TO_INT.get(ch)
        if i is None:
            return -1
        v = v * BASE + i
    return v'''

# Un-optimized
'''def kmers_encoded_set(seq: str, k: int) -> Set[int]:
    """
    Return set of encoded k-mers. Skips kmers containing non-canonical AAs.
    """
    out: Set[int] = set()
    n = len(seq)
    for i in range(n - k + 1):
        code = encode_kmer_aa(seq[i:i+k])
        if code != -1:
            out.add(code)
    return out'''

def kmers_encoded_set(seq: str, k: int) -> Set[int]:
    """
    Return set of encoded k-mers. Skips kmers containing non-canonical AAs.

    Optimized: rolling base-BASE encoding over bytes (no slicing, no dict lookups).
    """
    n = len(seq)
    if k <= 0 or n < k:
        return set()

    # FASTA protein sequences should be ASCII; this also makes lookup cheap
    s = seq.encode("ascii", "ignore")

    out: Set[int] = set()
    out_add = out.add  # local bind is slightly faster
    tbl = AA_TABLE
    base = BASE
    high = base ** (k - 1)  # BASE^(k-1) for removing the oldest digit

    code = 0
    valid_len = 0  # number of consecutive valid chars seen in current run

    for i, b in enumerate(s):
        x = tbl[b]
        if x == 255:
            code = 0
            valid_len = 0
            continue

        if valid_len < k:
            code = code * base + x
            valid_len += 1
            if valid_len == k:
                out_add(code)
        else:
            old = tbl[s[i - k]]  # guaranteed valid when valid_len==k
            code = (code - old * high) * base + x
            out_add(code)

    return out