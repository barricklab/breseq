# Rewrites every fourth FASTQ read as a read through a mate-pair junction: its own first 22 bases,
# the 38-base Nextera junction adapter, and the first 22 bases of the read before it. Every eighth
# has one base of the adapter changed, as a sequencing error would. See testcmd.sh.
NR % 4 == 1 { name = $0 }
NR % 4 == 2 { seq = $0 }
NR % 4 == 0 {
  n++
  qual = $0
  if (n % 4 == 0 && length(seq) >= 22 && length(last) >= 22) {
    adapter = (n % 8 == 0) ? "CTGTCTCTTAAACACATCTAGATGTGTATAAGAGACAG" : "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"
    adapter_qual = adapter
    gsub(/./, "I", adapter_qual)
    print name
    print substr(seq, 1, 22) adapter substr(last, 1, 22)
    print "+"
    print substr(qual, 1, 22) adapter_qual substr(lastqual, 1, 22)
  } else {
    print name
    print seq
    print "+"
    print qual
  }
  last = seq
  lastqual = qual
}
