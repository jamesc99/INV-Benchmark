#!/usr/bin/env python3
import sys, os
import pysam

def add_svlen(in_vcf):
    # determine output file names based on input name
    dirname = os.path.dirname(in_vcf)
    basename = os.path.basename(in_vcf)
    root, ext = os.path.splitext(basename)
    # output VCF with SVLEN added
    out_vcf = os.path.join(dirname, f"addLen_{root}{ext}")
    # output no-header VCF (records only)
    out_noheader = os.path.join(dirname, f"addLen_{root}_woheader{ext}")

    # open input VCF
    vfin = pysam.VariantFile(in_vcf, "r")

    # copy header and add SVLEN if missing
    hdr = vfin.header.copy()
    if "SVLEN" not in hdr.info:
        hdr.add_line(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">'
        )

    # open output VCF for writing
    vout = pysam.VariantFile(out_vcf, "w", header=hdr)
    # open plain-text file for records only
    voh = open(out_noheader, "w")

    for rec in vfin:
        # if SVLEN not already set, compute and add
        if "SVLEN" not in rec.info:
            svlen = rec.stop - rec.pos  # stop is end (0-based exclusive)
            rec.info["SVLEN"] = svlen

        vout.write(rec)
        # write record line (no header) exactly as VCF record
        voh.write(str(rec).rstrip("\n") + "\n")

    # clean up
    vfin.close()
    vout.close()
    voh.close()

    print(f"Wrote: {out_vcf}\nWrote records-only: {out_noheader}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} <input.vcf>")
    add_svlen(sys.argv[1])
