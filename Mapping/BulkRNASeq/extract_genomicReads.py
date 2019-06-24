def getNonGenomicRead(infile):
    import pysam
    bamfile =  pysam.AlignmentFile(infile,'rb')
    #fout = open(outfile,"w")
    for read in bamfile.fetch():
	    if not read.get_tag("XF") in ["INTERGENIC","INTRONIC"] and not read.is_secondary :
                print "%s" % read.qname
    #fout.close()



if __name__ == "__main__":

    import sys
    from sys import argv

    if len(argv) ==0 or len(argv) > 2:
        print "Usage:"
        print "python extract_genomicReads.py input.bam out.bed"
        print "Error: No argument was passed"
        sys.exit(1)

    infile = argv[1]
    #outfile = argv[2]

    getNonGenomicRead(infile)
    #print "%s generated" % outfile
