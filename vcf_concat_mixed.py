#!/usr/bin/env python
import sys
import os
import pysam


class ChromIter(object):

   # __slots__ = ["handle", "current", "stop"]

    def __init__(self, vcf, chrom):
        self.vcf = vcf
        self.stop = False
        self.current = None
        try:
            self.handle = vcf.fetch(chrom)
        except ValueError:  #not in Tabixfile
            self.handle = None
            self.stop = True
            return
        try:
            self.current = next(self.handle)
        except StopIteration:
            self.stop = True

    def __iter__(self):
        return self

    def next(self):
        if self.stop:
            raise StopIteration
        n = self.current
        try:
            self.current = next(self.handle)
        except StopIteration:
            self.stop = True
        return n


class MultiIter(object):
    ''' For multiple ChromIters retrieve records in order '''
    
    def __init__(self, chrom_iters):
        self.iters = chrom_iters
        self.n_recs = []
        for i in self.iters:
            self.n_recs.append( (i.next() or None) )

    def __iter__(self):
        return self

    def next(self):
        nxt = self.get_next_record()
        if nxt is None:
            raise StopIteration
        else:
            return nxt

    def get_next_record(self):
        nxt = None
        iterated = None
        for i in range(len(self.n_recs)):
            if self.n_recs[i] is not None:
                if nxt is None:
                    nxt = self.n_recs[i]
                    iterated = i
                elif self.n_recs[i].pos < nxt.pos:
                    nxt = self.n_recs[i]
                    iterated = i
                elif self.n_recs[i].pos == nxt.pos:
                    if self.n_recs[i].ref < nxt.ref:
                        nxt = self.n_recs[i]
                        iterated = i
                    elif self.n_recs[i].alt < nxt.alt:
                        nxt = self.n_recs[i]
                        iterated = i
        if iterated is not None:
            try:
                self.n_recs[iterated] = self.iters[iterated].next()
            except StopIteration:
                self.n_recs[iterated] = None
        return nxt
            
    __next__ = next #python3            
    

def get_iters(vcfs, c):
    iters = []
    for v in vcfs:
        i = ChromIter(v, c)
        if not i.stop:
            iters.append(i)
    return iters

def is_int(x):
    try:
        int(x)
        return True
    except ValueError:
        return False

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def chrom_cmp(x, y):
    if x == y:
        return 0 
    x = x.replace('chr', '')
    y = y.replace('chr', '')
    try:
        if int(x) < int(y):
            return -1
        return 1
    except ValueError:
        if is_int(x):
            return -1
        elif is_int(y):
            return 1
        else:
            for xym in ['X', 'Y', 'M', 'MT']:
                if x == xym:
                    return -1
                elif y == xym: 
                    return 1
            if x < y:
                return -1
            return 1
        
def get_contigs(fns):
    contigs = set()
    i = 0
    vcfs = []
    for v in fns:
        tbx = pysam.TabixFile(v, parser=pysam.asVCF())
        cs = tbx.contigs
        contigs.update(cs) 
        vcfs.append(tbx)
    return (sorted(contigs, key=cmp_to_key(chrom_cmp)), vcfs)

if __name__ == '__main__':
    if (len(sys.argv) < 3):                                                    
        sys.exit("Usage: {} 1.vcf.gz 2.vcf.gz [3.vcf.gz ...N.vcf.gz]" 
                 .format(sys.argv[0]))
    vf = pysam.VariantFile(sys.argv[1])
    sys.stdout.write(str(vf.header))
    vf.close()
    #does not check header matches - assumes identical for all
    contigs, tbxs = get_contigs(sys.argv[1:])
    n = 0
    prog_string = ''
    for ch in contigs:
        iters = get_iters(tbxs, ch)
        sys.stderr.write(" " * len(prog_string) + "\r")
        sys.stderr.write("For contig {} got {} VCFs.\n".format(ch, len(iters)))
        multi_iter = MultiIter(iters)
        for rec in multi_iter:
            n += 1
            print(rec)
            prog_string = "Written {:,} records\r".format(n)
            sys.stderr.write(prog_string)
    sys.stderr.write("\nFinished writing output.\n")
