GFF3.jl
====

Description
-----------

GFF3 is a text-based file format for representing genomic annotations. The major
difference from BED is that GFF3 is more structured and can include sequences in
the FASTA file format.

I/O tools for GFF3 are provided by the `GFF3` module,
which exports following three types:
* Reader type: `GFF3.Reader`
* Writer type: `GFF3.Writer`
* Element type: `GFF3.Record`

A GFF3 file may contain directives and/or comments in addition to genomic
features. These lines are skipped by default but you can control the behavior by
passing keyword arguments to `GFF3.Reader`. See the docstring for details.

Installation
------------

Install GFF3 from the Julia REPL:

```julia
julia> Pkg.add("GFF3")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

Examples
--------

Here is a common workflow to iterate over all records in a GFF3 file:
```julia
# Import the GFF3 module.
using GFF3

# Open a GFF3 file.
reader = open(GFF3.Reader, "data.gff3")

# Iterate over records.
for record in reader
    # Do something on record (see Accessors section).
    seqid = GFF3.seqid(reader)
    # ...
end

# Finally, close the reader
close(reader)
```

If you are interested in directives (which starts with '#') in addition to
genomic features, you need to pass `skip_directives=false` when initializing a
GFF3 constructor:
```julia
# Set skip_directives to true (this is set to false by default).
reader = GFF3.Reader(open("data.gff3"), skip_directives=false)
for record in record
    # Branch by record type.
    if GFF3.isfeature(record)
        # ...
    elseif GFF3.isdirective(record)
        # ...
    else
        # This never happens.
        assert(false)
    end
end
close(reader)
```

jl supports [tabix](http://www.htslib.org/doc/tabix.html) to
retrieve records overlapping with a specific interval. First you need to create
a block compression file from a GFF3 file using bgzip and then index it using
tabix.
```
cat data.gff3 | grep -v "^#" | sort -k1,1 -k4,4n | bgzip >data.gff3.bgz
tabix data.gff3.bgz  # this creates data.gff3.bgz.tbi
```

Then you can read the block compression file as follows:
```julia
# Read the block compression gzip file.
reader = GFF3.Reader("data.gff3.bgz")
for record in eachoverlap(reader, Interval("chr1", 250_000, 300_000))
    # Each record overlap the query interval.
    # ...
end
```


API
---

```@docs
GFF3.Reader
GFF3.directives
GFF3.hasfasta
GFF3.getfasta
GFF3.Writer
GFF3.Record
GFF3.isfeature
GFF3.isdirective
GFF3.iscomment
GFF3.seqid
GFF3.source
GFF3.featuretype
GFF3.seqstart
GFF3.seqend
GFF3.score
GFF3.strand
GFF3.phase
GFF3.attributes
GFF3.content
```
