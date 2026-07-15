# GFF3 File Format
# ================

module GFF3

using BioGenerics
using Indexes
using FASTX.FASTA #TODO: move responsibility to FASTX.jl.
using TranscodingStreams

import Automa
import Automa: @re_str, onenter!, onexit!, rep, rep1, opt

import BGZFStreams
import BioGenerics.Exceptions: missingerror
import GenomicFeatures: GenomicFeatures, GenomicInterval, GenomicIntervalCollection
import URIParser

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
