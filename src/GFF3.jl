# GFF3 File Format
# ================

module GFF3

using BioGenerics
using Indexes
using FASTX.FASTA #TODO: move responsibility to FASTX.jl.
using TranscodingStreams

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos

import BGZFStreams
import BioGenerics.Exceptions: missingerror
import GenomicFeatures: GenomicFeatures, GenomicInterval, GenomicIntervalCollection
import URIParser

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
