# GFF3 File Format
# ================

module GFF3

using BioCore
using Indexes
using FASTX.FASTA #TODO: move responsibility to FASTX.jl.

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioCore.Exceptions: missingerror
import BioSequences
import BufferedStreams
import GenomicFeatures: GenomicFeatures, Interval, IntervalCollection
import URIParser

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
