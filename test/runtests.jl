using GFF3
using Test
using Documenter

using FASTX.FASTA
using FormatSpecimens
using GenomicFeatures

import BioSequences: @dna_str

import BGZFStreams

import BioCore.Exceptions: MissingFieldException

import BioCore:
    hasleftposition,
    leftposition,
    hasrightposition,
    rightposition,
    hasseqname,
    isfilled,
    seqname



@testset "GFF3" begin
    record = GFF3.Record()
    @test !isfilled(record)
    @test repr(record) == "GFF3.Record: <not filled>"

    record = GFF3.Record("CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tNAME=LINC00115")
    @test isfilled(record)
    @test GFF3.isfeature(record)
    @test hasseqname(record)
    @test GFF3.hasseqid(record)
    @test seqname(record) == GFF3.seqid(record) == "CCDS1.1"
    @test GFF3.hassource(record)
    @test GFF3.source(record) == "CCDS"
    @test GFF3.hasfeaturetype(record)
    @test GFF3.featuretype(record) == "gene"
    @test GFF3.hasseqstart(record) === hasleftposition(record) === true
    @test GFF3.seqstart(record) === leftposition(record) === 801943
    @test GFF3.hasseqend(record) === hasrightposition(record) === true
    @test GFF3.seqend(record) === rightposition(record) === 802434
    @test !GFF3.hasscore(record)
    @test_throws MissingFieldException GFF3.score(record)
    @test GFF3.hasstrand(record)
    @test strand(record) === GFF3.strand(record) === STRAND_NEG
    @test !GFF3.hasphase(record)
    @test_throws MissingFieldException GFF3.phase(record)
    @test GFF3.attributes(record) == ["NAME" => ["LINC00115"]]
    @test GFF3.attributes(record, "NAME") == ["LINC00115"]
    @test GFF3.content(record) == "CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tNAME=LINC00115"
    @test startswith(repr(record), "GFF3.Record:\n")
    @test string(record) == "CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tNAME=LINC00115"

    record = GFF3.Record("##gff-version 3")
    @test isfilled(record)
    @test GFF3.isdirective(record)
    @test GFF3.content(record) == "gff-version 3"
    @test convert(String, record) == "##gff-version 3"

    record = GFF3.Record("#comment")
    @test isfilled(record)
    @test GFF3.iscomment(record)
    @test GFF3.content(record) == "comment"
    @test convert(String, record) == "#comment"

    function check_gff3_parse(filename)
        # Reading from a stream
        num_intervals = 0
        for interval in GFF3.Reader(open(filename))
            num_intervals += 1
        end

        # Reading from a regular file
        for interval in open(GFF3.Reader, filename)
        end

        collection = IntervalCollection(open(GFF3.Reader, filename))
        @test length(collection) == num_intervals

        # in-place parsing
        stream = open(GFF3.Reader, filename)
        entry = eltype(stream)()
        while !eof(stream)
            try
                read!(stream, entry)
            catch ex
                if isa(ex, EOFError)
                    break
                end
            end
        end
        close(stream)

        # copy
        records = GFF3.Record[]
        reader = open(GFF3.Reader, filename)
        output = IOBuffer()
        writer = GFF3.Writer(output)
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        close(reader)
        flush(writer)

        records2 = GFF3.Record[]
        for record in GFF3.Reader(IOBuffer(take!(output)))
            push!(records2, record)
        end
        return records == records2
    end

    dir_gff3 = path_of_format("GFF3")

    for specimen in list_valid_specimens("GFF3")

        if hastag(specimen, "gzip")
            # skip compressed files
            continue
        end

        filepath = joinpath(dir_gff3, filename(specimen))

        @test check_gff3_parse(filepath)

    end

    for specimen in list_invalid_specimens("GFF3")

        if "gzip" ∈ split(get(specimen, "tags", ""))
            # skip compressed files
            continue
        end

        filepath = joinpath(dir_gff3, filename(specimen))

        @test_throws Exception check_gff3_parse(filepath)
    end

    # no fasta
    test_input = """
1	havana	exon	870086	870201	.	-	.	Parent=transcript:ENST00000432963;Name=ENSE00001791782;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001791782;rank=1;version=1
1	havana	lincRNA	868403	876802	.	-	.	ID=transcript:ENST00000427857;Parent=gene:ENSG00000230368;Name=FAM41C-002;biotype=lincRNA;havana_transcript=OTTHUMT00000007022;havana_version=1;transcript_id=ENST00000427857;transcript_support_level=3;version=1
"""
    stream = GFF3.Reader(IOBuffer(test_input))
    collect(stream)
    @test !GFF3.hasfasta(stream)
    @test_throws Exception GFF3.getfasta(stream)


    # implicit fasta
    test_input2 = string(test_input, """
>seq1
ACGTACGT
>seq2
TGCATGCA
""")
    stream = GFF3.Reader(IOBuffer(test_input2))
    collect(stream)
    @test GFF3.hasfasta(stream)
    @test collect(GFF3.getfasta(stream)) ==
        [FASTA.Record("seq1", dna"ACGTACGT")
         FASTA.Record("seq2", dna"TGCATGCA")]

    # explicit fasta
    test_input3 = string(test_input, """
##FASTA
>seq1
ACGTACGT
>seq2
TGCATGCA
""")
    stream = GFF3.Reader(IOBuffer(test_input3))
    collect(stream)
    @test GFF3.hasfasta(stream)
    @test collect(GFF3.getfasta(stream)) ==
        [FASTA.Record("seq1", dna"ACGTACGT")
         FASTA.Record("seq2", dna"TGCATGCA")]


    test_input4 = """
##directive1
#comment1
##directive2
1	havana	exon	869528	869575	.	-	.	Parent=transcript:ENST00000432963;Name=ENSE00001605362;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001605362;rank=2;version=1
1	havana	exon	870086	870201	.	-	.	Parent=transcript:ENST00000432963;Name=ENSE00001791782;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001791782;rank=1;version=1
##directive3
#comment2
##directive4
1	havana	lincRNA	868403	876802	.	-	.	ID=transcript:ENST00000427857;Parent=gene:ENSG00000230368;Name=FAM41C-002;biotype=lincRNA;havana_transcript=OTTHUMT00000007022;havana_version=1;transcript_id=ENST00000427857;transcript_support_level=3;version=1
##directive5
#comment3
##directive6
"""
    stream = GFF3.Reader(IOBuffer(test_input4), save_directives=true)
    read(stream)
    @test GFF3.directives(stream) == ["directive1", "directive2"]
    read(stream)
    @test isempty(GFF3.directives(stream))
    read(stream)
    @test GFF3.directives(stream) == ["directive3", "directive4"]
    @test_throws EOFError read(stream)
    @test eof(stream)
    close(stream)
    @test GFF3.directives(stream) == ["directive5", "directive6"]

    test_input5 = """
    ##directive1
    feature1\t.\t.\t.\t.\t.\t.\t.\t
    #comment1
    feature2\t.\t.\t.\t.\t.\t.\t.\t
    ##directive2
    feature3\t.\t.\t.\t.\t.\t.\t.\t
    """
    @test [r.kind for r in GFF3.Reader(IOBuffer(test_input5))] == [:feature, :feature, :feature]
    @test [r.kind for r in GFF3.Reader(IOBuffer(test_input5), skip_directives=false)] == [:directive, :feature, :feature, :directive, :feature]
    @test [r.kind for r in GFF3.Reader(IOBuffer(test_input5), skip_directives=false, skip_comments=false)] == [:directive, :feature, :comment, :feature, :directive, :feature]

    @testset "eachoverlap" begin
        path = joinpath(path_of_format("GFF3"), "TAIR10.part.gff.bgz")
        stream = BGZFStreams.BGZFStream(path)
        reader = GFF3.Reader(stream, index=string(path, ".tbi"))
        for (interval, n_records) in [
                (Interval("ChrC", 1:70_000), 68),
                (Interval("Chr4", 1:10_000),  2),
                (Interval("ChrM", 1:30_000), 10),]
            n = 0
            for record in eachoverlap(reader, interval)
                n += 1
            end
            @test n == n_records
        end
        @test isa(GFF3.Reader(path), GFF3.Reader)
    end

    # Include doctests.
    DocMeta.setdocmeta!(GFF3, :DocTestSetup, :(using GFF3); recursive=true)
    doctest(GFF3; manual = false)
end
