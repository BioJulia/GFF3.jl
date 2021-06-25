# GFF3 Reader
# ===========

mutable struct Reader{S <: TranscodingStream} <: BioGenerics.IO.AbstractReader
    state::BioGenerics.Automa.State{S}
    index::Union{Indexes.Tabix, Nothing}
    save_directives::Bool
    targets::Vector{Symbol}
    found_fasta::Bool
    directives::Vector{Record}
    directive_count::Int
    preceding_directive_count::Int

    function Reader(input::S, index=nothing, save_directives::Bool=false, skip_features::Bool=false, skip_directives::Bool=true, skip_comments::Bool=true) where S <: TranscodingStream

        if isa(index, Indexes.Tabix) && !isa(input.stream, BGZFStreams.BGZFStream)
            throw(ArgumentError("not a BGZF stream"))
        end
        targets = Symbol[]
        if !skip_features
            push!(targets, :feature)
        end
        if !skip_directives
            push!(targets, :directive)
        end
        if !skip_comments
            push!(targets, :comment)
        end
        return new{S}(BioGenerics.Automa.State(input, body_machine.start_state, 1, false), index, save_directives, targets, false, Record[], 0, 0)
    end
end

"""
    GFF3.Reader(input::IO;
                index=nothing,
                save_directives::Bool=false,
                skip_features::Bool=false,
                skip_directives::Bool=true,
                skip_comments::Bool=true)

    GFF3.Reader(input::AbstractString;
                index=:auto,
                save_directives::Bool=false,
                skip_features::Bool=false,
                skip_directives::Bool=true,
                skip_comments::Bool=true)

Create a reader for data in GFF3 format.

The first argument specifies the data source. When it is a filepath that ends with *.bgz*, it is considered to be block compression file format (BGZF) and the function will try to find a tabix index file (<filename>.tbi) and read it if any.
See <http://www.htslib.org/doc/tabix.html> for bgzip and tabix tools.

Arguments
---------
- `input`: data source (`IO` object or filepath)
- `index`: path to a tabix file
- `save_directives`: flag to save directive records (which can be accessed with `GFF3.directives`)
- `skip_features`: flag to skip feature records
- `skip_directives`: flag to skip directive records
- `skip_comments`:  flag to skip comment records
"""
function Reader(input::IO; index=nothing, save_directives::Bool=false, skip_features::Bool=false, skip_directives::Bool=true, skip_comments::Bool=true)

    if isa(index, AbstractString)
        index = Indexes.Tabix(index)
    end

    if isa(input, TranscodingStream)
        return Reader(input, index, save_directives, skip_features, skip_directives, skip_comments)
    end

    stream = TranscodingStreams.NoopStream(input)

    return Reader(stream, index, save_directives, skip_features, skip_directives, skip_comments)

end

function Reader(filepath::AbstractString; index=:auto, save_directives::Bool=false, skip_features::Bool=false, skip_directives::Bool=true, skip_comments::Bool=true)
    if isa(index, Symbol) && index != :auto
        throw(ArgumentError("invalid index argument: ':$(index)'"))
    end
    if endswith(filepath, ".bgz")
        input = BGZFStreams.BGZFStream(filepath)
        if index == :auto
            index = Indexes.findtabix(filepath)
        end
    else
        input = open(filepath)
    end
    return Reader(input, index=index, save_directives=save_directives, skip_features=skip_features, skip_directives=skip_directives, skip_comments=skip_comments)
end

function Base.eltype(::Type{<:Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.eof(reader::Reader)
    return reader.state.filled || eof(reader.state.stream)
end

function Base.close(reader::Reader)
    # make trailing directives accessable
    reader.directive_count = reader.preceding_directive_count
    reader.preceding_directive_count = 0
    close(BioGenerics.IO.stream(reader))
end

function Base.read!(reader::Reader, record::Record)
    return readrecord!(reader.state.stream, reader, record)
end

function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    return index!(stream, record)
end

function IntervalCollection(reader::Reader)
    intervals = collect(Interval{Record}, reader)
    return IntervalCollection(intervals, true)
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::Interval)
    if reader.index === nothing
        throw(ArgumentError("index is null"))
    end
    return Indexes.TabixOverlapIterator(reader, interval)
end


"""
Return all directives that preceded the last GFF entry parsed as an array of strings.

Directives at the end of the file can be accessed by calling `close(reader)` and then `directives(reader)`.
"""
function directives(reader::Reader)
    ret = String[]
    for i in lastindex(reader.directives)-reader.directive_count+1:lastindex(reader.directives)
        push!(ret, content(reader.directives[i]))
    end
    return ret
end

"""
Return true if the GFF3 stream is at its end and there is trailing FASTA data.
"""
function hasfasta(reader::Reader)
    if eof(reader)
        return reader.found_fasta
    end
    error("GFF3 file must be read until the end before any FASTA sequences can be accessed")
end

"""
Return a FASTA.Reader initialized to parse trailing FASTA data.

Throws an exception if there is no trailing FASTA, which can be checked using `hasfasta`.
"""
function getfasta(reader::Reader) #TODO: move responsibility to FASTX.jl.
    if !hasfasta(reader)
        error("GFF3 file has no FASTA data")
    end
    return FASTA.Reader(reader.state.stream)
end

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

const record_machine, body_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    rep1 = Automa.RegExp.rep1
    alt = Automa.RegExp.alt
    opt = Automa.RegExp.opt

    feature = let
        seqid = re"[a-zA-Z0-9.:^*$@!+_?\-|%]*"
        seqid.actions[:enter] = [:pos]
        seqid.actions[:exit]  = [:feature_seqid]

        source = re"[ -~]*"
        source.actions[:enter] = [:pos]
        source.actions[:exit]  = [:feature_source]

        type_ = re"[ -~]*"
        type_.actions[:enter] = [:pos]
        type_.actions[:exit]  = [:feature_type_]

        start = re"[0-9]+|\."
        start.actions[:enter] = [:pos]
        start.actions[:exit]  = [:feature_start]

        end_ = re"[0-9]+|\."
        end_.actions[:enter] = [:pos]
        end_.actions[:exit]  = [:feature_end_]

        score = re"[ -~]*[0-9][ -~]*|\."
        score.actions[:enter] = [:pos]
        score.actions[:exit]  = [:feature_score]

        strand = re"[+\-?]|\."
        strand.actions[:enter] = [:feature_strand]

        phase = re"[0-2]|\."
        phase.actions[:enter] = [:feature_phase]

        attributes = let
            char = re"[^=;,\t\r\n]"
            key = rep1(char)
            key.actions[:enter] = [:pos]
            key.actions[:exit]  = [:feature_attribute_key]
            val = rep(char)
            attr = cat(key, '=', val, rep(cat(',', val)))

            cat(rep(cat(attr, ';')), opt(attr))
        end

        cat(seqid,  '\t',
            source, '\t',
            type_,  '\t',
            start,  '\t',
            end_,   '\t',
            score,  '\t',
            strand, '\t',
            phase,  '\t',
            attributes)
    end
    feature.actions[:exit] = [:feature]

    directive = re"##[^\r\n]*"
    directive.actions[:exit] = [:directive]

    comment = re"#([^#\r\n][^\r\n]*)?"
    comment.actions[:exit] = [:comment]

    record = alt(feature, directive, comment)
    record.actions[:enter] = [:mark]
    record.actions[:exit]  = [:record]

    blank = re"[ \t]*"

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    body = rep(cat(alt(record, blank), newline))
    body.actions[:exit] = [:body]

    # look-ahead of the beginning of FASTA
    body′ = cat(body, opt('>'))

    return map(Automa.compile, (record, body′))
end)()

const record_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :feature_seqid   => :(record.seqid  = pos:@relpos(p-1)),
    :feature_source  => :(record.source = pos:@relpos(p-1)),
    :feature_type_   => :(record.type_  = pos:@relpos(p-1)),
    :feature_start   => :(record.start  = pos:@relpos(p-1)),
    :feature_end_    => :(record.end_   = pos:@relpos(p-1)),
    :feature_score   => :(record.score  = pos:@relpos(p-1)),
    :feature_strand  => :(record.strand = @relpos(p)),
    :feature_phase   => :(record.phase  = @relpos(p)),
    :feature_attribute_key => :(push!(record.attribute_keys, pos:@relpos(p-1))),
    :feature         => :(record.kind = :feature),
    :directive       => :(record.kind = :directive),
    :comment         => :(record.kind = :comment),
    :record          => quote
        appendfrom!(record.data, 1, data, @markpos, p-@markpos)
        record.filled = 1:(p-@markpos)
    end
)

context = Automa.CodeGenContext(
    generator = :goto,
    checkbounds = false,
    loopunroll = 0
)

Automa.Stream.generate_reader(
    :index!,
    record_machine,
    arguments = (:(record::Record),),
    actions = record_actions,
    context = context,
    returncode = quote
        if cs == 0
            return record
        end
        throw(ArgumentError(string("failed to index ", eltype(record), " ~>", repr(String(data[p:min(p+7,p_end)])))))
    end
) |> eval


Automa.Stream.generate_reader(
    :readrecord!,
    body_machine,
    arguments = (:(reader::Reader), :(record::Record)),
    actions = merge(record_actions,
        Dict(
            :countline => :(linenum += 1),
            :record => quote
                appendfrom!(record.data, 1, data, @markpos, p-@markpos)
                record.filled = 1:(p-@markpos)
                if isfeature(record)
                    reader.directive_count = reader.preceding_directive_count
                    reader.preceding_directive_count = 0
                elseif isdirective(record)
                    reader.preceding_directive_count += 1
                    if reader.save_directives
                        push!(reader.directives, copy(record))
                    end
                    if is_fasta_directive(record)
                        reader.found_fasta = true
                        reader.state.filled = true
                    end
                end

                if record.kind ∈ reader.targets
                    found_record = true
                    @escape
                end
            end,
            :body => quote
                if data[p] == UInt8('>')
                    reader.found_fasta = true
                    reader.state.filled = true
                    # HACK: any better way?
                    cs = 0
                    @goto exit
                end
            end,
        )
    ),
    context = context,
    initcode = quote
        cs = reader.state.state
        linenum = reader.state.linenum
        found_record = false
    end,
    loopcode = quote
        if found_record
            @goto __return__
        end
    end,
    returncode = quote

        reader.state.state = cs
        # reader.state.filled |= cs == 0 # Note: if set to true, remains true.
        reader.state.linenum = linenum

        if found_record
            return record
        end

        if cs == 0 || eof(stream)
            throw(EOFError())
        end

        if cs < 0
            error(eltype(Reader), " file format error on line ", linenum, " ~>", repr(String(data[p:min(p+7,p_end)])))
        end

        if p > p_eof ≥ 0
            error("incomplete $(typeof(reader)) input on line ", linenum)
        end

    end
) |> eval
