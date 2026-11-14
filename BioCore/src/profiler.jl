# src/profiler.jl

const MAX_CYCLE_LEN = 512

mutable struct SeqProfile
    total_reads::Int
    total_bases::Int
    reads_with_n::Int
    q20_bases::Int
    q30_bases::Int
    gc_bases::Int
    poly_g_reads::Int
    adapter_hits::Int
    min_len::Int
    max_len::Int

    # Distribuciones por posición del ciclo (1 a MAX_CYCLE_LEN)
    pos_total::Vector{Int}
    pos_qual_sum::Vector{Int}
    pos_a::Vector{Int}
    pos_c::Vector{Int}
    pos_g::Vector{Int}
    pos_t::Vector{Int}
    pos_n::Vector{Int}

    # Distribución de longitudes de lectura
    length_hist::Vector{Int}

    function SeqProfile()
        new(0, 0, 0, 0, 0, 0, 0, 0, typemax(Int), 0,
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN),
            zeros(Int, MAX_CYCLE_LEN))
    end
end

function merge_profile!(dest::SeqProfile, src::SeqProfile)
    dest.total_reads  += src.total_reads
    dest.total_bases  += src.total_bases
    dest.reads_with_n += src.reads_with_n
    dest.q20_bases    += src.q20_bases
    dest.q30_bases    += src.q30_bases
    dest.gc_bases     += src.gc_bases
    dest.poly_g_reads += src.poly_g_reads
    dest.adapter_hits += src.adapter_hits
    dest.min_len       = min(dest.min_len, src.min_len)
    dest.max_len       = max(dest.max_len, src.max_len)

    dest.pos_total    .+= src.pos_total
    dest.pos_qual_sum .+= src.pos_qual_sum
    dest.pos_a        .+= src.pos_a
    dest.pos_c        .+= src.pos_c
    dest.pos_g        .+= src.pos_g
    dest.pos_t        .+= src.pos_t
    dest.pos_n        .+= src.pos_n
    dest.length_hist  .+= src.length_hist
    return dest
end

# Semillas de 12 nt de los adaptadores universales más comunes para detección automática
const DEFAULT_ADAPTER_SEEDS = [
    Vector{UInt8}(codeunits("AGATCGGAAGAG")), # Illumina Universal / TruSeq
    Vector{UInt8}(codeunits("TGGAATTCTCGG")), # Illumina Small RNA 3'
    Vector{UInt8}(codeunits("CTGTCTCTTATA")), # Nextera Transposase
    Vector{UInt8}(codeunits("GATCGTCGGACT"))  # Illumina Small RNA alternativo
]

@inline function has_adapter_seed(seq::AbstractVector{UInt8}, seeds::Vector{Vector{UInt8}})
    len = length(seq)
    @inbounds for seed in seeds
        slen = length(seed)
        len < slen && continue
        for start_pos in 1:(len - slen + 1)
            match = true
            for j in 1:slen
                if (seq[start_pos + j - 1] & 0xDF) != seed[j]
                    match = false
                    break
                end
            end
            match && return true
        end
    end
    return false
end

@inline function inspect_read!(prof::SeqProfile, rec::FastqRecordView,
                               offset::UInt8, seeds::Vector{Vector{UInt8}})
    seq = rec.seq
    qual = rec.qual
    len = length(seq)
    len == 0 && return

    prof.total_reads += 1
    prof.total_bases += len
    len < prof.min_len && (prof.min_len = len)
    len > prof.max_len && (prof.max_len = len)

    if len <= MAX_CYCLE_LEN
        @inbounds prof.length_hist[len] += 1
    end

    q20_ascii = offset + 0x14
    q30_ascii = offset + 0x1e
    has_n = false
    max_c = min(len, MAX_CYCLE_LEN)

    @inbounds for i in 1:max_c
        b = seq[i] & 0xDF
        q_raw = qual[i]
        q_val = Int(q_raw - offset)

        prof.pos_total[i] += 1
        prof.pos_qual_sum[i] += q_val
        prof.q20_bases += (q_raw >= q20_ascii)
        prof.q30_bases += (q_raw >= q30_ascii)

        if b == UInt8('A')
            prof.pos_a[i] += 1
        elseif b == UInt8('C')
            prof.pos_c[i] += 1
            prof.gc_bases += 1
        elseif b == UInt8('G')
            prof.pos_g[i] += 1
            prof.gc_bases += 1
        elseif b == UInt8('T') || b == UInt8('U')
            prof.pos_t[i] += 1
        else
            prof.pos_n[i] += 1
            has_n = true
        end
    end

    has_n && (prof.reads_with_n += 1)

    # Detectar cola PolyG (>= 10 G en el extremo 3')
    if len >= 10 && trim_poly_g_tail(seq, 1, len, 10) < len
        prof.poly_g_reads += 1
    end

    # Detectar presencia de adaptadores
    if has_adapter_seed(seq, seeds)
        prof.adapter_hits += 1
    end
end

"""
    run_profiler_parallel(input_fastq::String; phred_offset::UInt8=0x21, custom_adapters=Vector{UInt8}[]) -> SeqProfile

Analiza un archivo FASTQ en paralelo usando Mmap sin modificar ni escribir lecturas.
"""
function run_profiler_parallel(input_fastq::String;
                               phred_offset::UInt8 = 0x21,
                               custom_adapters::Vector{Vector{UInt8}} = Vector{UInt8}[])
    reader = FastqMmapReader(input_fastq)
    data = reader.data
    total_len = length(data)
    n_threads = Threads.nthreads()

    seeds = isempty(custom_adapters) ? DEFAULT_ADAPTER_SEEDS : custom_adapters

    boundaries = zeros(Int, n_threads + 1)
    boundaries[1] = sync_fastq_record(data, 1)
    for t in 2:n_threads
        approx_pos = cld(total_len * (t - 1), n_threads)
        boundaries[t] = sync_fastq_record(data, approx_pos)
    end
    boundaries[n_threads + 1] = total_len + 1

    thread_profiles = [SeqProfile() for _ in 1:n_threads]

    Threads.@threads for t in 1:n_threads
        pos = boundaries[t]
        limit = boundaries[t + 1]
        local_prof = thread_profiles[t]

        while pos < limit
            rec_and_next = iterate(reader, pos)
            rec_and_next === nothing && break
            rec, next_pos = rec_and_next
            pos = next_pos
            inspect_read!(local_prof, rec, phred_offset, seeds)
        end
    end

    global_prof = SeqProfile()
    for p in thread_profiles
        merge_profile!(global_prof, p)
    end
    global_prof.min_len == typemax(Int) && (global_prof.min_len = 0)
    return global_prof
end