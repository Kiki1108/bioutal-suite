# src/trimmer.jl

const MAX_READ_LEN = 512 # Soporta lecturas cortas (sRNA) y estándar de Illumina

Base.@kwdef struct FastpConfig
    cut_front::Bool = false
    cut_tail::Bool = true
    cut_right::Bool = false
    window_size::Int = 4
    cut_mean_qual::Int = 28
    adapters::Vector{Vector{UInt8}} = Vector{UInt8}[]
    min_adapter_overlap::Int = 5
    max_adapter_err::Float64 = 0.15
    trim_poly_g::Bool = true
    poly_min_len::Int = 10
    qualified_phred::Int = 28
    unqualified_pct_limit::Float64 = 20.0
    max_n_count::Int = 2
    min_len::Int = 18
    max_len::Int = 30
    phred_offset::UInt8 = 0x21
end

mutable struct FastpStats
    # Conteos de filtrado
    total_reads::Int
    passed_reads::Int
    failed_quality::Int
    failed_n::Int
    failed_too_short::Int
    failed_too_long::Int
    adapter_trimmed::Int
    
    # Métricas de bases ANTES del filtro
    raw_bases::Int
    raw_q20::Int
    raw_q30::Int
    raw_gc::Int
    
    # Métricas de bases DESPUÉS del filtro
    clean_bases::Int
    clean_q20::Int
    clean_q30::Int
    clean_gc::Int
    
    # Distribución de longitudes de lecturas limpias (clave para ver el pico de 22 pb en miRNAs)
    length_hist::Vector{Int}

    function FastpStats()
        new(0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            zeros(Int, MAX_READ_LEN))
    end
end

"""
    merge_stats!(dest::FastpStats, src::FastpStats)

Suma las estadísticas recolectadas por un hilo secundario al acumulador global.
"""
function merge_stats!(dest::FastpStats, src::FastpStats)
    dest.total_reads += src.total_reads
    dest.passed_reads += src.passed_reads
    dest.failed_quality += src.failed_quality
    dest.failed_n += src.failed_n
    dest.failed_too_short += src.failed_too_short
    dest.failed_too_long += src.failed_too_long
    dest.adapter_trimmed += src.adapter_trimmed
    dest.raw_bases += src.raw_bases
    dest.raw_q20 += src.raw_q20
    dest.raw_q30 += src.raw_q30
    dest.raw_gc += src.raw_gc
    dest.clean_bases += src.clean_bases
    dest.clean_q20 += src.clean_q20
    dest.clean_q30 += src.clean_q30
    dest.clean_gc += src.clean_gc
    dest.length_hist .+= src.length_hist
    return dest
end

# Acumulador rápido de métricas por base (Q20, Q30 y GC)
@inline function accumulate_base_metrics(seq::AbstractVector{UInt8}, qual::AbstractVector{UInt8},
                                         left::Int, right::Int, offset::UInt8)
    q20 = 0; q30 = 0; gc = 0
    q20_ascii = offset + 0x14 # Phred 20
    q30_ascii = offset + 0x1e # Phred 30
    @inbounds @simd for i in left:right
        b = seq[i] & 0xDF # Convierte minúscula a mayúscula limpiando el bit 5
        q = qual[i]
        gc += (b == UInt8('G') || b == UInt8('C'))
        q20 += (q >= q20_ascii)
        q30 += (q >= q30_ascii)
    end
    return q20, q30, gc
end

# (Aquí se mantienen las funciones trim_sliding_window, trim_poly_g_tail y trim_adapters_3p)
@inline function trim_sliding_window(qual::AbstractVector{UInt8}, left::Int, right::Int, cfg::FastpConfig)
    W = cfg.window_size
    target_sum = cfg.cut_mean_qual * W
    offset = cfg.phred_offset

    if cfg.cut_front && (right - left + 1) >= W
        win_sum = 0
        @inbounds for i in left:(left + W - 1); win_sum += Int(qual[i] - offset); end
        @inbounds while left + W <= right && win_sum < target_sum
            win_sum -= Int(qual[left] - offset)
            win_sum += Int(qual[left + W] - offset)
            left += 1
        end
        @inbounds while left <= right && Int(qual[left] - offset) < cfg.cut_mean_qual; left += 1; end
    end

    if cfg.cut_tail && (right - left + 1) >= W
        win_sum = 0
        @inbounds for i in (right - W + 1):right; win_sum += Int(qual[i] - offset); end
        @inbounds while right - W >= left && win_sum < target_sum
            win_sum -= Int(qual[right] - offset)
            win_sum += Int(qual[right - W] - offset)
            right -= 1
        end
        @inbounds while right >= left && Int(qual[right] - offset) < cfg.cut_mean_qual; right -= 1; end
    end

    if cfg.cut_right && (right - left + 1) >= W
        win_sum = 0
        @inbounds for i in left:(left + W - 1); win_sum += Int(qual[i] - offset); end
        cut_pos = right; found_drop = false
        @inbounds for pos in left:(right - W)
            if win_sum < target_sum
                cut_pos = pos - 1; found_drop = true; break
            end
            win_sum -= Int(qual[pos] - offset)
            win_sum += Int(qual[pos + W] - offset)
        end
        found_drop && (right = cut_pos)
    end
    return left, right
end

@inline function trim_poly_g_tail(seq::AbstractVector{UInt8}, left::Int, right::Int, min_len::Int)
    (right - left + 1) < min_len && return right
    mismatches = 0; g_count = 0; cut_pos = right
    @inbounds for i in right:-1:left
        b = seq[i] & 0xDF
        if b == UInt8('G')
            g_count += 1; cut_pos = i - 1
        else
            mismatches += 1
            mismatches > 1 && break
        end
    end
    return g_count >= min_len ? cut_pos : right
end

@inline function trim_adapters_3p(seq::AbstractVector{UInt8}, left::Int, right::Int, cfg::FastpConfig)
    isempty(cfg.adapters) && return right, false
    best_right = right; trimmed = false
    @inbounds for adp in cfg.adapters
        adp_len = length(adp)
        max_start = right - cfg.min_adapter_overlap + 1
        for start_pos in left:max_start
            overlap = min(right - start_pos + 1, adp_len)
            max_mm = floor(Int, overlap * cfg.max_adapter_err)
            mm = 0
            for j in 1:overlap
                if (seq[start_pos + j - 1] & 0xDF) != (adp[j] & 0xDF)
                    mm += 1
                    mm > max_mm && break
                end
            end
            if mm <= max_mm
                best_right = start_pos - 1; trimmed = true; break
            end
        end
    end
    return best_right, trimmed
end

@inline function process_read(rec::FastqRecordView, cfg::FastpConfig, stats::FastpStats)
    stats.total_reads += 1
    seq = rec.seq; qual = rec.qual
    raw_len = length(seq)
    left = 1; right = raw_len

    # Métricas crudas
    rq20, rq30, rgc = accumulate_base_metrics(seq, qual, 1, raw_len, cfg.phred_offset)
    stats.raw_bases += raw_len
    stats.raw_q20 += rq20
    stats.raw_q30 += rq30
    stats.raw_gc += rgc

    # Recortes
    cfg.trim_poly_g && (right = trim_poly_g_tail(seq, left, right, cfg.poly_min_len))
    if !isempty(cfg.adapters) && right >= left
        right, was_trimmed = trim_adapters_3p(seq, left, right, cfg)
        stats.adapter_trimmed += was_trimmed
    end
    if (cfg.cut_front || cfg.cut_tail || cfg.cut_right) && right >= left
        left, right = trim_sliding_window(qual, left, right, cfg)
    end

    # Filtros de longitud
    new_len = right - left + 1
    if new_len < cfg.min_len
        stats.failed_too_short += 1; return false, rec
    end
    if cfg.max_len > 0 && new_len > cfg.max_len
        stats.failed_too_long += 1; return false, rec
    end

    # Filtros de N y calidad porcentual
    unqualified_bases = 0; n_bases = 0
    min_ascii_qual = UInt8(cfg.qualified_phred) + cfg.phred_offset
    @inbounds for i in left:right
        b = seq[i] & 0xDF
        n_bases += (b == UInt8('N'))
        unqualified_bases += (qual[i] < min_ascii_qual)
    end
    if n_bases > cfg.max_n_count
        stats.failed_n += 1; return false, rec
    end
    if (unqualified_bases * 100.0 / new_len) > cfg.unqualified_pct_limit
        stats.failed_quality += 1; return false, rec
    end

    # Métricas post-filtro
    cq20, cq30, cgc = accumulate_base_metrics(seq, qual, left, right, cfg.phred_offset)
    stats.passed_reads += 1
    stats.clean_bases += new_len
    stats.clean_q20 += cq20
    stats.clean_q30 += cq30
    stats.clean_gc += cgc
    if new_len <= MAX_READ_LEN
        @inbounds stats.length_hist[new_len] += 1
    end

    return true, FastqRecordView(rec.id, @view(seq[left:right]), @view(qual[left:right]))
end

"""
    sync_fastq_record(data::Vector{UInt8}, pos::Int) -> Int

Encuentra el byte exacto donde inicia el siguiente registro FASTQ válido a partir de `pos`,
evitando confundir líneas de calidad que comienzan con '@' (Phred 31).
"""
function sync_fastq_record(data::Vector{UInt8}, pos::Int)
    total_len = length(data)
    # Si no estamos en el byte 1, avanzar hasta después del próximo salto de línea
    if pos > 1
        pos = find_newline(data, pos, total_len) + 1
    end

    @inbounds while pos <= total_len
        if data[pos] == UInt8('@')
            nl1 = find_newline(data, pos, total_len)
            pos2 = nl1 + 1
            if pos2 <= total_len
                nl2 = find_newline(data, pos2, total_len)
                seq_len = trim_cr(data, pos2, nl2) - pos2 + 1
                pos3 = nl2 + 1
                if pos3 <= total_len && data[pos3] == UInt8('+')
                    nl3 = find_newline(data, pos3, total_len)
                    pos4 = nl3 + 1
                    if pos4 <= total_len
                        nl4 = find_newline(data, pos4, total_len)
                        qual_len = trim_cr(data, pos4, nl4) - pos4 + 1
                        if seq_len == qual_len && seq_len > 0
                            return pos # Confirmado: es un inicio real de registro FASTQ
                        end
                    end
                end
            end
        end
        pos = find_newline(data, pos, total_len) + 1
    end
    return total_len + 1
end

function load_adapters_fasta(fasta_path::String)
    adapters = Vector{UInt8}[]
    isempty(fasta_path) && return adapters
    for line in eachline(fasta_path)
        s = strip(line)
        if !isempty(s) && !startswith(s, '>')
            push!(adapters, Vector{UInt8}(codeunits(s)))
        end
    end
    return adapters
end

"""
    run_fastp_parallel(input_fastq::String, output_fastq::String, cfg::FastpConfig) -> FastpStats

Divide el archivo Mmap en bloques por hilo, procesa las lecturas en paralelo y escribe en bloques de 4 MB.
"""
function run_fastp_parallel(input_fastq::String, output_fastq::String, cfg::FastpConfig)
    reader = FastqMmapReader(input_fastq)
    data = reader.data
    total_len = length(data)
    n_threads = Threads.nthreads()

    # 1. Calcular los límites exactos de bytes [boundaries[t] : boundaries[t+1]-1] para cada hilo
    boundaries = zeros(Int, n_threads + 1)
    boundaries[1] = sync_fastq_record(data, 1)
    for t in 2:n_threads
        approx_pos = cld(total_len * (t - 1), n_threads)
        boundaries[t] = sync_fastq_record(data, approx_pos)
    end
    boundaries[n_threads + 1] = total_len + 1

    thread_stats = [FastpStats() for _ in 1:n_threads]
    io_lock = ReentrantLock()

    open(output_fastq, "w") do out_io
        Threads.@threads for t in 1:n_threads
            pos = boundaries[t]
            limit = boundaries[t + 1]
            local_stats = thread_stats[t]
            local_buf = IOBuffer(maxsize = 4 * 1024 * 1024) # Buffer RAM de 4 MB por hilo

            while pos < limit
                rec_and_next = iterate(reader, pos)
                rec_and_next === nothing && break
                rec, next_pos = rec_and_next
                pos = next_pos

                passed, clean_rec = process_read(rec, cfg, local_stats)
                if passed
                    write(local_buf, UInt8('@'), clean_rec.id, UInt8('\n'))
                    write(local_buf, clean_rec.seq, UInt8('\n'))
                    write(local_buf, UInt8('+'), UInt8('\n'))
                    write(local_buf, clean_rec.qual, UInt8('\n'))

                    # Volcar a disco solo cuando el buffer del hilo supera 3.5 MB
                    if local_buf.size >= 3_500_000
                        lock(io_lock) do
                            write(out_io, take!(local_buf))
                        end
                    end
                end
            end

            # Volcar el remanente del hilo
            if local_buf.size > 0
                lock(io_lock) do
                    write(out_io, take!(local_buf))
                end
            end
        end
    end

    # Consolidar estadísticas de todos los hilos
    global_stats = FastpStats()
    for st in thread_stats
        merge_stats!(global_stats, st)
    end
    return global_stats
end