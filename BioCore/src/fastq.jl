# src/fastq.jl
using Mmap

"""
    FastqRecordView

Vista sin asignaciones en el heap (zero-allocation) sobre un registro FASTQ en memoria mapeada.
"""
struct FastqRecordView{V<:SubArray{UInt8, 1, Vector{UInt8}}}
    id::V
    seq::V
    qual::V
end

struct FastqMmapReader
    data::Vector{UInt8}
end

function FastqMmapReader(filepath::AbstractString)
    io = open(filepath, "r")
    data = Mmap.mmap(io)
    close(io) # El mapeo permanece activo aunque se cierre el descriptor
    return FastqMmapReader(data)
end

# Función auxiliar inline para encontrar el próximo salto de línea (\n = 0x0a)
@inline function find_newline(data::Vector{UInt8}, start_pos::Int, len::Int)
    pos = start_pos
    @inbounds while pos <= len && data[pos] != 0x0a
        pos += 1
    end
    return pos
end

# Ignorar retorno de carro (\r = 0x0d) en caso de archivos generados en Windows
@inline function trim_cr(data::Vector{UInt8}, start_pos::Int, newline_pos::Int)
    end_pos = newline_pos - 1
    @inbounds if end_pos >= start_pos && data[end_pos] == 0x0d
        end_pos -= 1
    end
    return end_pos
end

# Implementación del protocolo de iteración nativo de Julia
function Base.iterate(reader::FastqMmapReader, pos::Int = 1)
    data = reader.data
    total_len = length(data)
    
    # Saltar líneas vacías al final del archivo
    @inbounds while pos <= total_len && (data[pos] == 0x0a || data[pos] == 0x0d)
        pos += 1
    end
    pos > total_len && return nothing
    
    # Línea 1: @ID
    nl1 = find_newline(data, pos, total_len)
    id_end = trim_cr(data, pos, nl1)
    id_view = @view data[(pos + 1):id_end] # Omitimos el carácter '@'
    
    # Línea 2: Secuencia
    pos2 = nl1 + 1
    nl2 = find_newline(data, pos2, total_len)
    seq_end = trim_cr(data, pos2, nl2)
    seq_view = @view data[pos2:seq_end]
    
    # Línea 3: + (Separador)
    pos3 = nl2 + 1
    nl3 = find_newline(data, pos3, total_len)
    
    # Línea 4: Calidad Phred
    pos4 = nl3 + 1
    nl4 = find_newline(data, pos4, total_len)
    qual_end = trim_cr(data, pos4, nl4)
    qual_view = @view data[pos4:qual_end]
    
    return FastqRecordView(id_view, seq_view, qual_view), nl4 + 1
end

Base.IteratorSize(::Type{FastqMmapReader}) = Base.SizeUnknown()
Base.eltype(::Type{FastqMmapReader}) = FastqRecordView

"""
    mean_phred(rec::FastqRecordView, offset::UInt8 = 0x21) -> Float64

Calcula la calidad Phred promedio al vuelo (por defecto Phred+33, 0x21) usando vectorización SIMD.
"""
@inline function mean_phred(rec::FastqRecordView, offset::UInt8 = 0x21)
    len = length(rec.qual)
    len == 0 && return 0.0
    total_q = 0
    @inbounds @simd for i in 1:len
        total_q += Int(rec.qual[i] - offset)
    end
    return total_q / len
end

"""
    count_ns(rec::FastqRecordView) -> Int

Cuenta cuántas bases ambiguas ('N' o 'n') tiene la lectura sin convertir a String.
"""
@inline function count_ns(rec::FastqRecordView)
    n_count = 0
    @inbounds @simd for i in 1:length(rec.seq)
        b = rec.seq[i]
        n_count += (b == UInt8('N') || b == UInt8('n'))
    end
    return n_count
end