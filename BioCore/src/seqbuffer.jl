# src/seqbuffer.jl

struct LongSeq2Bit
    chunks::Vector{UInt64} # Cada UInt64 guarda 32 bases
    len::Int               # Longitud total en nucleótidos
end

function encode_long(bytes::AbstractVector{UInt8})
    len = length(bytes)
    num_chunks = cld(len, 32)
    chunks = zeros(UInt64, num_chunks)
    
    @inbounds for i in 1:len
        code = ASCII_TO_2BIT[bytes[i] + 0x01]
        val = code == BASE_INVALID ? UInt64(0) : UInt64(code)
        
        chunk_idx = ((i - 1) >> 5) + 1          # División entera por 32
        bit_offset = ((i - 1) & 31) << 1        # Posición dentro del UInt64 (0 a 62)
        chunks[chunk_idx] |= (val << bit_offset)
    end
    return LongSeq2Bit(chunks, len)
end

"""
    extract_kmer(seq::LongSeq2Bit, pos::Int, k::Int) -> UInt64

Extrae un k-mer de longitud `k` (<= 32) comenzando en el índice `pos` (1-based) en tiempo O(1).
"""
@inline function extract_kmer(seq::LongSeq2Bit, pos::Int, k::Int)
    chunk_idx = ((pos - 1) >> 5) + 1
    bit_offset = ((pos - 1) & 31) << 1
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (k << 1)) - UInt64(1)
    
    @inbounds begin
        low = seq.chunks[chunk_idx] >> bit_offset
        if bit_offset + (k << 1) <= 64
            return low & mask
        else
            high = seq.chunks[chunk_idx + 1] << (64 - bit_offset)
            return (low | high) & mask
        end
    end
end