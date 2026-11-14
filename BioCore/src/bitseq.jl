# src/bitseq.jl

# Códigos en 2 bits: A=00 (0), C=01 (1), G=10 (2), T/U=11 (3), N/Inválido=0xFF
const BASE_A = 0x00
const BASE_C = 0x01
const BASE_G = 0x02
const BASE_T = 0x03
const BASE_INVALID = 0xFF

# Tabla de búsqueda (LUT) en memoria constante de 256 bytes para conversión sin condicionales (branchless)
const ASCII_TO_2BIT = let tbl = fill(BASE_INVALID, 256)
    tbl[Int('A') + 1] = BASE_A; tbl[Int('a') + 1] = BASE_A
    tbl[Int('C') + 1] = BASE_C; tbl[Int('c') + 1] = BASE_C
    tbl[Int('G') + 1] = BASE_G; tbl[Int('g') + 1] = BASE_G
    tbl[Int('T') + 1] = BASE_T; tbl[Int('t') + 1] = BASE_T
    tbl[Int('U') + 1] = BASE_T; tbl[Int('u') + 1] = BASE_T
    Tuple(tbl)
end

const TWOBIT_TO_ASCII = ('A', 'C', 'G', 'T')

"""
    ShortSeq2Bit

Representa una secuencia corta o k-mer (longitud <= 32 nt) empaquetada en un solo registro UInt64.
Almacena además una máscara `n_mask` de 32 bits para recordar posiciones con bases ambiguas 'N'.
"""
struct ShortSeq2Bit
    data::UInt64    # 2 bits por base (empacados en los bits menos significativos)
    n_mask::UInt32  # Bit 1 si la posición tenía una 'N'
    len::UInt8      # Longitud de la secuencia (0 a 32)
end

"""
    encode_short(bytes::AbstractVector{UInt8}) -> ShortSeq2Bit

Empaqueta un vector de bytes ASCII (hasta 32 nt) directamente en registros de CPU sin alocar memoria.
"""
@inline function encode_short(bytes::AbstractVector{UInt8})
    len = length(bytes)
    @assert len <= 32 "ShortSeq2Bit soporta hasta 32 nt. Usa LongSeq2Bit para secuencias mayores."
    
    data = UInt64(0)
    n_mask = UInt32(0)
    
    @inbounds for i in 1:len
        code = ASCII_TO_2BIT[bytes[i] + 0x01]
        data <<= 2
        n_mask <<= 1
        if code == BASE_INVALID
            # Si es 'N', guardamos 'A' (00) en data y encendemos el bit en n_mask
            n_mask |= UInt32(1)
        else
            data |= UInt64(code)
        end
    end
    return ShortSeq2Bit(data, n_mask, UInt8(len))
end

# Constructor conveniente desde String (para pruebas)
encode_short(s::AbstractString) = encode_short(codeunits(s))

"""
    reverse_complement(seq::ShortSeq2Bit) -> ShortSeq2Bit

Calcula el reverso complementario en ~4 instrucciones de reloj usando bswap y máscaras de bits.
"""
@inline function reverse_complement(seq::ShortSeq2Bit)
    seq.len == 0 && return seq
    
    # 1. Complementar todas las bases con NOT bit a bit (00 <-> 11, 01 <-> 10)
    c = ~seq.data
    
    # 2. Invertir el orden de los 8 bytes usando la instrucción nativa de CPU BSWAP
    c = bswap(c)
    
    # 3. Invertir pares de 4 bits (nibbles) y luego pares de 2 bits (crumbs) dentro de cada byte
    c = ((c & 0x0F0F0F0F0F0F0F0F) << 4) | ((c >> 4) & 0x0F0F0F0F0F0F0F0F)
    c = ((c & 0x3333333333333333) << 2) | ((c >> 2) & 0x3333333333333333)
    
    # 4. Alinear a la derecha según la longitud real de la secuencia
    shift = 64 - (Int(seq.len) << 1)
    rc_data = c >> shift
    
    # 5. Invertir la máscara de 'N's usando bitreverse nativo
    rc_nmask = bitreverse(seq.n_mask) >> (32 - Int(seq.len))
    
    return ShortSeq2Bit(rc_data, rc_nmask, seq.len)
end

"""
    canonical(seq::ShortSeq2Bit) -> ShortSeq2Bit

Retorna el k-mer canónico (el menor lexicográficamente entre la hebra directa y su reverso complementario).
"""
@inline function canonical(seq::ShortSeq2Bit)
    rc = reverse_complement(seq)
    return seq.data <= rc.data ? seq : rc
end

"""
    hamming_distance(a::ShortSeq2Bit, b::ShortSeq2Bit) -> Int

Calcula el número exacto de mismatches entre dos secuencias del mismo largo usando XOR y POPCNT.
"""
@inline function hamming_distance(a::ShortSeq2Bit, b::ShortSeq2Bit)
    diff = a.data ⊻ b.data
    # Colapsar cada par de 2 bits en su bit inferior si hubo diferencia en cualquiera de los dos bits
    mismatches = (diff | (diff >> 1)) & 0x5555555555555555
    # Sumar cualquier posición donde haya una 'N'
    n_diff = UInt64(a.n_mask | b.n_mask)
    return count_ones(mismatches) + count_ones(n_diff)
end

# Decodificación a String (solo cuando se necesita imprimir o exportar)
function Base.String(seq::ShortSeq2Bit)
    len = Int(seq.len)
    buf = Vector{UInt8}(undef, len)
    @inbounds for i in 1:len
        shift_2bit = (len - i) << 1
        shift_1bit = len - i
        if ((seq.n_mask >> shift_1bit) & 0x01) == 0x01
            buf[i] = UInt8('N')
        else
            code = (seq.data >> shift_2bit) & 0x03
            buf[i] = UInt8(TWOBIT_TO_ASCII[code + 1])
        end
    end
    return String(buf)
end