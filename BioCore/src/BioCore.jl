# src/BioCore.jl
module BioCore

include("bitseq.jl")
include("seqbuffer.jl")
include("fastq.jl")
include("trimmer.jl")
include("profiler.jl")
include("report.jl")

export ShortSeq2Bit, LongSeq2Bit,
       encode_short, encode_long, reverse_complement, canonical, hamming_distance, extract_kmer,
       FastqMmapReader, FastqRecordView, sync_fastq_record,
       FastpConfig, FastpStats, load_adapters_fasta, run_fastp_parallel,
       SeqProfile, run_profiler_parallel,
       print_profile_summary, write_profile_json, write_profile_html

end 