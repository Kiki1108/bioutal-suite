#!/usr/bin/env julia
# bin/jtrimmer.jl
using BioCore

function print_trimmer_help()
    println("""
    jtrimmer.jl - Recortador y Filtrador Multihilo con Ejecución Automática de jfastp

    USO:
      julia --project=. -t <hilos> bin/jtrimmer.jl -i <crudo.fastq> -o <limpio.fastq> [OPCIONES]

    ENTRADA / SALIDA:
      -i, --input <archivo>           Archivo FASTQ de entrada.
      -o, --output <archivo>          Archivo FASTQ limpio de salida.

    CORTE POR VENTANA DESLIZANTE:
      -5, --cut-front                 Recortar bases de baja calidad en el extremo 5'.
      -3, --cut-tail                  Recortar bases de baja calidad en el extremo 3'.
      -r, --cut-right                 Cortar desde 5' hacia 3' apenas una ventana baje del umbral.
      -W, --cut-window-size <int>     Tamaño de la ventana deslizante (por defecto: 4).
      -M, --cut-mean-quality <int>    Calidad Phred media mínima de la ventana (por defecto: 28).

    RECORTE DE ADAPTADORES Y COLAS POLYG:
      -a, --adapter <secuencia>       Secuencia de adaptador 3' en texto plano.
      --adapter-fasta <archivo.fa>    Archivo FASTA con lista de adaptadores a recortar.
      --min-overlap <int>             Solapamiento mínimo del adaptador (por defecto: 5).
      --err-rate <float>              Tasa máxima de error en el adaptador (por defecto: 0.15).
      -g, --trim-poly-g               Recortar colas artificiales PolyG en 3'.

    FILTROS DE CALIDAD Y LONGITUD:
      -q, --qualified-phred <int>     Calidad Phred mínima para considerar una base buena (por defecto: 28).
      -u, --unqualified-limit <float> % máximo permitido de bases bajo -q (por defecto: 20.0).
      -n, --n-base-limit <int>        Número máximo de bases ambiguas 'N' permitidas (por defecto: 2).
      -l, --length-required <int>     Longitud mínima requerida tras el corte (por defecto: 18).
      --length-limit <int>            Longitud máxima permitida tras el corte, 0=sin límite (por defecto: 0).

    CONTROL DE POST-DIAGNÓSTICO (jfastp AUTOMÁTICO):
      --no-post-qc                    Desactiva la ejecución automática de jfastp sobre el archivo limpio.
      -h, --help                      Muestra este menú de ayuda y termina.

    EJEMPLO (sRNA / miRNA):
      julia --project=. -t 8 bin/jtrimmer.jl -i crudo.fastq -o limpio.fastq \\
            --cut-tail -M 28 -q 28 -u 20 -l 18 --length-limit 30 --adapter-fasta adaptadores.fasta
    """)
end

function parse_trimmer_args(args::Vector{String})
    isempty(args) && (print_trimmer_help(); exit(0))

    in_fq = ""; out_fq = ""
    cut_front = false; cut_tail = false; cut_right = false
    win_size = 4; cut_mean_q = 28
    adapters = Vector{UInt8}[]
    min_ovlp = 5; err_rate = 0.15
    poly_g = false
    qual_phred = 28; unqual_lim = 20.0; max_n = 2
    min_len = 18; max_len = 0
    run_post_qc = true

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            print_trimmer_help(); exit(0)
        elseif arg in ("-i", "--input")
            i += 1; in_fq = args[i]
        elseif arg in ("-o", "--output")
            i += 1; out_fq = args[i]
        elseif arg in ("-5", "--cut-front")
            cut_front = true
        elseif arg in ("-3", "--cut-tail")
            cut_tail = true
        elseif arg in ("-r", "--cut-right")
            cut_right = true
        elseif arg in ("-W", "--cut-window-size")
            i += 1; win_size = parse(Int, args[i])
        elseif arg in ("-M", "--cut-mean-quality")
            i += 1; cut_mean_q = parse(Int, args[i])
        elseif arg in ("-a", "--adapter")
            i += 1; push!(adapters, Vector{UInt8}(codeunits(args[i])))
        elseif arg == "--adapter-fasta"
            i += 1; append!(adapters, load_adapters_fasta(args[i]))
        elseif arg == "--min-overlap"
            i += 1; min_ovlp = parse(Int, args[i])
        elseif arg == "--err-rate"
            i += 1; err_rate = parse(Float64, args[i])
        elseif arg in ("-g", "--trim-poly-g")
            poly_g = true
        elseif arg in ("-q", "--qualified-phred")
            i += 1; qual_phred = parse(Int, args[i])
        elseif arg in ("-u", "--unqualified-limit")
            i += 1; unqual_lim = parse(Float64, args[i])
        elseif arg in ("-n", "--n-base-limit")
            i += 1; max_n = parse(Int, args[i])
        elseif arg in ("-l", "--length-required")
            i += 1; min_len = parse(Int, args[i])
        elseif arg == "--length-limit"
            i += 1; max_len = parse(Int, args[i])
        elseif arg == "--no-post-qc"
            run_post_qc = false
        else
            println(stderr, "Error: Opción desconocida '$arg'. Usa --help para ver el manual.")
            exit(1)
        end
        i += 1
    end

    if isempty(in_fq) || isempty(out_fq)
        println(stderr, "Error: Debes indicar -i <entrada.fastq> y -o <salida.fastq>.")
        exit(1)
    end

    cfg = FastpConfig(
        cut_front = cut_front,
        cut_tail = cut_tail,
        cut_right = cut_right,
        window_size = win_size,
        cut_mean_qual = cut_mean_q,
        adapters = adapters,
        min_adapter_overlap = min_ovlp,
        max_adapter_err = err_rate,
        trim_poly_g = poly_g,
        qualified_phred = qual_phred,
        unqualified_pct_limit = unqual_lim,
        max_n_count = max_n,
        min_len = min_len,
        max_len = max_len
    )
    return in_fq, out_fq, cfg, run_post_qc
end

function main()
    in_fq, out_fq, cfg, run_post_qc = parse_trimmer_args(ARGS)

    println("[1/2] Ejecutando jtrimmer sobre $in_fq con $(Threads.nthreads()) hilos...")
    t0 = time()
    trim_stats = run_fastp_parallel(in_fq, out_fq, cfg)
    t_trim = round(time() - t0, digits=2)

    println("  -> Recorte finalizado en $t_trim s.")
    println("  -> Sobrevivieron: $(trim_stats.passed_reads) / $(trim_stats.total_reads) lecturas.")
    println("  -> Adaptadores cortados: $(trim_stats.adapter_trimmed) | Muy cortas: $(trim_stats.failed_too_short) | Muy largas: $(trim_stats.failed_too_long) | Baja calidad: $(trim_stats.failed_quality)")

    # Ejecución automática de jfastp sobre el archivo limpio recién generado
    if run_post_qc && trim_stats.passed_reads > 0
        println("\n[2/2] Ejecutando jfastp automáticamente sobre el archivo limpio ($out_fq)...")
        t1 = time()
        post_prof = run_profiler_parallel(out_fq; phred_offset=cfg.phred_offset, custom_adapters=cfg.adapters)
        t_qc = round(time() - t1, digits=2)

        print_profile_summary(stdout, post_prof, "POST-QC AUTOMÁTICO (jfastp): $(basename(out_fq)) ($t_qc s)")
        write_profile_json(out_fq * "_jfastp.json", post_prof)
        write_profile_html(out_fq * "_jfastp.html", post_prof, basename(out_fq))
        println("  -> Reportes post-recorte generados: $(out_fq)_jfastp.html y $(out_fq)_jfastp.json")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end