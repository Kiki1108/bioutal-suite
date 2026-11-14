#!/usr/bin/env julia
# bin/jfastp.jl
using BioCore

function print_help()
    println("""
    jfastp.jl - Perfilador y Diagnóstico Ultrarrápido de Secuencias FASTQ (Modo Solo Lectura)

    USO:
      julia --project=. -t <hilos> bin/jfastp.jl -i <archivo.fastq> [OPCIONES]

    OPCIONES OBLIGATORIAS:
      -i, --input <archivo>         Archivo FASTQ de entrada a diagnosticar (no será modificado).

    OPCIONES DE REPORTE Y DIAGNÓSTICO:
      -o, --prefix <texto>          Prefijo para los reportes de salida (por defecto: nombre del FASTQ).
      --html <archivo.html>         Ruta personalizada para el reporte visual HTML.
      --json <archivo.json>         Ruta personalizada para el reporte de métricas JSON.
      --adapter-fasta <archivo.fa>  Archivo FASTA opcional con adaptadores específicos a buscar.
      --phred64                     Indica que las calidades usan Phred+64 (por defecto: Phred+33).
      -q, --quiet                   No imprimir la tabla de diagnóstico en la terminal.
      -h, --help                    Muestra este menú de ayuda y termina.

    EJEMPLO:
      julia --project=. -t 8 bin/jfastp.jl -i SRR12334146.fastq -o qc_SRR12334146
    """)
end

function parse_jfastp_args(args::Vector{String})
    isempty(args) && (print_help(); exit(0))
    opts = Dict{String, Any}(
        "input" => "",
        "prefix" => "",
        "html" => "",
        "json" => "",
        "adapter_fasta" => "",
        "phred_offset" => 0x21,
        "quiet" => false
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            print_help(); exit(0)
        elseif arg in ("-i", "--input")
            i += 1; opts["input"] = args[i]
        elseif arg in ("-o", "--prefix")
            i += 1; opts["prefix"] = args[i]
        elseif arg == "--html"
            i += 1; opts["html"] = args[i]
        elseif arg == "--json"
            i += 1; opts["json"] = args[i]
        elseif arg == "--adapter-fasta"
            i += 1; opts["adapter_fasta"] = args[i]
        elseif arg == "--phred64"
            opts["phred_offset"] = 0x40
        elseif arg in ("-q", "--quiet")
            opts["quiet"] = true
        else
            println(stderr, "Error: Opción desconocida '$arg'. Usa --help para ver las opciones.")
            exit(1)
        end
        i += 1
    end

    if isempty(opts["input"])
        println(stderr, "Error: Debes especificar un archivo de entrada con -i <archivo.fastq>.")
        exit(1)
    end
    return opts
end

function main()
    opts = parse_jfastp_args(ARGS)
    in_fq = opts["input"]
    prefix = isempty(opts["prefix"]) ? in_fq : opts["prefix"]
    html_out = isempty(opts["html"]) ? (prefix * "_jfastp.html") : opts["html"]
    json_out = isempty(opts["json"]) ? (prefix * "_jfastp.json") : opts["json"]

    custom_adp = load_adapters_fasta(opts["adapter_fasta"])

    t0 = time()
    prof = run_profiler_parallel(in_fq; phred_offset=opts["phred_offset"], custom_adapters=custom_adp)
    elapsed = round(time() - t0, digits=2)

    !opts["quiet"] && print_profile_summary(stdout, prof, "DIAGNÓSTICO jfastp: $(basename(in_fq)) ($elapsed s)")
    write_profile_json(json_out, prof)
    write_profile_html(html_out, prof, basename(in_fq))
    !opts["quiet"] && println("Reportes guardados en: $html_out y $json_out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end