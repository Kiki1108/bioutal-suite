# src/report.jl
using Printf

function print_profile_summary(io::IO, prof::SeqProfile, title::String="DIAGNÓSTICO DE SECUENCIAS")
    reads = max(1, prof.total_reads)
    bases = max(1, prof.total_bases)
    mean_len = prof.total_bases / reads

    println(io, "================================================================")
    println(io, "  \$title")
    println(io, "================================================================")
    @printf(io, "  Lecturas totales        : %d\n", prof.total_reads)
    @printf(io, "  Bases totales           : %d\n", prof.total_bases)
    @printf(io, "  Rango de longitudes     : %d - %d nt (Media: %.1f nt)\n", prof.min_len, prof.max_len, mean_len)
    @printf(io, "  Bases >= Q20            : %.2f%%\n", prof.q20_bases * 100.0 / bases)
    @printf(io, "  Bases >= Q30            : %.2f%%\n", prof.q30_bases * 100.0 / bases)
    @printf(io, "  Contenido GC            : %.2f%%\n", prof.gc_bases * 100.0 / bases)
    @printf(io, "  Lecturas con bases 'N'  : %d (%.2f%%)\n", prof.reads_with_n, prof.reads_with_n * 100.0 / reads)
    @printf(io, "  Lecturas con cola PolyG : %d (%.2f%%)\n", prof.poly_g_reads, prof.poly_g_reads * 100.0 / reads)
    @printf(io, "  Lecturas con adaptador  : %d (%.2f%%)\n", prof.adapter_hits, prof.adapter_hits * 100.0 / reads)
    println(io, "----------------------------------------------------------------")

    # Sugerencias automáticas basadas en los datos observados
    println(io, "  SUGERENCIAS PARA TOMAR DECISIONES:")
    if prof.adapter_hits * 100.0 / reads > 1.0
        println(io, "  [!] Se detectaron adaptadores (>1%). Usa jtrimmer.jl con -a o --adapter-fasta.")
    else
        println(io, "  [OK] Nivel de adaptadores bajo o ya recortado.")
    end
    if prof.max_len > 0 && prof.pos_total[prof.max_len] > 0
        tail_q = prof.pos_qual_sum[prof.max_len] / prof.pos_total[prof.max_len]
        if tail_q < 28.0
            @printf(io, "  [!] La calidad media al final del ciclo (%d nt) cae a Phred %.1f. Usa --cut-tail -M 28.\n", prof.max_len, tail_q)
        else
            @printf(io, "  [OK] Calidad en el extremo 3' estable (Phred %.1f).\n", tail_q)
        end
    end
    println(io, "================================================================")
end

function write_profile_json(filepath::String, prof::SeqProfile)
    bases = max(1, prof.total_bases)
    reads = max(1, prof.total_reads)
    max_c = min(prof.max_len, MAX_CYCLE_LEN)
    mean_q_per_pos = [prof.pos_total[i] > 0 ? round(prof.pos_qual_sum[i] / prof.pos_total[i], digits=2) : 0.0 for i in 1:max_c]

    open(filepath, "w") do io
        @printf(io, "{\n")
        @printf(io, "  \"total_reads\": %d,\n  \"total_bases\": %d,\n", prof.total_reads, prof.total_bases)
        @printf(io, "  \"min_len\": %d,\n  \"max_len\": %d,\n", prof.min_len, prof.max_len)
        @printf(io, "  \"q20_rate\": %.4f,\n  \"q30_rate\": %.4f,\n  \"gc_rate\": %.4f,\n",
                prof.q20_bases/bases, prof.q30_bases/bases, prof.gc_bases/bases)
        @printf(io, "  \"reads_with_n\": %d,\n  \"poly_g_reads\": %d,\n  \"adapter_hits\": %d,\n",
                prof.reads_with_n, prof.poly_g_reads, prof.adapter_hits)
        @printf(io, "  \"cycle_mean_phred\": [%s]\n}\n", join(mean_q_per_pos, ", "))
    end
end

function write_profile_html(filepath::String, prof::SeqProfile, sample_name::String="Muestra")
    bases = max(1, prof.total_bases)
    reads = max(1, prof.total_reads)
    max_c = max(1, min(prof.max_len, MAX_CYCLE_LEN))

    # 1. Construir curva SVG de calidad media por posición (Ciclo 1..max_c)
    qual_points = IOBuffer()
    for i in 1:max_c
        mq = prof.pos_total[i] > 0 ? (prof.pos_qual_sum[i] / prof.pos_total[i]) : 0.0
        x = 50 + round(Int, (i - 1) * (620 / max(1, max_c - 1)))
        y = 240 - round(Int, clamp(mq, 0.0, 42.0) * (200 / 42.0))
        @printf(qual_points, "%d,%d ", x, y)
    end
    y_q28 = 240 - round(Int, 28.0 * (200 / 42.0))

    # 2. Construir barras SVG de distribución de longitudes
    min_l = max(1, prof.min_len)
    max_l = max(min_l, min(prof.max_len, MAX_CYCLE_LEN))
    max_cnt = maximum(prof.length_hist[min_l:max_l]; init=1)
    svg_bars = IOBuffer()
    num_bins = max_l - min_l + 1
    bar_w = max(8, div(620, max(1, num_bins)))
    for (idx, len) in enumerate(min_l:max_l)
        cnt = prof.length_hist[len]
        bar_h = round(Int, (cnt / max_cnt) * 200)
        x = 50 + (idx - 1) * (bar_w + 3)
        y = 240 - bar_h
        @printf(svg_bars, "<rect x='%d' y='%d' width='%d' height='%d' fill='#0284c7' rx='2'><title>%d nt: %d lecturas</title></rect>\n", x, y, bar_w, bar_h, len, cnt)
        if num_bins <= 35 || idx % 5 == 1
            @printf(svg_bars, "<text x='%d' y='258' font-size='10' text-anchor='middle' fill='#334155'>%d</text>\n", x + div(bar_w, 2), len)
        end
    end
    svg_w = max(720, 70 + num_bins * (bar_w + 3))

    open(filepath, "w") do io
        println(io, """
        <!DOCTYPE html>
        <html lang="es">
        <head>
        <meta charset="UTF-8">
        <title>Reporte jfastp - \$sample_name</title>
        <style>
            body { font-family: system-ui, sans-serif; background: #f8fafc; color: #1e293b; max-width: 960px; margin: 2rem auto; padding: 0 1rem; }
            .card { background: white; border-radius: 8px; padding: 1.5rem; margin-bottom: 1.5rem; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
            table { width: 100%; border-collapse: collapse; }
            th, td { padding: 0.6rem; text-align: left; border-bottom: 1px solid #e2e8f0; }
            th { background: #f1f5f9; }
        </style>
        </head>
        <body>
            <h1>Diagnóstico de Secuencias jfastp (\$sample_name)</h1>
            <div class="card">
                <h2>Resumen General de la Muestra</h2>
                <table>
                    <tr><th>Métrica</th><th>Valor</th></tr>
                    <tr><td>Lecturas Totales</td><td>\$(prof.total_reads)</td></tr>
                    <tr><td>Bases Totales</td><td>\$(prof.total_bases)</td></tr>
                    <tr><td>Rango de Longitud</td><td>$(prof.min_len) a$(prof.max_len) nt</td></tr>
                    <tr><td>Bases Q20 / Q30</td><td>$(round(prof.q20_bases*100/bases, digits=2))\% /$(round(prof.q30_bases*100/bases, digits=2))%</td></tr>
                    <tr><td>Contenido GC</td><td>\$(round(prof.gc_bases*100/bases, digits=2))%</td></tr>
                    <tr><td>Lecturas con 'N'</td><td>$(prof.reads_with_n) ($(round(prof.reads_with_n*100/reads, digits=2))%)</td></tr>
                    <tr><td>Colas PolyG detectadas</td><td>$(prof.poly_g_reads) ($(round(prof.poly_g_reads*100/reads, digits=2))%)</td></tr>
                    <tr><td>Lecturas con semillas de adaptador</td><td>$(prof.adapter_hits) ($(round(prof.adapter_hits*100/reads, digits=2))%)</td></tr>
                </table>
            </div>
            <div class="card">
                <h2>Calidad Phred Promedio por Posición (Ciclo 1 a \$max_c)</h2>
                <svg width="720" height="270">
                    <line x1="50" y1="$y_q28" x2="670" y2="$y_q28" stroke="#16a34a" stroke-dasharray="4" stroke-width="1.5"/>
                    <text x="675" y="\$(y_q28+4)" font-size="11" fill="#16a34a">Q28</text>
                    <line x1="50" y1="240" x2="670" y2="240" stroke="#94a3b8" stroke-width="1.5"/>
                    <polyline fill="none" stroke="#dc2626" stroke-width="2.5" points="\$(String(take!(qual_points)))"/>
                </svg>
            </div>
            <div class="card" style="overflow-x: auto;">
                <h2>Distribución de Longitudes de Lectura (nt)</h2>
                <svg width="\$svg_w" height="280">
                    <line x1="45" y1="240" x2="\$(svg_w-20)" y2="240" stroke="#94a3b8" stroke-width="1.5"/>
                    \$(String(take!(svg_bars)))
                </svg>
            </div>
        </body>
        </html>
        """)
    end
end