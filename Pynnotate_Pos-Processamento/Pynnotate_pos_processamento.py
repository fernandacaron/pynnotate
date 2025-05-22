#!/usr/bin/env python3

import os
import argparse
import shutil
import tempfile
from statistics import stdev
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import subprocess
from collections import defaultdict

def is_rRNA(gene_name):
    return gene_name in ["12S", "16S"]
    
def is_tRNA(gene_name):
    return gene_name.upper().startswith("TRNA")

def clean_sequence(seq):
    """Substitui qualquer caractere diferente de A, T, C, G por '-' (gap)."""
    return ''.join([nt if nt.upper() in 'ATCG' else '-' for nt in seq.upper()])

def trim_flanks(records, threshold):
    print(f"[INFO] Realizando trimming com threshold {threshold}...")
    if not records:
        return []
    seqs = [str(r.seq) for r in records]
    length = len(seqs[0])
    num_seqs = len(seqs)
    start, end = 0, length
    for i in range(length):
        non_gap = sum(1 for s in seqs if s[i] != '-')
        if non_gap / num_seqs >= threshold:
            start = i
            break
    for i in range(length - 1, -1, -1):
        non_gap = sum(1 for s in seqs if s[i] != '-')
        if non_gap / num_seqs >= threshold:
            end = i + 1
            break
    trimmed = [SeqRecord(Seq(r.seq[start:end]), id=r.id, description="") for r in records]
    return [r for r in trimmed if len(r.seq) > 0]

def detect_stop_codons_post_alignment(records, table_id, detect_best_frame=False):
    warnings = []
    consensus_frame = None

    def infer_consensus_frame(records, table_id):
        frame_counts = [0, 0, 0]
        for r in records:
            seq = str(r.seq).replace('-', '').upper()
            if len(seq) < 3:
                continue
            best_frame = None
            min_stops = float('inf')
            for frame in range(3):
                sub_seq = seq[frame:]
                sub_seq = sub_seq[:len(sub_seq) - (len(sub_seq) % 3)]
                try:
                    aa = Seq(sub_seq).translate(table=table_id)
                    stops = aa[:-1].count("*")
                    if stops < min_stops:
                        min_stops = stops
                        best_frame = frame
                except Exception:
                    continue
            if best_frame is not None:
                frame_counts[best_frame] += 1
        return frame_counts.index(max(frame_counts))

    if detect_best_frame:
        consensus_frame = infer_consensus_frame(records, table_id)
        for r in records:
            seq = str(r.seq).replace('-', '').upper()
            sub_seq = seq[consensus_frame:]
            sub_seq = sub_seq[:len(sub_seq) - (len(sub_seq) % 3)]
            try:
                aa = Seq(sub_seq).translate(table=table_id)
                stops = aa[:-1].count("*")
                if stops > 0:
                    warnings.append(f"{r.id} - {stops} STOP codons (frame consenso = {consensus_frame})")
            except Exception as e:
                warnings.append(f"{r.id} - Translation error: {e}")
    else:
        for r in records:
            seq = str(r.seq).replace('-', '').upper()
            trimmed = seq[:len(seq) - (len(seq) % 3)]
            try:
                aa = Seq(trimmed).translate(table=table_id)
                if "*" in aa[:-1]:
                    warnings.append(f"{r.id} - STOP codon detected (default frame 0)")
            except Exception as e:
                warnings.append(f"{r.id} - Translation error: {e}")

    return warnings, consensus_frame

def run_mafft(input_path, output_path, extra_args=None):
    print(f"[INFO] Rodando MAFFT em {input_path.name}...")
    cmd = ["mafft"]
    if extra_args:
        cmd.extend(extra_args)
    else:
        cmd.append("--auto")
    cmd.append(str(input_path))

    # 🔍 Mostrar o comando real no terminal
    print("[DEBUG] MAFFT command:", " ".join(cmd))

    with open(output_path, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL)

def run_gblocks(input_path, output_path, mode="moderate"):
    print(f"[INFO] Rodando Gblocks em {input_path.name} com modo: {mode}")

    if mode == "strict":
        params = ["-b4=10", "-b5=n"]
    elif mode == "permissive":
        params = ["-b4=2", "-b5=a"]
    else:  # moderate (default)
        params = ["-b4=2", "-b5=h"]

    subprocess.run(["Gblocks", str(input_path), "-t=d"] + params,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    gb_file = str(input_path) + "-gb"
    if os.path.exists(gb_file):
        shutil.move(gb_file, output_path)

def trim_to_met(alignment, min_freq=0.3):
    """
    Trima o alinhamento de aminoácidos até a primeira coluna com 'M'
    presente em pelo menos `min_freq` das sequências.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    for col_idx in range(alignment.get_alignment_length()):
        col = alignment[:, col_idx]
        if col.count("M") / len(col) >= min_freq:
            trimmed = [rec.seq[col_idx:] for rec in alignment]
            return MultipleSeqAlignment([
                SeqRecord(Seq(str(seq)), id=rec.id, description="") for rec, seq in zip(alignment, trimmed)
            ])
    return alignment  # se nenhuma coluna tiver M suficiente, retorna sem cortar   

def trim_flanks_by_occupancy(alignment, threshold):
    """
    Remove colunas nos flancos (início e fim) do alinhamento cuja
    ocupação (proporção de sequências não-gap) está abaixo do limiar.

    Funciona tanto para alinhamentos de aminoácidos quanto de nucleotídeos.

    Args:
        alignment (MultipleSeqAlignment): alinhamento de entrada
        threshold (float): proporção mínima de ocupação para manter coluna (ex: 0.2)

    Returns:
        MultipleSeqAlignment: alinhamento com flancos removidos
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    cols = list(zip(*[list(rec.seq) for rec in alignment]))  # transpor
    start, end = 0, len(cols) - 1

    # Trim do início
    while start < len(cols):
        non_gap = cols[start].count('-')
        occ = 1 - (non_gap / len(cols[start]))
        if occ >= threshold:
            break
        start += 1

    # Trim do fim
    while end > start:
        non_gap = cols[end].count('-')
        occ = 1 - (non_gap / len(cols[end]))
        if occ >= threshold:
            break
        end -= 1

    # Aplicar corte a cada sequência
    trimmed_seqs = [rec.seq[start:end+1] for rec in alignment]

    return MultipleSeqAlignment([
        SeqRecord(Seq(str(seq)), id=rec.id, description="") for rec, seq in zip(alignment, trimmed_seqs)
    ])

def remove_sparse_codons_in_nt_alignment(records, min_occupancy=0.3):
    """
    Remove blocos de 3 colunas (códons) com ocupação média abaixo de um limiar
    em um alinhamento de nucleotídeos já back-transladado.

    Garante preservação da estrutura de leitura.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aln_len = len(records[0].seq)
    num_seqs = len(records)

    keep_indices = []

    for i in range(0, aln_len - aln_len % 3, 3):
        occ = []
        for j in range(3):
            col = [rec.seq[i + j] for rec in records]
            occ.append(sum(1 for b in col if b not in ['-', 'N']) / num_seqs)
        mean_occ = sum(occ) / 3
        if mean_occ >= min_occupancy:
            keep_indices.extend([i, i+1, i+2])

    new_records = []
    for rec in records:
        new_seq = ''.join(rec.seq[i] for i in keep_indices)
        new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=""))

    return new_records

def align_proteico_e_backtranslate(records, gene, out_dir, table_id, trim_threshold=None, trim_to_first_met=False, mafft_args=None):

    # Etapa 1: tradução para aminoácidos
    nt_dict = {}
    translated = []
    for rec in records:
        nt_seq = str(rec.seq).replace("-", "")
        # 🔧 Corrigir para múltiplo de 3 com Ns
        remainder = len(nt_seq) % 3
        if remainder != 0:
            nt_seq += "N" * (3 - remainder)
        try:
            aa_seq = Seq(nt_seq).translate(table=table_id)
        except Exception as e:
            print(f"[ERRO] Falha ao traduzir {rec.id} no gene {gene}: {e}")
            continue
        translated.append(SeqRecord(aa_seq, id=rec.id, description=""))
        nt_dict[rec.id] = nt_seq

    # Etapa 2: alinhamento em AA
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_in:
        SeqIO.write(translated, tmp_in.name, "fasta")
        tmp_out = tmp_in.name + "_aligned.faa"
        run_mafft(Path(tmp_in.name), Path(tmp_out), extra_args=mafft_args)
    aa_aln = MultipleSeqAlignment(list(SeqIO.parse(tmp_out, "fasta")))
    os.remove(tmp_in.name)
    os.remove(tmp_out)

    # Etapa 3: trims em AA
    if trim_to_first_met:
        aa_aln = trim_to_met(aa_aln)
    if trim_threshold is not None:
        aa_aln = trim_flanks_by_occupancy(aa_aln, trim_threshold)

    # Etapa 4: back-translation
    nt_aln = []
    for rec in aa_aln:
        aa_seq = str(rec.seq)
        ungapped_nt = nt_dict[rec.id]
        codons = [ungapped_nt[i*3:i*3+3] if i*3+3 <= len(ungapped_nt) else '---' for i in range(len(aa_seq.replace("-", "")))]
        nt_seq = ""
        codon_idx = 0
        for aa in aa_seq:
            if aa == "-":
                nt_seq += "---"
            else:
                nt_seq += codons[codon_idx]
                codon_idx += 1
        nt_aln.append(SeqRecord(Seq(nt_seq), id=rec.id, description=""))

    return nt_aln

def align_gene(gene, files, out_dir, table_id, trim_threshold, column_gap_threshold,
               min_seq_occupancy, save_prefiltered, prefilter_dir, detect_frame,
               gblocks_mode, mafft_args, trim_to_first_met=False,
               remove_sparse_triples=False, sparse_triple_occupancy=0.3, remove_stop_codons=False):

    EMPTY_RESULT = ([], 0, [], 0, 0, 0, None, [], 0)

    print(f"[INFO] Alinhando gene {gene}...")
    all_records = []
    for f in files:
        records = list(SeqIO.parse(f, "fasta"))
        if not records:
            print(f"[WARNING] Arquivo vazio ou inválido: {f}")
            continue
        for r in records:
            r.seq = Seq(clean_sequence(str(r.seq)))
            all_records.append(r)

    if not all_records:
        return EMPTY_RESULT

    lengths = [r.seq.count("A") + r.seq.count("T") + r.seq.count("C") + r.seq.count("G") for r in all_records]
    avg_len = sum(lengths) / len(lengths)
    std_len = stdev(lengths) if len(lengths) > 1 else 0
    max_len = avg_len + 10 * std_len

    filtered_all = []
    long_warnings = []
    for r, l in zip(all_records, lengths):
        if l <= max_len:
            filtered_all.append(r)
        else:
            long_warnings.append(f"{gene} - {r.id}: {l} bp (limite = {int(max_len)} bp)")
    all_records = filtered_all

    if not all_records:
        return EMPTY_RESULT

    stop_warnings = []
    col_removed = 0
    short_filtered = 0
    original_length = 0
    consensus_frame = None

    if is_rRNA(gene) or is_tRNA(gene):
        temp_input = out_dir / f"{gene}_temp_in.fasta"
        temp_output = out_dir / f"{gene}_aligned.fasta"
        SeqIO.write(all_records, temp_input, "fasta")
        run_mafft(temp_input, temp_output, extra_args=mafft_args)

        if temp_output.exists():
            pre_gblocks_records = list(SeqIO.parse(temp_output, "fasta"))
            if pre_gblocks_records:
                original_length = len(pre_gblocks_records[0].seq)

        if save_prefiltered and prefilter_dir:
            inter_path = prefilter_dir / f"{gene}_intermediate.fasta"
            SeqIO.write(SeqIO.parse(temp_output, "fasta"), inter_path, "fasta")
            if inter_path.exists():
                inter_records = list(SeqIO.parse(inter_path, "fasta"))
                if inter_records:
                    original_length = len(inter_records[0].seq)
                    
        if gblocks_mode:
            run_gblocks(temp_output, temp_output, mode=gblocks_mode)

        final_records = list(SeqIO.parse(temp_output, "fasta"))
        Path(temp_input).unlink(missing_ok=True)

        final_records = trim_flanks(final_records, trim_threshold)
        if not final_records or any(len(r.seq) == 0 for r in final_records):
            print(f"[WARNING] {gene}: trimming resultou em sequência(s) vazia(s). Gene descartado.")
            return EMPTY_RESULT

    else:
        final_records = align_proteico_e_backtranslate(
            all_records, gene, out_dir, table_id,
            trim_threshold=trim_threshold,
            trim_to_first_met=trim_to_first_met,
            mafft_args=mafft_args
        )

        if len(final_records) == 0:
            return EMPTY_RESULT
            
        if save_prefiltered and prefilter_dir:
            inter_path = prefilter_dir / f"{gene}_intermediate.fasta"
            SeqIO.write(final_records, inter_path, "fasta")
            if inter_path.exists():
                inter_records = list(SeqIO.parse(inter_path, "fasta"))
                if inter_records:
                    original_length = len(inter_records[0].seq)

        if remove_sparse_triples:
            print(f"[INFO] {gene}: removendo códons esparsos pós back-translation...")
            final_records = remove_sparse_codons_in_nt_alignment(final_records, min_occupancy=sparse_triple_occupancy)

        print(f"[INFO] {gene}: realinhando NT após filtros...")
        temp_realign_in = out_dir / f"{gene}_realign_input.fasta"
        temp_realign_out = out_dir / f"{gene}_realign_output.fasta"
        SeqIO.write(final_records, temp_realign_in, "fasta")
        run_mafft(temp_realign_in, temp_realign_out, extra_args=mafft_args)
        final_records = list(SeqIO.parse(temp_realign_out, "fasta"))
        Path(temp_realign_in).unlink(missing_ok=True)
        Path(temp_realign_out).unlink(missing_ok=True)

        final_records = trim_flanks(final_records, threshold=0.25)
        if not final_records or any(len(r.seq) == 0 for r in final_records):
            print(f"[WARNING] {gene}: trim_flanks final removeu todas as bases. Gene descartado.")
            return EMPTY_RESULT
        print(f"[INFO] {gene}: trim_flanks final aplicado (0.25)")
        if not final_records:
            return EMPTY_RESULT
    
    if save_prefiltered and prefilter_dir:
        pre_out = prefilter_dir / f"{gene}_post_filtered.fasta"
        #SeqIO.write(final_records, pre_out, "fasta")
        SeqIO.write([SeqRecord(Seq(str(r.seq).upper()), id=r.id, description="") for r in final_records],pre_out,"fasta")

    if column_gap_threshold is not None:
        print(f"[INFO] Removendo colunas com > {(1 - column_gap_threshold)*100:.1f}% gaps em {gene}...")
        num_seqs = len(final_records)
        aln_len = len(final_records[0].seq) if final_records else 0
        to_keep = []
        for i in range(aln_len):
            n_gaps = sum(1 for r in final_records if r.seq[i] == '-')
            if (n_gaps / num_seqs) <= (1 - column_gap_threshold):
                to_keep.append(i)
        col_removed = aln_len - len(to_keep)
        final_records = [
            SeqRecord(Seq(''.join(r.seq[i] for i in to_keep)), id=r.id, description="")
            for r in final_records
        ]
        print(f"[INFO] {col_removed} colunas removidas de {gene}")

    aln_len = len(final_records[0].seq) if final_records else 0
    filtered_records = []
    for r in final_records:
        occupancy = sum(1 for b in r.seq if b != '-') / aln_len
        if occupancy >= min_seq_occupancy:
            filtered_records.append(r)
        else:
            short_filtered += 1
    final_records = filtered_records

    if not final_records:
        return EMPTY_RESULT

    if detect_frame and not is_rRNA(gene) and not is_tRNA(gene):
        frame_counts = {1: 0, 2: 0, 3: 0}
        for r in final_records:
            best_frame = None
            fewest_stops = float("inf")
            for frame in [0, 1, 2]:
                seq = str(r.seq[frame:]).replace("-", "")
                seq = seq[:len(seq) - len(seq) % 3]
                try:
                    aa = Seq(seq).translate(table=table_id)
                    n_stops = aa.count("*")
                    if n_stops < fewest_stops:
                        best_frame = frame
                        fewest_stops = n_stops
                except Exception:
                    continue
            if best_frame is not None:
                frame_counts[best_frame + 1] += 1
        consensus_frame = max(frame_counts, key=frame_counts.get)

    if not is_rRNA(gene) and not is_tRNA(gene):
        stop_warnings, _ = detect_stop_codons_post_alignment(final_records, table_id)
        n_internal_stops = sum(1 for w in stop_warnings if "STOP codon" in w)
        
        if remove_stop_codons:
            ids_com_stop = {w.split()[0] for w in stop_warnings if "STOP codon" in w}
            n_stop_removed = len(ids_com_stop)
            final_records = [r for r in final_records if r.id not in ids_com_stop]
        else:
            n_stop_removed = 0
        
    else:
        stop_warnings = []
        n_internal_stops = 0
        n_stop_removed = 0

    if long_warnings:
        with open(out_dir.parent / "long_sequences.txt", "a") as longf:
            for line in long_warnings:
                longf.write(line + "\n")

    final_path = out_dir / f"{gene}_final.fasta"
    final_records_n = [
        SeqRecord(Seq(str(r.seq).upper()), id=r.id, description="")
        for r in final_records
    ]
    SeqIO.write(final_records_n, final_path, "fasta")

    length = len(final_records[0].seq) if final_records else 0

    return (
        final_records,
        length,
        stop_warnings,
        col_removed,
        short_filtered,
        original_length,
        consensus_frame,
        long_warnings,
        n_internal_stops,
        n_stop_removed
    )
                                                                
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", required=True, help="Pasta com subpastas contendo os FASTAs por gene")
    parser.add_argument("-o", "--output", required=True, help="Pasta base de saída")
    parser.add_argument("-T", "--threads", type=int, default=2, help="Número de threads (default = 2)")
    parser.add_argument("-t", "--translation_table", type=int, default=1, help="Tabela de tradução genética")
    parser.add_argument("-f", "--trim_flanks", type=float, default=None, help="Trimming dos flancos")
    parser.add_argument("-c", "--column_gap_threshold", type=float, default=None, help="Threshold para colunas com gaps")
    parser.add_argument("-m", "--min_seq_occupancy", type=float, default=0.2, help="Ocupação mínima por sequência")
    parser.add_argument("-s", "--save_prefiltered", action="store_true", help="Salvar alinhamentos intermediários")
    parser.add_argument("-d", "--detect_best_frame", action="store_true", help="Detectar frame consenso")
    parser.add_argument("-g", "--gblocks_mode", choices=["strict", "moderate", "permissive"], default="moderate", help="Modo do Gblocks")
    parser.add_argument("--sparse_triple_occupancy", type=float, default=0.3, help="Ocupação mínima para manter códons (default = 0.3)")
    parser.add_argument("--remove_stop_codons", action="store_true", help="Remove sequências com códons de parada internos após o alinhamento (somente para CDS)")
    parser.add_argument("--mafft_args", nargs=argparse.REMAINDER, default=[], help="Parâmetros adicionais para MAFFT")
    parser.add_argument("--disable_sparse_triples", action="store_true", help="Desativa remoção de códons esparsos")
    parser.add_argument("--disable_trim_met", action="store_true", help="Desativa trimming até a primeira metionina")

    args = parser.parse_args()
    
    #Ativa por padrão (usuário pode desativar com --disable_trim_met ou --disable_sparse_triples)
    args.trim_to_first_met = not args.disable_trim_met
    args.remove_sparse_triples = not args.disable_sparse_triples

    input_path = Path(args.inputs)
    output_base = Path(args.output)
    aln_dir = output_base / "aligned_by_gene"
    prefilter_dir = output_base / "pre_filtered_alignments" if args.save_prefiltered else None
    aln_dir.mkdir(parents=True, exist_ok=True)
    if prefilter_dir:
        prefilter_dir.mkdir(parents=True, exist_ok=True)

    gene_files = defaultdict(list)
    for folder in input_path.iterdir():
        if folder.is_dir():
            fasta_files = list(folder.glob("*.fasta"))
            if not fasta_files:
                print(f"[WARNING] Nenhum arquivo .fasta encontrado em {folder}")
            for file in fasta_files:
                gene = file.stem
                gene_files[gene].append(file)

    log_lines = []
    concat = []
    sample_ids = set()
    all_stop_warnings = []
    gap_log = []
    short_log = []
    long_log = []
    stats_log = []
    frame_table = {}
    internal_stop_table = {}
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {}
        for gene, files in gene_files.items():
            futures[executor.submit(
                align_gene, gene, files, aln_dir,
                args.translation_table, args.trim_flanks,
                args.column_gap_threshold, args.min_seq_occupancy,
                args.save_prefiltered, prefilter_dir,
                args.detect_best_frame, args.gblocks_mode,
                args.mafft_args,
                args.trim_to_first_met,
                args.remove_sparse_triples,
                args.sparse_triple_occupancy,
                args.remove_stop_codons
            )] = gene

        for fut in futures:
            gene = futures[fut]
            try:
                result = fut.result()
                if not result or not result[0]:
                    print(f"[WARNING] Nenhuma sequência válida processada para {gene}.")
                    continue
            except Exception as e:
                print(f"[ERROR] Falha ao processar {gene}: {e}")
                continue

            records, length, warnings, col_removed, short_filtered, original_length, consensus_frame, gene_long_log, n_internal_stops, n_stop_removed = result
            if length == 0:
                continue

            log_lines.append(f"{gene}: {len(records)} seqs, length = {length}")
            for r in records:
                sample_ids.add(r.id)
            concat.append((gene, records))
            all_stop_warnings += [f"{gene} - {w}" for w in warnings]
            gap_log.append(f"{gene} - {col_removed} columns removed due to gap filtering")
            short_log.append(f"{gene} - {short_filtered} sequences removed due to short length")
            long_log += gene_long_log
            true_removed = original_length - length
            percent_removed = (100 * true_removed / original_length) if original_length > 0 else 0
            stats_log.append(
                f"[{gene}]\n"
                f"  - Original length: {original_length}\n"
                f"  - Columns removed (gap filter): {col_removed}\n"
                f"  - Final length: {length}\n"
                f"  - % columns removed (incl. trimming): {percent_removed:.2f}%\n"
                f"  - Sequences retained: {len(records)}\n"
                f"  - Sequences removed (short): {short_filtered}\n"
                f"  - Sequences removed (stop codon): {n_stop_removed}\n"
            )
            if consensus_frame:
                frame_table[gene] = consensus_frame
            internal_stop_table[gene] = n_internal_stops
            
            gene_order = [
        "tRNA-Leu1", "tRNA-Thr", "tRNA-Pro", "tRNA-Phe", "12S", "tRNA-Val", "16S", "tRNA-Leu2", "ND1", "tRNA-Ile", "tRNA-Gln", "tRNA-Met", "ND2", 
        "tRNA-Trp", "tRNA-Ala", "tRNA-Asn", "tRNA-Cys", "tRNA-Tyr", "COI", "tRNA-Ser1", "tRNA-Asp", "COII", "tRNA-Lys", "ATP8", "ATP6", "COIII", "tRNA-Gly",
        "ND3", "tRNA-Arg", "ND4L", "ND4", "tRNA-His", "tRNA-Ser2", "ND5", "ND6", "tRNA-Glu", "CYTB"
        ]
        
    ordered_concat = sorted(
        concat,
        key=lambda x: gene_order.index(x[0]) if x[0] in gene_order else 999
    )

    partition_log = {}
    current_pos = 1
    for gene, records in ordered_concat:
        length = len(records[0].seq)
        partition_log[gene] = (current_pos, current_pos + length - 1)
        current_pos += length

    all_ids = sorted(sample_ids)
    supermatrix = []
    for sid in all_ids:
        seq = ""
        for gene, records in ordered_concat:
            rec_dict = {r.id: str(r.seq) for r in records}
            seq += rec_dict.get(sid, '-' * (partition_log[gene][1] - partition_log[gene][0] + 1))
        supermatrix.append(SeqRecord(Seq(seq.replace('-', 'N').upper()), id=sid, description=""))

    super_path = output_base / "supermatrix.fasta"
    SeqIO.write(supermatrix, super_path, "fasta")

    # Logs
    with open(output_base / "log.txt", "w") as logf:
        logf.write("\n".join(log_lines))

    with open(output_base / "codon_warnings.txt", "w") as codonlog:
        codonlog.write("\n".join([w for w in all_stop_warnings if "STOP codon" in w or "Translation error" in w]))

    if args.column_gap_threshold is not None and any(gap_log):
        with open(output_base / "gap_filtering.txt", "w") as gaplogfile:
            gaplogfile.write("\n".join(gap_log))

    with open(output_base / "short_sequences.txt", "w") as shortlogfile:
        shortlogfile.write("\n".join(short_log))

    with open(output_base / "long_sequences.txt", "w") as longlogfile:
        longlogfile.write("\n".join(long_log))

    with open(output_base / "alignment_stats.txt", "w") as statfile:
        statfile.write("\n".join(stats_log))

    if frame_table:
        with open(output_base / "reading_frames.tsv", "w") as framelog:
            framelog.write("Gene\tConsensusFrame\tN_Sequences_With_StopCodon\tProporcao\n")
            for gene, frame in frame_table.items():
                n_stops = internal_stop_table.get(gene, 0)
                total = len([r for g, rlist in concat if g == gene for r in rlist])
                framelog.write(f"{gene}\t{frame}\t{n_stops}\t{n_stops}/{total}\n")

    # Arquivos de partições
    ordered_genes = sorted(
        partition_log.items(),
        key=lambda x: gene_order.index(x[0]) if x[0] in gene_order else 999
    )

    with open(output_base / "partitions.txt", "w") as pf:
        for gene, (start, end) in ordered_genes:
            pf.write(f"{gene} = {start}-{end}\n")

    with open(output_base / "partitions_iqtree-raxml.txt", "w") as bf:
        for gene, (start, end) in ordered_genes:
            bf.write(f"DNA, {gene} = {start}-{end}\n")

    with open(output_base / "partitions_nexus.txt", "w") as nf:
        nf.write("begin assumptions;\n")
        for gene, (start, end) in ordered_genes:
            nf.write(f"  charset {gene} = {start}-{end};\n")
        nf.write("end;\n")

    # Resumo
    print("\n[RESUMO FINAL]")
    print(f"Genes processados: {len(partition_log)}")
    for gene, (start, end) in ordered_genes:
        print(f"  - {gene}: {end - start + 1} bp")

    print("\n[INFO] Limpando arquivos temporários...")
    for f in aln_dir.glob("*_temp_in.fasta"):
        f.unlink(missing_ok=True)
    for f in aln_dir.glob("*_aligned.fasta"):
        f.unlink(missing_ok=True)
    for f in aln_dir.glob("*.fasta-gb.htm"):
        f.unlink(missing_ok=True)

    print("\n[INFO] Pós-processamento finalizado com sucesso!")

if __name__ == "__main__":
    main()