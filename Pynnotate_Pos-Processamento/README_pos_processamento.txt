
GenBank Toolkit – Pós-processamento de Genes
============================================

Este script automatiza o alinhamento, filtragem e concatenação de genes extraídos de genomas mitocondriais (ou similares), agrupados por gene em arquivos FASTA.

Ele trata:
- Genes codificantes (com tradução e retrotradução)
- Genes rRNA (12S, 16S)
- Genes tRNA (qualquer gene começando com tRNA-)

--------------------------------------------

Como rodar
----------

Exemplo de uso:

python Pynnotate_pos_processamento.py \
  --inputs inputs \
  --output resultados \
  --threads 4 \
  --translation_table 2 \
  --trim_flanks 0.25 \
  --min_seq_occupancy 0.2 \
  --sparse_triple_occupancy 0.3 \
  --save_prefiltered \
  --detect_best_frame \
  --remove_stop_codons \
  --gblocks_mode moderate \
  --mafft_args --maxiterate 500 --op 2.0 --ep 0.5

--------------------------------------------

Estrutura esperada da pasta de entrada
--------------------------------------

A pasta --inputs deve conter subpastas com arquivos .fasta por gene, por exemplo:

inputs/
├── tRNA-Leu/
│   └── tRNA-Leu.fasta
├── 16S/
│   └── 16S.fasta
├── COI/
│   └── COI.fasta

--------------------------------------------

Opções disponíveis
------------------

-i, --inputs                 Pasta com subpastas contendo os FASTAs por gene (obrigatório)
-o, --output                 Caminho base da saída (obrigatório)
-T, --threads                Número de threads (default = 2)
-t, --translation_table      Tabela de tradução genética (ex: 2 = vertebrado mitocondrial)
-f, --trim_flanks            Trimming dos flancos baseado na ocupação
-m, --min_seq_occupancy      Ocupação mínima por sequência (default = 0.2)
-s, --save_prefiltered       Salva alinhamentos intermediários
-d, --detect_best_frame      Detecta frame consenso (para genes codificantes)
-g, --gblocks_mode           Modo do Gblocks: strict, moderate, permissive
--sparse_triple_occupancy    Ocupação mínima dos códons (default = 0.3)
--mafft_args                 Parâmetros adicionais para MAFFT

--------------------------------------------

Pipeline por tipo de gene
-------------------------

Codificantes (COI, ND1, etc):
    Traduz → Alinha (AA) → Back-translate → Filtra → Realinha → Trim

rRNA (12S, 16S):
    MAFFT → Gblocks → trim_flanks

tRNA (tRNA-Leu, etc):
    Igual ao pipeline de rRNA

--------------------------------------------

Saídas geradas
--------------

aligned_by_gene/             FASTAs finais por gene
supermatrix.fasta            Concatenado com todos os genes
partitions.txt               Intervalos por gene (uso geral)
partitions_iqtree-raxml.txt  Formato IQ-TREE/RAxML
partitions_nexus.txt         Formato NEXUS
log.txt                      Resumo do processamento
alignment_stats.txt          Estatísticas por gene
short_sequences.txt          Sequências descartadas
long_sequences.txt           Sequências muito longas
codon_warnings.txt           Códons de parada detectados
reading_frames.tsv           Frames de leitura detectados (se ativado)

--------------------------------------------

Requisitos de Software
----------------------

Este script assume que os comandos `mafft` e `gblocks` estejam disponíveis no terminal como comandos nativos do sistema.

Você pode instalá-los facilmente via Conda com os seguintes comandos:

conda install -c bioconda mafft
conda install -c bioconda gblocks

Ou, alternativamente, utilizando o Mamba para instalações mais rápidas:

mamba install -c bioconda mafft gblocks

Verifique se os comandos `mafft` e `Gblocks` funcionam no terminal antes de executar o script.
