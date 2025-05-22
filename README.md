
# 🧬 Pynnotate

**Pynnotate** é uma ferramenta gráfica (GUI) em Python para buscar, baixar e anotar sequências genéticas do GenBank com base em genes, organismos e metadados definidos pelo usuário. Ideal para pesquisas em filogenia, sistemática molecular e bioinformática.

---

## 🚀 Funcionalidades

- Busca por termos personalizados, genes e organismos no GenBank
- Download de sequências genéticas (em lote ou por ID)
- Extração de genes anotados automaticamente
- Agrupamento por sinônimos (via alias_map personalizado)
- Geração de logs de processamento, arquivos FASTA e planilhas Excel
- Interface gráfica intuitiva com suporte a campos de cabeçalho customizáveis

---

## 🖼️ Interface

![screenshot opcional aqui]

---

## 🛠️ Requisitos

- Python 3.8 ou superior
- Bibliotecas:
  - `Biopython`
  - `openpyxl`
  - `tkinter` (nativo)
  - `pillow` (para exibir logos)

Instale via pip:

```bash
pip install biopython openpyxl pillow
```

---

## ▶️ Executando o programa

```bash
python pynnotate.py
```

O programa abrirá uma interface gráfica com todas as opções.

---

## 📦 Compilando para `.exe` (Windows)

Para criar um executável com ícone personalizado:

1. Certifique-se de que o arquivo `pynnotate.ico` está na mesma pasta do script.
2. Use o PyInstaller:

```bash
pyinstaller --onefile --noconsole --icon=pynnotate.ico --name=Pynnotate pynnotate.py
```

O executável final será criado em `dist/Pynnotate.exe`.

---

## 📁 Estrutura dos arquivos gerados

- `GenBank_Output/` → Pasta principal de saída
  - `*_sequencias.fasta` → FASTA com os registros baixados
  - `*_metadados.xlsx` → Planilha com metadados por amostra
  - `*_log.txt` → Log geral do processamento
  - `*_genes/` → Pasta com genes extraídos individualmente
    - `log_genes.txt` → Número de sequências por gene
    - `processamento.log` → Genes ausentes por amostra
    - `aliases_desconhecidos.txt` → Sinônimos não reconhecidos

---

## 🧬 Exemplo de uso

1. Defina um gene (ex: *COI*) e um organismo (ex: *Anura*)
2. Clique em "🔽 Buscar e baixar sequências"
3. O programa irá buscar, baixar e extrair os dados automaticamente
4. Veja os arquivos gerados na pasta `GenBank_Output`

---

## 💻 Autores

(INCLUIR)

**Felipe de Medeiros Magalhães**  
Universidade Federal da Paraíba (UFPB)  


---

## 📄 Licença

Este projeto é de código aberto e livre para fins acadêmicos.  
Contribuições são bem-vindas!
