# Compilando o Pynnotate 1.0.2rc1

O PyInstaller precisa ser executado no sistema operacional de destino. Por
isso, compile a versão Windows no Windows, a versão macOS no macOS e a versão
Linux no Linux.

## 1. Preparar o ambiente

Crie um ambiente virtual limpo:

```bash
python -m venv .venv
```

Ative no macOS ou Linux:

```bash
source .venv/bin/activate
```

Ative no Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Se você continuar usando o ambiente antigo do macOS, a sequência também é
válida:

```bash
conda deactivate
source ~/venv_pynnotate/bin/activate
```

Confirme que o Python ativo é o esperado:

```bash
which python
python -c "import sys; print(sys.executable)"
```

No Windows, use `where python` no lugar de `which python`.

Instale as dependências:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt -r requirements-build.txt
```

No Debian/Ubuntu, instale também o Tk antes de compilar:

```bash
sudo apt update
sudo apt install python3-tk
```

Não é necessário instalar `pandas`, `requests` ou `lxml`: eles não são usados
por esta versão do Pynnotate.

## 2. Compilar no macOS

Na raiz da matriz, onde está `Pynnotate_mac.spec`:

```bash
python -m PyInstaller --clean --noconfirm Pynnotate_mac.spec
open dist/Pynnotate.app
```

O aplicativo terá a arquitetura do Python ativo. Use Python ARM para Apple
Silicon e Python x86_64 para Macs Intel.

Se o macOS bloquear somente a cópia que você acabou de compilar localmente:

```bash
xattr -cr dist/Pynnotate.app
```

O comando `chmod` normalmente não é necessário em builds normais do
PyInstaller. Para criar um DMG de teste depois da validação:

```bash
hdiutil create -volname "Pynnotate" \
  -srcfolder dist/Pynnotate.app \
  -ov -format UDZO Pynnotate-mac.dmg
```

Para distribuição pública, o aplicativo deve ser assinado e notarizado.

## 3. Compilar no Windows

Primeiro gere a versão em pasta, mais fácil de testar e diagnosticar:

```powershell
python -m PyInstaller --clean --noconfirm Pynnotate_win.spec
```

Saída: `dist\Pynnotate\Pynnotate.exe`. A pasta inteira deve acompanhar o
executável.

Depois dos testes, gere opcionalmente um único `.exe`, como no procedimento
antigo:

```powershell
python -m PyInstaller --clean --noconfirm Pynnotate_win_onefile.spec
```

Saída: `dist\Pynnotate.exe`. O modo de arquivo único demora um pouco mais para
abrir porque extrai os componentes em uma pasta temporária a cada execução.

Não acrescente `--onefile` ao comando que usa um `.spec`: essa característica
precisa estar definida dentro do próprio arquivo `.spec`.

## 4. Compilar no Linux

```bash
python -m PyInstaller --clean --noconfirm Pynnotate_linux.spec
```

Saída: `dist/Pynnotate/Pynnotate`. Distribua a pasta `dist/Pynnotate` completa.

## 5. Teste mínimo antes da distribuição

1. Abra o editor de sinônimos, salve um pequeno dicionário JSON e feche a janela.
2. Busque dois acessos GenBank manualmente; cancele uma vez e confirme na segunda.
3. Confira que o cancelamento não criou uma pasta de resultados vazia.
4. Repita com **Generate log of missing genes per sample** ligado.
5. Confira `sequences.fasta`, `metadata.xlsx`, `log.txt` e a pasta `genes`.
6. Teste uma busca sem API key e outra com uma NCBI API key.
7. No Windows, execute esses testes na versão em pasta antes de gerar o `.exe` único.
8. No Step 1, use **Check for updates** e confirme que uma falha de internet mantém a base atual.
9. Confira se `synonyms_used.json` foi criado junto dos outros resultados.
