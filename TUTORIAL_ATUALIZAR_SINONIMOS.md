# Tutorial: atualizar a base de sinônimos pelo GitHub

Este tutorial mostra como publicar uma nova versão da base oficial de
sinônimos do Pynnotate usando apenas o navegador e o GitHub Codespaces. Não é
necessário recompilar o Pynnotate nem gerar novos executáveis.

## O que acontece durante uma atualização

Existem dois arquivos relacionados:

```text
pynnotate/data/synonyms_database.json
pynnotate/data/synonyms_manifest.json
```

O primeiro contém os genes e seus sinônimos. O segundo informa ao aplicativo
qual revisão está disponível e contém o SHA-256 usado para verificar se o banco
foi baixado corretamente.

Você edita apenas `synonyms_database.json`. O comando de reconstrução valida o
banco e atualiza `synonyms_manifest.json` automaticamente.

## Antes de editar: auditar aliases desconhecidos

Quando um `log.txt` contém a seção `[UNKNOWN ALIASES]`, use o auditor para
separar nomes realmente novos de nomes que já pertencem a outro dicionário:

```bash
python scripts/audit_unknown_synonyms.py /caminho/para/log.txt
```

O dicionário selecionado na busca é detectado pelo próprio log. Para verificar
o mesmo log contra outro tipo de genoma, informe-o explicitamente:

```bash
python scripts/audit_unknown_synonyms.py /caminho/para/log.txt \
  --genome-type animal_mito \
  --output-dir synonym_audit
```

No Windows PowerShell, o mesmo comando pode ser escrito em uma linha:

```powershell
python scripts/audit_unknown_synonyms.py "C:\caminho\log.txt" --genome-type animal_mito --output-dir synonym_audit
```

O comando cria três arquivos:

- `audit_report.md`: resumo e alerta de possível dicionário incorreto;
- `unknown_aliases_audit.csv`: alias, frequência, exemplos e classificação;
- `synonym_proposals.json`: sugestões agrupadas por gene para revisão.

O auditor nunca modifica a base. Revise as sugestões antes de copiá-las para
`synonyms_database.json`. Casos de pseudogene e nomes ambíguos ficam separados.
Diferenças apenas entre maiúsculas e minúsculas não exigem novos sinônimos,
pois a comparação do Pynnotate não diferencia caixa.

## Antes da primeira atualização

Confirme que estes arquivos estão no repositório:

```text
pynnotate/data/synonyms_database.json
pynnotate/data/synonyms_manifest.json
pynnotate/main/synonym_updates.py
scripts/rebuild_synonym_manifest.py
```

Se algum estiver ausente, copie-o da pasta `Pynnotate-1.0.2rc1` para o mesmo
caminho no GitHub. Isso é necessário somente na primeira configuração.

## Parte 1 — abrir o repositório no Codespaces

1. Acesse `https://github.com/fernandacaron/pynnotate`.
2. Clique no botão verde **Code**.
3. Abra a aba **Codespaces**.
4. Clique em **Create codespace on main**.
5. Aguarde a abertura do editor. Ele se parece com o Visual Studio Code e roda
   inteiramente no navegador.

Documentação oficial: [Creating a codespace for a
repository](https://docs.github.com/en/codespaces/developing-in-a-codespace/creating-a-codespace-for-a-repository).

## Parte 2 — editar os sinônimos

1. No explorador de arquivos à esquerda, abra `pynnotate`.
2. Abra a pasta `data`.
3. Clique em `synonyms_database.json`.
4. Localize um dos quatro dicionários:

```json
"animal_mito": {}
"animal_nuclear": {}
"plant_mito": {}
"plant_chloro": {}
```

O conteúdo real aparece dentro das chaves `{}`.

### Adicionar um sinônimo a um gene existente

Localize o gene. Por exemplo:

```json
"COI": [
  "COI",
  "CO1",
  "COX1"
]
```

Adicione uma vírgula ao item anterior e inclua o novo sinônimo:

```json
"COI": [
  "COI",
  "CO1",
  "COX1",
  "NOVO SINÔNIMO"
]
```

Não coloque vírgula depois do último item da lista.

### Adicionar um novo gene canônico

Um novo gene deve ter um nome principal e uma lista de nomes reconhecidos:

```json
"NOVO_GENE": [
  "NOVO_GENE",
  "PRIMEIRO SINÔNIMO",
  "SEGUNDO SINÔNIMO"
]
```

Separe esse bloco do gene anterior com uma vírgula.

## Parte 3 — aumentar a revisão

No início de `synonyms_database.json`, localize:

```json
"database_revision": 1,
"dictionary_version": "2026.08.23",
"schema_version": 1
```

Em cada atualização publicada:

1. Aumente `database_revision` em uma unidade.
2. Altere `dictionary_version` para a data ou identificação da nova versão.
3. Não altere `schema_version`.

Exemplo da segunda atualização:

```json
"database_revision": 2,
"dictionary_version": "2026.09.15",
"schema_version": 1
```

Use `Ctrl+S` no Windows/Linux ou `Command+S` no macOS para salvar.

## Parte 4 — reconstruir o manifest

1. No menu superior do Codespaces, clique em **Terminal**.
2. Clique em **New Terminal**.
3. Confira que o terminal está na raiz do repositório. A linha deve terminar
   aproximadamente em `/workspaces/pynnotate`.
4. Execute:

```bash
python scripts/rebuild_synonym_manifest.py
```

Se o comando `python` não existir, tente:

```bash
python3 scripts/rebuild_synonym_manifest.py
```

Resultado esperado:

```text
Validated revision 2
Updated .../pynnotate/data/synonyms_manifest.json
```

Esse comando:

- verifica a estrutura do JSON;
- detecta aliases atribuídos a genes diferentes no mesmo dicionário;
- copia revisão e versão para o manifest;
- calcula o novo SHA-256;
- atualiza `synonyms_manifest.json`.

Sempre execute esse comando depois da última alteração no banco. Se editar o
banco novamente, execute o comando novamente.

## Parte 5 — conferir as alterações

No terminal, execute:

```bash
git status
```

Os dois arquivos esperados são:

```text
modified: pynnotate/data/synonyms_database.json
modified: pynnotate/data/synonyms_manifest.json
```

Para visualizar exatamente o que mudou:

```bash
git diff -- pynnotate/data/synonyms_database.json pynnotate/data/synonyms_manifest.json
```

Confirme que:

- os novos sinônimos estão corretos;
- a revisão aumentou;
- a versão/data está correta;
- o manifest ganhou um novo `sha256`.

## Parte 6 — publicar no GitHub

No mesmo terminal, execute:

```bash
git add pynnotate/data/synonyms_database.json pynnotate/data/synonyms_manifest.json
```

Crie o commit, ajustando o número da revisão na mensagem:

```bash
git commit -m "Update synonym database to revision 2"
```

Publique no branch principal:

```bash
git push origin main
```

Depois do `git push`, as alterações ficam disponíveis no repositório oficial.
O Codespace pode ser fechado; as alterações publicadas não serão perdidas.

## Parte 7 — testar no Pynnotate

1. Abra o Pynnotate.
2. No Step 1, clique em **Check for updates**.
3. O aplicativo deve mostrar a revisão nova e um resumo das alterações.
4. Confirme a instalação.
5. Verifique se o status mostra a nova `dictionary_version`.
6. Opcionalmente, selecione o gene atualizado e faça uma busca pequena.

Não é necessário recompilar o executável. Todos os usuários que clicarem em
**Check for updates** receberão a nova base oficial.

## Erros comuns

### `The synonym database is already up to date`

Possíveis causas:

- `database_revision` não foi aumentada;
- o manifest não foi reconstruído;
- os arquivos ainda não foram publicados com `git push`.

### `The synonym database checksum does not match the manifest`

O banco foi alterado depois que o manifest foi reconstruído. Execute novamente:

```bash
python scripts/rebuild_synonym_manifest.py
```

Depois faça um novo commit dos dois JSONs.

### `Invalid JSON`

Geralmente significa uma vírgula, aspas ou colchete ausente. Verifique a região
editada. O editor costuma sublinhar o ponto do erro em vermelho.

### `Conflicting alias`

O mesmo nome foi colocado em dois genes diferentes dentro do mesmo dicionário.
Remova o alias de um deles ou corrija o nome.

### `No such file or directory: scripts/rebuild_synonym_manifest.py`

O terminal não está na raiz do repositório ou o script ainda não foi enviado
ao GitHub. Volte à raiz com:

```bash
cd /workspaces/pynnotate
```

Se o script estiver ausente, copie-o do pacote `Pynnotate-1.0.2rc1`.

## Regras de segurança da base

- Nunca diminua ou reutilize `database_revision`.
- Nunca altere `schema_version` sem atualizar também o código do aplicativo.
- Não edite o `sha256` manualmente.
- Sempre reconstrua o manifest depois da última alteração.
- Sempre publique `synonyms_database.json` e `synonyms_manifest.json` juntos.
- Não execute `scripts/build_synonym_database.py` para atualizações rotineiras;
  ele existe apenas para reconstruir a base inicial a partir dos mapas antigos.

## Resumo curto

```text
Abrir Codespaces
→ editar synonyms_database.json
→ aumentar database_revision e dictionary_version
→ executar rebuild_synonym_manifest.py
→ conferir git diff
→ commit e push dos dois JSONs
→ testar Check for updates no Pynnotate
```
