def create_directory_safely(directory_path, root=None):
	import os
	
	try:
		os.makedirs(directory_path, exist_ok=True)
		
		test_file = os.path.join(directory_path, '.write_test')
		try:
			with open(test_file, 'w') as f:
				f.write('test')
			os.remove(test_file)
			return True
		except PermissionError:
			report_message(f"Permission denied to write in directory: {directory_path}. Try running as administrator or choose a different folder.", "error", "Permission Error", root)
			return False
		except Exception as e:
			report_message(f"Cannot write to directory {directory_path}: {e}", "error", "Directory Error", root)
			return False
			
	except PermissionError:
		report_message(f"Permission denied to create directory: {directory_path}. Try running as administrator or choose a different folder.", "error", "Permission Error", root)
		return False
	except FileExistsError:
		return True
	except Exception as e:
		report_message(f"Error creating directory {directory_path}: {e}", "error", "Directory Error", root)
		return False

def load_synonyms_from_str(input_str, root=None):
	import yaml

	try:
		syn_dict = yaml.safe_load(input_str)
		if not isinstance(syn_dict, dict):
			report_message("Parsed synonyms must be a dictionary.", "error", "Error", root)
			return {}
		return syn_dict
	except Exception as e:
		report_message(f"Failed to parse synonyms: {e}", "error", "Error", root)
		return {}

def merge_synonyms(base_dict, user_dict):
	merged = base_dict.copy()
	for key, synonyms in user_dict.items():
		if isinstance(synonyms, str):
			synonyms = [synonyms]
		elif not isinstance(synonyms, list):
			synonyms = list(synonyms) if synonyms else []
		
		if key in merged:
			merged[key] = list(dict.fromkeys(merged[key] + synonyms))
		else:
			merged[key] = synonyms
	return merged


def normalize(text):
	import unicodedata
	return unicodedata.normalize("NFKC", text).strip().upper()

def standardize_organism(name):

	import re
	return re.sub(r'\s+', ' ', (name or '').strip())

def search_genbank(email, query, retmax, usehistory=None, download=None, api_key=None):
	
	from Bio import Entrez

	Entrez.email = email
	Entrez.api_key = api_key or None
	with Entrez.esearch(db="nucleotide", term=query, retmax=retmax, usehistory=usehistory) as search_handle:

		search_results = Entrez.read(search_handle)

		if download:
			count = int(search_results["Count"])
			webenv = search_results["WebEnv"]
			query_key = search_results["QueryKey"]
			return count, webenv, query_key
		else: 
			return search_results["IdList"]

def count_genes(rec):
	genes_found = set()
	for feature in rec.features:
		if feature.type in ["gene", "CDS", "rRNA", "tRNA"]:
			for key in ["gene", "product", "note"]:
				if key in feature.qualifiers:
					genes_found.add(feature.qualifiers[key][0].strip().upper())
					break
	return len(genes_found)

def filter_ind_by_spp(all_records, alias_map, genes_selected, removal_reasons, unique, prioritize):
	
	from collections import defaultdict

	filtered_records = []
	groups_by_spp = defaultdict(list)

	for record in all_records:
		organism, identity_source, ident = get_individual_identity(record)
		record.annotations["pynnotate_individual_id"] = ident
		record.annotations["pynnotate_individual_id_source"] = identity_source
		groups_by_spp[(organism, identity_source, ident)].append(record)

	candidates_by_spp = defaultdict(list)
	for (organism, identity_source, ident), records in groups_by_spp.items():
		candidates_by_spp[organism].append(((identity_source, ident), records))

	for spp, candidates in candidates_by_spp.items():
		ranked_candidates = []
		for ident, records in candidates:
			genes_temp = set()
			seen_temp = set()
			dummy_dict = {}
			dummy_seen_leu_ser = set()
			dummy_leu = defaultdict(list)
			dummy_ser = defaultdict(list)
			for r in records:
				r.annotations["original_id"] = r.id.split(".")[0]
				extract_genes(
					r, "dummy", alias_map, genes_selected, seen_temp,
					dummy_dict, set(), set(), dummy_leu, dummy_ser, dummy_seen_leu_ser, unique, prioritize
				)
			genes_temp.update(dummy_dict.keys())
			total_bp = sum(
				max((len(sequence) for _, sequence, _ in entries), default=0)
				for entries in dummy_dict.values()
			)
			accession_key = tuple(sorted(r.id.split(".")[0] for r in records))
			ranked_candidates.append(
				(ident, records, len(genes_temp), total_bp, accession_key)
			)

		if not ranked_candidates:
			continue

		ranked_candidates.sort(
			key=lambda item: (
				-item[2],
				-item[3],
				item[4],
				str(item[0][0]).casefold(),
				str(item[0][1]).casefold(),
			)
		)
		best_id, best, best_gene_count, best_total_bp, _ = ranked_candidates[0]

		filtered_records.extend(best)

		for ident, records, gene_count, total_bp, _ in ranked_candidates[1:]:
			if gene_count < best_gene_count:
				reason = "Removed by strict prioritization: another individual contained more recognized target genes"
			elif total_bp < best_total_bp:
				reason = "Removed by strict prioritization: gene count tied, but another individual had a greater cumulative target-gene length"
			else:
				reason = "Removed by strict prioritization: gene count and cumulative length tied; stable GenBank-accession tie-breaker applied"
			if ident != best_id:
				for r in records:
					removal_reasons[r.id].append(reason)

	return filtered_records

def filter_largest_fragment(gene_dict, header_to_species=None):

	from collections import defaultdict

	new_dict = {}
	for gene, entries in gene_dict.items():
		by_key = defaultdict(list)
		for header, seq, original_id in entries:
			key = (
				header_to_species.get(header)
				if header_to_species
				else header
			)
			if key:
				by_key[key].append((header, seq, original_id))

		new_dict[gene] = []
		for group in by_key.values():
			largest = max(group, key=lambda x: len(x[1]))
			new_dict[gene].append(largest)

	return new_dict

def fetch_with_retry(ids=None, rettype="gb", retmode="text", **kwargs):

	from http.client import IncompleteRead
	from Bio import Entrez, SeqIO
	import time

	last_error = None
	for try_n in range(3):
		try:
			if ids:
				with Entrez.efetch(db="nucleotide", id=ids, rettype=rettype, retmode=retmode, **kwargs) as handle:
					return list(SeqIO.parse(handle, "genbank"))
			else:
				with Entrez.efetch(db="nucleotide", rettype=rettype, retmode=retmode, **kwargs) as handle:
					return list(SeqIO.parse(handle, "genbank"))
		except IncompleteRead as e:
			last_error = e
			report_message(f"Incomplete response from GenBank, try {try_n+1}/3...", "warning", "Warning")
		except Exception as e:
			last_error = e
			report_message(f"GenBank request failed, try {try_n+1}/3: {e}", "warning", "Warning")
		if try_n < 2:
			time.sleep(2 * (try_n + 1))
	raise RuntimeError(f"Unable to retrieve a GenBank batch after 3 attempts: {last_error}")
	
def get_standardized_organism(record):

	organism = record.annotations.get("organism", "")
	if not organism:
		for feature in record.features:
			if feature.type == "source" and "organism" in feature.qualifiers:
				organism = feature.qualifiers["organism"][0]
				break
	return standardize_organism(organism or "Unknown")

def classify_positional_tRNA(start, ref_pos_all, name_trna, tolerancia_bp=10):
	
	anchors = sorted(
		[(gene, pos) for gene, pos_list in ref_pos_all.items()
		 if not gene.upper().startswith("TRNA") for pos in pos_list],
		key=lambda x: x[1]
	)

	for gene, pos in anchors:
		if start < (pos + tolerancia_bp):
			return f"{name_trna}_before_start_{gene}"

	return f"{name_trna}_unclassified"

def extract_genes(record, header_id, alias_map, genes_selected, seen_genes, gene_dict, aliases_unknown, aliases_duplicated, grouper_leu, grouper_ser, seen_leu_ser, unique, prioritize, species_to_genes=None, organism=None):

	from collections import defaultdict

	original_id = record.annotations.get("original_id", record.id)

	alias_lookup = {}
	for main_name, aliases in alias_map.items():
		alias_lookup[main_name.upper()] = main_name
		for alias in aliases:
				alias_lookup[alias.upper()] = main_name

	genes_selected_set = {g.upper() for g in genes_selected} if genes_selected else None

	ref_pos = {}
	ref_pos_all = {}
	for feature in record.features:
		try:
			if feature.type not in ["gene", "CDS", "rRNA", "tRNA"]:
				continue

			raw_name = None
			for key in ["gene", "product", "note"]:
				if key in feature.qualifiers:
					raw_name = feature.qualifiers[key][0].strip()
					break
				if not raw_name:
					continue

			raw_name_upper = raw_name.upper()
			gene = None
			for main_name, aliases in alias_map.items():
				all_names = [main_name.upper()] + [a.upper() for a in aliases]
				if raw_name_upper in all_names:
					gene = main_name
					break
				
			if not gene:
				continue

			pos = int(feature.location.start)
			if gene not in ref_pos_all:
				ref_pos_all[gene] = []
			ref_pos_all[gene].append(pos)

		except:
			pass

	for feature in record.features:
		try:
			if feature.type not in ["gene", "CDS", "rRNA", "tRNA"]:
				continue

			raw_name = None
			for key in ["gene", "product", "note"]:
				if key in feature.qualifiers:
					raw_name = feature.qualifiers[key][0]
					break
			if not raw_name:
				continue

			raw_name = normalize(raw_name)
			raw_upper = raw_name.upper()
			name = alias_lookup.get(raw_upper)
			if not name:
				aliases_unknown[raw_name].add(record.id)
				continue

			if genes_selected_set and name.upper() not in genes_selected_set:
				continue

			if unique:
				if species_to_genes and organism:
					if name.upper() in species_to_genes.get(organism, set()):
						continue

			start = int(feature.location.start)
			end = int(feature.location.end)
			strand = feature.location.strand
			seq = record.seq[start:end]
			if strand == -1:
				seq = seq.reverse_complement()

			if name in ["tRNA-Leu", "tRNA-Ser"]:
				tag = classify_positional_tRNA(start, ref_pos_all, name)

				if (original_id, tag) in seen_genes:
					continue
				
				if unique:
					if species_to_genes and organism:
						if tag in species_to_genes.get(organism, set()):
							continue

						gene_name = tag.split("_")[0]
						tags_for_same_gene = [t for t in gene_dict if t.startswith(gene_name)]

						tags_with_header = set()
						for t in tags_for_same_gene:
							for h, _, oid in gene_dict[t]:
								if h == header_id:
									tags_with_header.add(t)

						# Agora permite até 2 tags diferentes com esse header_id
						if len(tags_with_header) >= 2 and tag not in tags_with_header:
							# header_id já está em pelo menos 2 tags diferentes e não está na tag atual,
							# então ignora essa entrada
							continue

				if tag not in gene_dict:
					gene_dict[tag] = []
				gene_dict[tag].append((header_id, str(seq), original_id))
				seen_genes.add((original_id, tag))

				key = (header_id, name, tag)
				if key not in seen_leu_ser:
					if name.upper() == "TRNA-LEU":
						grouper_leu[tag].append((header_id, str(seq), original_id))
					else:
						grouper_ser[tag].append((header_id, str(seq), original_id))
					seen_leu_ser.add(key)
				continue

			if (original_id, name) in seen_genes:
				aliases_duplicated.add(name)
				continue

			if name not in gene_dict:
				gene_dict[name] = []
			gene_dict[name].append((header_id, str(seq), original_id))
			seen_genes.add((original_id, name))
		except:
			pass

	#gene_dict.pop("tRNA-Leu", None)
	#gene_dict.pop("tRNA-Ser", None)

def get_gene_extraction_failure_reason(record, alias_map, genes_selected):

	feature_types = {"gene", "CDS", "rRNA", "tRNA"}
	annotated_features = []
	annotated_names = []
	for feature in record.features:
		if feature.type not in feature_types:
			continue
		annotated_features.append(feature)
		for key in ["gene", "product", "note"]:
			values = feature.qualifiers.get(key, [])
			if values and str(values[0]).strip():
				annotated_names.append(str(values[0]).strip())
				break

	if not annotated_features or not annotated_names:
		return "Gene extraction failed: no annotated gene features were found in the GenBank record"

	alias_lookup = {}
	for main_name, aliases in alias_map.items():
		alias_lookup[normalize(main_name)] = main_name
		for alias in aliases:
			alias_lookup[normalize(alias)] = main_name

	recognized_genes = {
		alias_lookup[normalize(raw_name)]
		for raw_name in annotated_names
		if normalize(raw_name) in alias_lookup
	}
	
	if not recognized_genes:
		unrecognized_names = sorted(list(set(annotated_names)))
		unrecognized_str = ", ".join(f"'{name}'" for name in unrecognized_names)
		return f"Gene extraction failed: annotated gene names were not recognized by the selected synonym dictionary (found: {unrecognized_str})"

	selected_genes = {normalize(gene) for gene in genes_selected} if genes_selected else set()
	if selected_genes and not any(
		normalize(gene) in selected_genes for gene in recognized_genes
	):
		return "No requested genes extracted: annotated genes in this record do not match the selected target genes"

	return "Gene extraction failed: requested genes were recognized, but no sequence could be extracted"

def get_individual_identity(record):

	organism = get_standardized_organism(record)
	voucher = ""
	isolate = ""
	for feature in record.features:
		if feature.type == "source":
			voucher = str(feature.qualifiers.get("specimen_voucher", [""])[0]).strip()
			isolate = str(feature.qualifiers.get("isolate", [""])[0]).strip()
			break

	if voucher:
		return organism, "specimen_voucher", voucher
	if isolate:
		return organism, "isolate", isolate
	return organism, "genbank_id", record.id.split(".")[0]

def get_record_gene_profile(record, alias_map, genes_selected):

	from collections import defaultdict

	original_id = record.id.split(".")[0]
	record.annotations["original_id"] = original_id
	dummy_dict = {}
	extract_genes(
		record,
		f"profile_{original_id}",
		alias_map,
		genes_selected,
		set(),
		dummy_dict,
		defaultdict(set),
		set(),
		defaultdict(list),
		defaultdict(list),
		set(),
		False,
		False,
	)
	genes = set(dummy_dict.keys())
	total_bp = sum(
		max((len(sequence) for _, sequence, _ in entries), default=0)
		for entries in dummy_dict.values()
	)
	return genes, total_bp

def order_records_for_flexible(records, alias_map, genes_selected):

	from collections import defaultdict
	import re

	def normalized(value):
		return re.sub(r"\s+", " ", str(value or "").strip()).casefold()

	profiles = {}
	groups_by_species = defaultdict(lambda: defaultdict(list))
	for record in records:
		organism, identity_source, identity_value = get_individual_identity(record)
		record.annotations["pynnotate_individual_id"] = identity_value
		record.annotations["pynnotate_individual_id_source"] = identity_source
		profiles[id(record)] = get_record_gene_profile(record, alias_map, genes_selected)
		individual_key = (identity_source, normalized(identity_value))
		groups_by_species[organism][individual_key].append(record)

	ordered_records = []
	for organism in sorted(groups_by_species, key=normalized):
		individuals = []
		for individual_key, individual_records in groups_by_species[organism].items():
			genes = set()
			total_bp = 0
			for record in individual_records:
				record_genes, record_bp = profiles[id(record)]
				genes.update(record_genes)
				total_bp += record_bp
			individual_records.sort(
				key=lambda record: (
					-len(profiles[id(record)][0]),
					-profiles[id(record)][1],
					record.id,
				)
			)
			individuals.append((individual_key, genes, total_bp, individual_records))

		individuals.sort(
			key=lambda item: (-len(item[1]), -item[2], item[0][0], item[0][1])
		)
		for _, _, _, individual_records in individuals:
			ordered_records.extend(individual_records)

	return ordered_records

def assign_flexible_taxon_labels(records, header_fields=None):

	import re

	def component(value):
		return re.sub(r"\s+", "_", str(value or "").strip()).strip("_")

	for record in records:
		genbank_id = record.id.split(".")[0]
		organism = get_standardized_organism(record)
		if organism and organism.casefold() != "unknown":
			label = component(organism)
			label_source = "organism"
		else:
			label = genbank_id
			label_source = "genbank_id"
		record.annotations["pynnotate_taxon_label"] = label
		record.annotations["pynnotate_taxon_label_source"] = label_source

def overlap_animals(record, gene_dict, alias_map):

	original_id = record.annotations.get("original_id", record.id)
	organism_name = record.annotations.get("organism", "Unknown")
	header_id = f"{original_id}"
	alias_lookup = {}
	for main_name, aliases in alias_map.items():
		alias_lookup[normalize(main_name)] = main_name
		for alias in aliases:
			alias_lookup[normalize(alias)] = main_name

	gene_coords = []
	for feature in record.features:
		if feature.type not in ["gene", "CDS", "rRNA"]:
			continue

		raw_name = None
		for key in ["gene", "product", "note"]:
			if key in feature.qualifiers:
				raw_name = feature.qualifiers[key][0]
				break

		if not raw_name:
			continue

		name = alias_lookup.get(normalize(raw_name))
		if not name:
			continue

		start = int(feature.location.start)
		end = int(feature.location.end)
		strand = feature.location.strand
		gene_coords.append((name, start, end, strand, feature))

	gene_coords.sort(key=lambda x: x[1])
	log_records = []

	for i in range(len(gene_coords) - 1):
		name, start, end, strand, _ = gene_coords[i]
		next_name, next_start, next_end, _, _ = gene_coords[i + 1]

		if name != next_name and end > next_start:
			overlap = end - next_start
			log_records.append((name, start, end, next_name, next_start, next_end, overlap))

			new_entries = []
			for header, seq, orig_id in gene_dict.get(name, []):
				if orig_id != original_id:
					new_entries.append((header, seq, orig_id))
					continue

				if strand == 1:
					seq_corrected = seq[:next_start - start]
				else:
					cut = end - next_start
					seq_corrected = seq[cut:]

				new_entries.append((header, seq_corrected, orig_id))

			gene_dict[name] = new_entries

	return header_id, organism_name, log_records
		
def metadata_by_gene(gene_dict, excel_data, fields_excel, output_dir, all_records, extract_all_genes):
	
	from openpyxl import Workbook
	import os
	
	wb = Workbook()

	ws = wb.active
	ws.title = "Metadata"

	header_map = dict(excel_data)

	ws.append(["GenBank_ID"] + fields_excel)
	excel_data.sort(key=lambda x: x[1].get("organism", "").lower())
	for rid, info in excel_data:
		ws.append([rid] + [info.get(f, "") for f in fields_excel])

	if extract_all_genes and len(gene_dict) > 1 and len(all_records) > 1:
		for gene, entries in sorted(gene_dict.items()):
			ws = wb.create_sheet(title=gene[:31])  # Excel limits tab name
			ws.append(["GenBank_ID"] + fields_excel)
				
			entries = sorted(entries, key=lambda x: header_map.get(x[2], {}).get("organism", "").lower())
				 
			for head, _, original_id in entries:
				meta = header_map.get(original_id, {}).copy()
				if "organism" in meta:
					meta["organism"] = standardize_organism(meta["organism"])
				ws.append([original_id] + [meta.get(f, "") for f in fields_excel])

	out_path = os.path.join(output_dir, "metadata.xlsx")
	wb.save(out_path)

def metadata_gene_vs_sample(gene_dict, excel_data, header_fields, output_dir, config):

	from openpyxl import Workbook
	from openpyxl.styles import Font, PatternFill
	from openpyxl.utils import get_column_letter
	import os
	import re

	header_map = dict(excel_data)
	genes = sorted(gene_dict.keys())
	flexible_mode = config.get("unique_species", False)
	show_internal_identity = not flexible_mode and not any(
		field in header_fields for field in ("specimen_voucher", "isolate")
	)
	matrix_fields = ["organism"] if flexible_mode else list(header_fields)
	output_metadata_fields = list(matrix_fields)
	if show_internal_identity:
		output_metadata_fields.extend(["individual_id", "individual_id_source"])
	samples = {}

	def normalize_identifier(value):
		return re.sub(r"\s+", " ", str(value or "").strip()).casefold()

	def add_unique(values, value):
		value = str(value or "").strip()
		if value and value not in values:
			values.append(value)

	for gene, entries in gene_dict.items():
		for _, _, genbank_id in entries:
			meta = header_map.get(genbank_id, {})
			organism = standardize_organism(meta.get("organism", ""))
			if flexible_mode:
				identifier_value = str(
					meta.get("_taxon_label", "") or genbank_id
				).strip()
				identifier_field = str(
					meta.get("_taxon_label_source", "") or "organism"
				).strip()
				key = ("species", normalize_identifier(identifier_value))
			else:
				identifier_field = str(
					meta.get("_individual_id_source", "") or ""
				).strip()
				identifier_value = str(
					meta.get("_individual_id", "") or ""
				).strip()

				if not identifier_value:
					voucher = str(meta.get("specimen_voucher", "") or "").strip()
					isolate = str(meta.get("isolate", "") or "").strip()
					if voucher:
						identifier_field, identifier_value = "specimen_voucher", voucher
					elif isolate:
						identifier_field, identifier_value = "isolate", isolate
					else:
						identifier_field, identifier_value = "genbank_id", genbank_id

				key = (
					"individual",
					normalize_identifier(organism),
					identifier_field,
					normalize_identifier(identifier_value),
				)

			if key not in samples:
				samples[key] = {
					"metadata": {field: [] for field in matrix_fields},
					"genes": {selected_gene: [] for selected_gene in genes},
					"organism": organism,
					"identifier": identifier_value or genbank_id,
					"identifier_source": identifier_field,
				}

			sample = samples[key]
			for field in matrix_fields:
				if flexible_mode and field == "organism":
					value = organism or identifier_value
				else:
					value = organism if field == "organism" else meta.get(field, "")
				add_unique(sample["metadata"][field], value)
			add_unique(sample["genes"][gene], genbank_id)

	wb = Workbook()
	ws = wb.active
	ws.title = "metadata_gene_vs_sample"
	ws.append(output_metadata_fields + genes)

	ordered = sorted(
		samples.values(),
		key=lambda sample: (
			normalize_identifier(sample["organism"]),
			normalize_identifier(sample["identifier"]),
		),
	)

	conflicts = []
	for sample in ordered:
		row = ["; ".join(sample["metadata"][field]) for field in matrix_fields]
		if show_internal_identity:
			row.extend([sample["identifier"], sample["identifier_source"]])
		for gene in genes:
			accessions = sample["genes"][gene]
			row.append("; ".join(accessions))
			if len(accessions) > 1:
				conflicts.append({
					"individual": sample["identifier"],
					"organism": sample["organism"],
					"gene": gene,
					"accessions": accessions.copy(),
				})
		ws.append(row)

	def format_data_sheet(sheet, widths):
		header_fill = PatternFill("solid", fgColor="D9EAD3")
		for cell in sheet[1]:
			cell.font = Font(bold=True)
			cell.fill = header_fill
		sheet.freeze_panes = "B2" if sheet.max_column > 1 else "A2"
		sheet.auto_filter.ref = sheet.dimensions
		for column, width in widths.items():
			sheet.column_dimensions[column].width = width

	main_widths = {"A": 34}
	for column_index in range(2, ws.max_column + 1):
		main_widths[get_column_letter(column_index)] = 16
	format_data_sheet(ws, main_widths)

	if flexible_mode:
		provenance_ws = wb.create_sheet(title="gene_provenance")
		provenance_ws.append([
			"organism", "gene", "GenBank_ID", "isolate", "specimen_voucher"
		])
		provenance_rows = []
		for gene in genes:
			for _, _, genbank_id in gene_dict.get(gene, []):
				meta = header_map.get(genbank_id, {})
				organism = standardize_organism(meta.get("organism", ""))
				if not organism:
					organism = str(meta.get("_taxon_label", "") or genbank_id)
				provenance_rows.append([
					organism,
					gene,
					genbank_id,
					meta.get("isolate", ""),
					meta.get("specimen_voucher", ""),
				])
		for provenance_row in sorted(
			provenance_rows,
			key=lambda row: (
				normalize_identifier(row[0]),
				normalize_identifier(row[1]),
				normalize_identifier(row[2]),
			),
		):
			provenance_ws.append(provenance_row)
		format_data_sheet(
			provenance_ws,
			{"A": 34, "B": 16, "C": 16, "D": 22, "E": 28},
		)

	out_path = os.path.join(output_dir, "genes_matrix.xlsx")
	wb.save(out_path)
	return conflicts

def make_query(email, alias_map, genes, organisms, mito, mitogenome, chloroplast, title, additional, min_len, max_len, add_unverified_exclusion, api_key=None):

	import logging

	# Generates gene synonyms
	gene_query = ""

	if isinstance(genes, list):
		selected_genes = [g.strip() for g in genes if g.strip()]
	elif isinstance(genes, str):
		selected_genes = [g.strip() for g in genes.split(",") if g.strip()]
	else:
		selected_genes = []
	grouped_aliases = []
	for g in selected_genes:
		match_found = False
		for main_name, aliases in alias_map.items():
			all_names = [main_name, *aliases]
			if normalize(g) in {normalize(name) for name in all_names}:
				unique_names = list(dict.fromkeys(all_names))
				group = " OR ".join(f'"{a}"' if " " in a else a for a in unique_names)
				grouped_aliases.append(f"({group})")
				match_found = True
				break
		if not match_found:
			grouped_aliases.append(f'"{g}"' if " " in g else g)
			
		if grouped_aliases:
			gene_query = " OR ".join(grouped_aliases)
		
	query_parts = []
	if gene_query:
		query_parts.append(f"({gene_query})")

	if organisms:
		if isinstance(organisms, str):
			organism_terms = [o.strip() for o in organisms.split(",") if o.strip()]
		else:
			organism_terms = [str(o).strip() for o in organisms if str(o).strip()]
		org_query = " OR ".join([f'"{o}"[Organism]' for o in organism_terms])
		query_parts.append(f"({org_query})")

	if mito:
		query_parts.append("mitochondrial")

	if mitogenome:
		query_parts.append("mitochondrion[Filter]")

	if chloroplast:
		query_parts.append("chloroplast[Filter]")

	if title:
		query_parts.append(title)
	 
	if min_len or max_len:
		min_val = int(min_len) if min_len else ""
		max_val = int(max_len) if max_len else ""
			
		if min_val and max_val:
			query_parts.append(f"{min_val}:{max_val}[SLEN]")
		elif min_val:
			query_parts.append(f"{min_val}:99999999[SLEN]")
		elif max_val:
			query_parts.append(f"2:{max_val}[SLEN]")

	query = " AND ".join(query_parts)
	query += " NOT wgs[prop] NOT tsa[prop]"
	if additional:
		query += f" {additional}"
	
	if add_unverified_exclusion:
		query += " NOT UNVERIFIED"

	try:
		ids = search_genbank(email, query, retmax=999999, api_key=api_key)
	except Exception as e:
		try:
			import tkinter
			from tkinter import messagebox
			root = tkinter._default_root
			if root is not None:
				report_message(f"An error occurred: {e}", "error", "Error", root)
			else:
				raise RuntimeError
		except:
			report_message(f"An error occurred: {e}", "error", "Error")
			ids = []

	return ids, query

def report_message(msg, level="info", title=None, root=None):

	import logging

	title = title or level.capitalize()
	if root is not None:
		try:
			import tkinter.messagebox as messagebox

			def show_dialog():
				if level == "error":
					messagebox.showerror(title, msg, parent=root)
				elif level == "warning":
					messagebox.showwarning(title, msg, parent=root)
				else:
					messagebox.showinfo(title, msg, parent=root)

			root.after(0, show_dialog)
			return
		except Exception:
			pass 
	
	if level == "error":
		logging.error(f"{title}: {msg}")
	elif level == "warning":
		logging.warning(f"{title}: {msg}")
	else:
		logging.info(f"{title}: {msg}")

def sanitize_filename(filename):

	import re
	import platform
	
	if not filename:
		return "unnamed"
	
	invalid_chars = r'[<>:"|?*\\\/]'
	filename = re.sub(invalid_chars, '_', filename)
	
	if platform.system() == 'Windows':
		reserved = {'CON', 'PRN', 'AUX', 'NUL', 'COM1', 'COM2', 'COM3', 
				   'COM4', 'COM5', 'COM6', 'COM7', 'COM8', 'COM9',
				   'LPT1', 'LPT2', 'LPT3', 'LPT4', 'LPT5', 'LPT6', 
				   'LPT7', 'LPT8', 'LPT9'}
		
		name_only = filename.split('.')[0].upper()
		if name_only in reserved:
			filename = f"_{filename}"
	
	if len(filename) > 200:
		name, ext = filename.rsplit('.', 1) if '.' in filename else (filename, '')
		filename = name[:200-len(ext)-1] + ('.' + ext if ext else '')
	
	return filename.strip()

def begin_search(config, root=None):

	from datetime import datetime
	import threading
	import logging
	import os
	import re

	logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

	from pynnotate.main.synonym_updates import export_effective_synonyms

	report_message("WELCOME TO PYNNOTATE!\n", "info", "Info")

	email = config.get("email", "")
	title = config.get("title", "")
	additional = config.get("additional", "")
	organisms = config.get("organisms", "")
	genes = config.get("genes", "")
	fields = config.get("fields", "")
	output_name = config.get("output_name", "")
	folder_name = config.get("folder_name", "")
	if folder_name:
		folder_name = sanitize_filename(folder_name)
	include_id = config.get("include_id", True)
	mito = config.get("mito", False)
	mitogenome = config.get("mitogenome", False)
	ids_text = config.get("ids_text", "")
	chloroplast = config.get("chloroplast", False)
	min_len = config.get("min_len", "")
	max_len = config.get("max_len", "")
	delete_unverified = config.get("delete_unverified", False)
	extraction = config.get("extraction", False)
	overlap = config.get("overlap", False)
	unique = config.get("unique", False)
	prioritize = config.get("prioritize_more_genes", False)
	logmissing = config.get("logmissing", False)
	org_type = config.get("org_type", "animal_mito")
	api_key = config.get("api_key", "")
	add_synonyms = config.get("add_synonyms", "")

	if unique and not extraction:
		report_message(
			"Flexible mode requires 'Separate gene extraction' to identify and compare genes.",
			"warning", "Required option", root
		)
		return

	if not extraction:
		if overlap:
			report_message("⚠️  '--overlap' ignored because '--extraction' was not enabled.", "warning", "Warning")
		if logmissing:
			report_message("⚠️  '--logmissing' ignored because '--extraction' was not enabled.", "warning", "Warning")
		overlap = False
		logmissing = False
	elif org_type != "animal_mito":
		if overlap:
			report_message("⚠️  '--overlap' is only applicable to 'animal_mito'. Ignoring.", "warning", "Warning")
		overlap = False

	if add_synonyms:
		try:
			if isinstance(add_synonyms, dict):
				user_synonyms = add_synonyms
			elif isinstance(add_synonyms, str):
				user_synonyms = load_synonyms_from_str(add_synonyms, root)
			else:
				raise ValueError("add_synonyms must be a string or a dictionary.")
		except ValueError as e:
			report_message(f"Failed to parse synonyms: {e}", "error", "Error", root)
			return
	else:
		user_synonyms = {}

	try:
		synonym_snapshot = export_effective_synonyms(org_type, user_synonyms)
	except Exception as e:
		report_message(f"Failed to load the synonym database: {e}", "error", "Error", root)
		return
	alias_map = synonym_snapshot["synonyms"]

	valid_fields = ["organism", "isolate", "specimen_voucher", "geo_loc_name",
		"country", "collection_date", "host", "strain", "isolation_source",
		"lat_lon"]

	add_unverified_deletion = delete_unverified

	if not fields or not output_name or not email:
		report_message("Fill in the required fields: email, output, and header", "warning", "Required fields", root)
		return

	if isinstance(fields, list):
		header_fields = [i.strip() for i in fields if i.strip()]
	elif isinstance(fields, str):
		header_fields = [i.strip() for i in re.split(r"[\s,]+", fields) if i.strip()]
	else:
		header_fields = []

	invalid_fields = [f for f in header_fields if f not in valid_fields]
	
	if invalid_fields:
		report_message(f"Invalid fields: {', '.join(invalid_fields)}", "error", "Error", root)
		return

	if min_len or max_len:
		try:
			min_val = int(min_len) if min_len else ""
			max_val = int(max_len) if max_len else ""
		except ValueError:
			report_message("Min/max size values must be integers.", "error", "Error", root)
			return

	date_string = datetime.now().strftime("%d%b%y").lower()
	time_string = datetime.now().strftime("%Hh%Mm%Ss").lower()
	folder = output_name if isinstance(output_name, str) else output_name[0]

	if folder_name:
		output_dir = os.path.join(folder, folder_name)
	else:
		output_dir = os.path.join(folder, f"pynnotate_{date_string}_{time_string}")

	def launch_download(download_config, download_genes, download_query):
		ui_state = None
		if root is not None:
			from pynnotate.graphic.utils import create_loading_window_with_progress
			ui_state = create_loading_window_with_progress(root)

		worker = threading.Thread(
			target=download_sequences,
			args=(
				email, ids, output_dir, header_fields, alias_map,
				download_config, download_genes, download_query, root, ui_state,
			),
			name="pynnotate-download",
		)
		worker.start()

	if ids_text:
		if isinstance(ids_text, str):
			ids = [i.strip() for i in re.split(r"[\s,]+", ids_text) if i.strip()]
		else:
			ids = [str(i).strip() for i in ids_text if str(i).strip()]

		if root is not None:
			from tkinter import messagebox
			proceed_man = messagebox.askyesno("Confirmation", f"🔢 {len(ids)} IDs provided manually.\n\nDo you want to start the download?")
		else:
			response_man = input(f"🔢 {len(ids)} IDs provided manually.\n\nDo you want to start the download? (y/n): ")
			proceed_man = response_man.lower() in ["y", "yes"]

		if proceed_man:
			config = {
				"include_id": include_id,
				"extract_all_genes": extraction,
				"delete_overlap": overlap,
				"unique_species": unique,
				"prioritize_more_genes": prioritize,
				"min_len": int(min_len or 0),
				"max_len": int(max_len or 999999),
				"log_missing": logmissing,
				"org_type": org_type,
				"api_key": api_key,
				"add_synonyms": add_synonyms,
				"mitogenome": mitogenome,
				"chloroplast": chloroplast,
				"synonym_database": synonym_snapshot,
			}
			launch_download(config, "", None)
		return

	if not any([genes, organisms, title, additional, mito, mitogenome, chloroplast]):
		report_message("Invalid search: Please enter genes, organism, title or select at least one mitochondria/mitogenome option to perform the search.", "warning", "Warning", root)
		return

	ids, query = make_query(
		email, alias_map, genes, organisms, mito, mitogenome, chloroplast,
		title, additional, min_len, max_len, add_unverified_deletion, api_key,
	)

	if ids:
		if root is not None:
			from pynnotate.graphic.utils import confirm_query_download
			proceed_down = confirm_query_download(root, len(ids), query)
		else:
			response_down = input(f"🔢 {len(ids)} records found.\n\nDo you want to start downloading? (y/n): ")
			proceed_down = response_down.lower() in ["y", "yes"]
		if proceed_down:
			config = {
				"include_id": include_id,
				"extract_all_genes": extraction,
				"delete_overlap": overlap,
				"unique_species": unique,
				"prioritize_more_genes": prioritize,
				"min_len": int(min_len or 0),
				"max_len": int(max_len or 999999),
				"log_missing": logmissing,
				"org_type": org_type,
				"api_key": api_key,
				"add_synonyms": add_synonyms,
				"mitogenome": mitogenome,
				"chloroplast": chloroplast,
				"synonym_database": synonym_snapshot,
			}
			launch_download(config, genes, query)
	else:
		report_message("No results: No ID found for search.", "info", "Info", root)

def download_sequences(email, genbank_ids, output_dir, header_fields, alias_map, config, genes, query=None, root=None, ui_state=None):

	import os
	import json
	import time
	from tqdm import tqdm
	from Bio import Entrez, SeqIO
	from openpyxl import Workbook
	from collections import defaultdict
	import traceback

	Entrez.email = email
	Entrez.api_key = config.get("api_key") or None

	valid_fields = ["organism", "isolate", "specimen_voucher", "geo_loc_name",
		"country", "collection_date", "host", "strain", "isolation_source",
		"lat_lon"]

	genbank_ids = genbank_ids or []
	include_id = config["include_id"]
	extract_all_genes = config["extract_all_genes"]
	delete_overlap = config["delete_overlap"]
	unique_species = config.get("unique_species", False)
	prioritize_more_genes = config.get("prioritize_more_genes", False)
	log_missing = config.get("log_missing", False)
	org_type = config["org_type"]
	add_synonyms = config["add_synonyms"]
	mitogenome = config["mitogenome"]
	chloroplast = config["chloroplast"]

	if not genbank_ids and not query:
		report_message("No ID found.", "warning", "Warning", root)
		return
	
	if isinstance(genes, str):
		genes_list = [g.strip() for g in genes.split(",") if g.strip()]
	elif isinstance(genes, list):
		genes_list = [g.strip() for g in genes if isinstance(g, str) and g.strip()]
	else:
		genes_list = []

	genes_selected = set()
	gene_alias_lookup = {}
	for main_name, aliases in alias_map.items():
		for alias in [main_name, *aliases]:
			gene_alias_lookup[normalize(alias)] = main_name

	for g in genes_list:
		match = gene_alias_lookup.get(normalize(g))
		genes_selected.add(match if match else g)
	genes_selected_original = set(genes_selected)
		
	removal_reasons = defaultdict(list)

	if extract_all_genes and not genes_selected:
		if mitogenome:
			all_keys = set(alias_map.keys())
			keys_to_remove = {"18S", "28S", "ITS1", "ITS2", "RAG1", "RAG2", "BDNF", "CMOS", "POMC", "TYR", "RHO", "SIAH1", "CRYB", "TNS3", "S7", "TYRP1", "EIF3", "CXCR4", "NT3", "FIB7", "MYH6", "CMYC", "CMYC2", "H3A", "ZEB2", "MC1R", "SNCAIP", "PTPRC", "OPA1"}
			genes_selected = all_keys - keys_to_remove
			report_message("No genes specified, but mitogenome flag enabled — extracting all mitogenome genes", "info", "Info", root)
		else:
			genes_selected = set(alias_map.keys())
			report_message("No genes specified — extracting all genes", "info", "Info", root)
		genes_selected_original = set(genes_selected)

	unknown_aliases = defaultdict(set)
	duplicated_aliases = set()
	grouper_leu = defaultdict(list)
	grouper_ser = defaultdict(list)
	loading_win = progress_var = None
	timer_state = None
	download_started_at = time.monotonic()
	completion_message = None
	completion_level = "info"
	completion_title = "Success"
	pbar = None
	if ui_state is not None:
		loading_win, progress_var, _ = ui_state
		timer_state = getattr(loading_win, "_pynnotate_timer_state", None)
		if timer_state is not None:
			download_started_at = timer_state["started_at"]

	try:
		if not create_directory_safely(output_dir, root):
			return

		gene_folder = os.path.join(output_dir, "genes")
		fasta_file = os.path.join(output_dir, "sequences.fasta")
		excel_file = os.path.join(output_dir, "metadata.xlsx")
		log_file = os.path.join(output_dir, "log.txt")
		synonym_file = os.path.join(output_dir, "synonyms_used.json")
		synonym_snapshot = config.get("synonym_database")
		if not synonym_snapshot:
			from pynnotate.main.synonym_updates import build_synonym_snapshot
			synonym_snapshot = build_synonym_snapshot(org_type, alias_map, bool(add_synonyms))
		with open(synonym_file, "w", encoding="utf-8", newline="\n") as synonym_out:
			json.dump(synonym_snapshot, synonym_out, ensure_ascii=False, indent=2, sort_keys=True)
			synonym_out.write("\n")

		all_records = []
		batch_size = 500
		ignored_ids = []

		if query:

			count, webenv, query_key = search_genbank(
				email=email, query=query, retmax=0, usehistory="y",
				download=True, api_key=config.get("api_key"),
			)

			if root is None:
				pbar = tqdm(total=count, desc="Downloading data from GenBank... Please wait a few minutes", unit="seq")

			for start in range(0, count, batch_size):

				records = fetch_with_retry(
					retstart=start, retmax=batch_size,
					webenv=webenv, query_key=query_key
				)

				for rec in records:
					try:
						seq_string = str(rec.seq).upper()
						if not rec.seq or len(seq_string.strip("N")) == 0 or "SEQUENCE CONTENT IS UNDEFINED" in seq_string:
							ignored_ids.append(rec.id)
							continue
						all_records.append(rec)
					except Exception:
						ignored_ids.append(rec.id)
						continue
				if root is not None:
					progress = int(((start + batch_size) / count) * 100)
					if progress_var is not None:
						root.after(0, lambda value=min(progress, 100): progress_var.set(value))
				else:
					pbar.update(batch_size)

		else:

			if root is None:
				pbar = tqdm(total=len(genbank_ids), desc="Downloading data from GenBank... Please wait a few minutes", unit="seq")
			
			for i in range(0, len(genbank_ids), batch_size):
				chunk = genbank_ids[i:i + batch_size]
				records = fetch_with_retry(ids=chunk)
				for rec in records:
					try:
						seq_string = str(rec.seq).upper()
						if not rec.seq or len(seq_string.strip("N")) == 0 or "SEQUENCE CONTENT IS UNDEFINED" in seq_string:
							ignored_ids.append(rec.id)
							continue
						all_records.append(rec)
					except Exception:
						ignored_ids.append(rec.id)
						continue
				if root is not None:
					progress = int(((i + batch_size) / len(genbank_ids)) * 100)
					if progress_var is not None:
						root.after(0, lambda value=min(progress, 100): progress_var.set(value))
				else:
					pbar.update(batch_size)

		all_records_orig = all_records.copy()
		min_len = config.get("min_len", 0)
		max_len = config.get("max_len", 999999)

		all_records = [rec for rec in all_records if min_len <= len(rec.seq) <= max_len]

		if not all_records:
			report_message(f"No sequences kept after aplying length filter: {min_len}-{max_len}.", "error", "Error", root)
			return

		use_genes_priority = False

		if prioritize_more_genes:
			if not genes_list:
				use_genes_priority = True
			else:
				num_genes = len(genes_list)
				if num_genes > 1:
					use_genes_priority = True
	
		if use_genes_priority:
			all_records = filter_ind_by_spp(
				all_records, alias_map, genes_selected, removal_reasons, unique_species, prioritize_more_genes
			)

		if unique_species:
			all_records = order_records_for_flexible(
				all_records, alias_map, genes_selected
			)
			assign_flexible_taxon_labels(all_records)

		filtered_records = []
		species_seen = set()
		seen_headers = set()
		excel_data = []
		fields_excel = valid_fields
		gene_dict = {}
		original_id_to_species = {}
		ignored_aliases = set()
		all_log_overlap = []
		seen_genes = set()
		seen_headers = set()
		seen_leu_ser = set()
		taxon_label_fallbacks = []

		species_to_genes = defaultdict(set)

		for record in all_records:

			genbank_id = record.id.split(".")[0]
			record.annotations["original_id"] = genbank_id
			
			organism = record.annotations.get("organism", "")
			if not organism:
				for feature in record.features:
					if feature.type == "source" and "organism" in feature.qualifiers:
						organism = feature.qualifiers["organism"][0]
						break
			if not organism:
				organism = "Unknown"

			organism = standardize_organism(organism)
			_, identity_source, identity_value = get_individual_identity(record)
			record.annotations["pynnotate_individual_id"] = identity_value
			record.annotations["pynnotate_individual_id_source"] = identity_source

			current_record_genes, _ = get_record_gene_profile(
				record, alias_map, genes_selected
			)
			current_record_genes = {gene.upper() for gene in current_record_genes}

			if unique_species:
				new_genes = current_record_genes - species_to_genes[organism]
				if not new_genes:
					removal_reasons[genbank_id].append("All genes already saved for species (flexible mode⎯unique_species enabled)")
					continue

			if prioritize_more_genes:
				if organism in species_seen:
					removal_reasons[genbank_id].append("Another sequence for species already added (strict mode⎯prioritize enabled)")
					continue
			species_seen.add(organism)

			values = []
			feature_info = {f: "" for f in fields_excel}
			feature_info["_individual_id"] = identity_value
			feature_info["_individual_id_source"] = identity_source
			feature_info["_taxon_label"] = record.annotations.get(
				"pynnotate_taxon_label", ""
			)
			feature_info["_taxon_label_source"] = record.annotations.get(
				"pynnotate_taxon_label_source", ""
			)
			for field in header_fields:
				val = ""
				for feature in record.features:
					if feature.type == "source":
						val = feature.qualifiers.get(field, [""])[0].strip()
						break
				values.append(val.replace(" ", "_") if val else "")
		
			for feature in record.features:
				if feature.type == "source":
					for f in fields_excel:
						if f in feature.qualifiers:
							feature_info[f] = feature.qualifiers[f][0]
					break
		
			base_header = "_".join([v for v in values if v]).strip("_")
			if not base_header:
				base_header = organism.replace(" ", "_")
			header_id = f"{base_header}_{genbank_id}" if include_id else base_header
			original_header = header_id
				
			if not prioritize_more_genes:
				suffix = 1
				while header_id in seen_headers:
					header_id = f"{original_header}_{genbank_id if include_id else suffix}"
					suffix += 1
				seen_headers.add(header_id)
			else:
				seen_headers.add(header_id)
		
			if not include_id and genbank_id in header_id:
				header_id = header_id.replace(f"_{genbank_id}", "")
		
			record.id = header_id
			record.description = ""
			gene_header_id = (
				record.annotations.get("pynnotate_taxon_label", header_id)
				if unique_species else header_id
			)

			original_id_to_species[gene_header_id] = organism

			genes_found = set()

			overlap_header=None
			overlap_name=None
			overlap_records=[]

			if extract_all_genes:

				extract_genes(
					record,
					gene_header_id,
					alias_map,
					genes_selected,
					seen_genes,
					gene_dict,
					unknown_aliases,
					duplicated_aliases,
					grouper_leu,
					grouper_ser,
					seen_leu_ser,
					unique_species,
					prioritize_more_genes,
					species_to_genes,
					organism
				)
				
				if delete_overlap and org_type == "animal_mito":
					overlap_header, overlap_name, overlap_records = overlap_animals(record, gene_dict, alias_map)
					if overlap_records:
						all_log_overlap.append((overlap_header, overlap_name, overlap_records))

				for gene, seqs in gene_dict.items():
					for head, seq, seq_original_id in seqs:
						if seq_original_id == record.annotations["original_id"]:
							genes_found.add(gene.upper())
							break

			species_to_genes[organism].update(genes_found)

			if not extract_all_genes or genes_found:
				filtered_records.append(record)
				excel_data.append((record.annotations["original_id"], feature_info))
				if (
					unique_species
					and feature_info.get("_taxon_label_source") == "genbank_id"
				):
					taxon_label_fallbacks.append((genbank_id, header_fields))
			else:
				removal_reasons[genbank_id].append(
					get_gene_extraction_failure_reason(record, alias_map, genes_selected)
				)
				continue

		if prioritize_more_genes or unique_species:
			gene_dict = filter_largest_fragment(gene_dict)
		
		if unique_species:
			gene_dict = filter_largest_fragment(gene_dict, original_id_to_species)
		
		if extract_all_genes and not genes:
			genes_selected = set(gene_dict.keys())

		if filtered_records:
			with open(fasta_file, "w", encoding='utf-8', newline='\n') as f_out:
				SeqIO.write(filtered_records, f_out, "fasta")

		for tag, entries in grouper_leu.items():
			if not create_directory_safely(gene_folder, root):
				return
			path = os.path.join(gene_folder, f"{tag}.fasta")
			with open(path, "w", encoding='utf-8', newline='\n') as out:
				for header, seq, oid in entries:
					out.write(f">{header}\n{seq}\n")
		 
		for tag, entries in grouper_ser.items():
			if not create_directory_safely(gene_folder, root):
				return
			path = os.path.join(gene_folder, f"{tag}.fasta")
			with open(path, "w", encoding='utf-8', newline='\n') as out:
				for header, seq, oid in entries:
					out.write(f">{header}\n{seq}\n")

		with open(log_file, "w", encoding='utf-8', newline='\n') as lf:
			lf.write("PYNNOTATE\n")
			lf.write("\n")
			lf.write(f"User email: {email}\n")
			lf.write(f"Output folder: {output_dir}\n")
			lf.write("\n")
			lf.write(f"FASTA file generated: {fasta_file}\n")
			if extract_all_genes:
				lf.write(f"FASTA files separated by genes generated: {gene_folder}\n")
			lf.write(f"Generated Excel file: {excel_file}\n")
			lf.write("\n")
			if query:
				lf.write(f"Search query used: {query}\n")
			else:
				lf.write(f"GenBank IDS search: {genbank_ids}\n")
			lf.write(f"Genome type chosen: {org_type}\n")
			lf.write(f"Synonym database version: {synonym_snapshot.get('dictionary_version', 'unknown')}\n")
			lf.write(f"Synonym database revision: {synonym_snapshot.get('database_revision', 'unknown')}\n")
			lf.write(f"Synonym database source: {synonym_snapshot.get('database_source', 'unknown')}\n")
			lf.write(f"Synonym database SHA256: {synonym_snapshot.get('database_sha256', 'unknown')}\n")
			lf.write(f"Effective synonym SHA256: {synonym_snapshot.get('effective_sha256', 'unknown')}\n")
			lf.write(f"Synonym snapshot saved: {synonym_file}\n")
			if add_synonyms:
				lf.write("Custom synonyms supplied by the user: Yes\n")
			lf.write("\n")

			if not unique_species and not prioritize_more_genes:
				lf.write("[🌐 UNCONSTRAINED MODE: INCLUDING ALL SEQUENCES]\n")
			elif unique_species:
				lf.write("[🌱 FLEXIBLE MODE ENABLED: ALLOWING MULTIPLE SEQUENCES PER SPECIES IF THEY HAVE NEW GENES]\n")
				lf.write("Individual identity priority: specimen_voucher > isolate > GenBank ID (independent of FASTA header fields)\n")
			elif prioritize_more_genes:
				lf.write("[🔒 STRICT MODE: INCLUDING ONLY ONE SEQUENCE PER SPECIES, PRIORITIZING INDIVIDUALS WITH MORE GENES]\n")
			lf.write("\n")

			lf.write(f"Total sequences found (before filters): {len(all_records_orig)}\n")
			lf.write(f"Total sequences saved: {len(filtered_records)}\n")
			lf.write(f"Total species saved: {len(species_seen)}\n")
			individuals_saved = {
				(
					standardize_organism(meta.get("organism", "")),
					meta.get("_individual_id_source", "genbank_id"),
					meta.get("_individual_id", genbank_id),
				)
				for genbank_id, meta in excel_data
			}
			lf.write(f"Total individuals represented: {len(individuals_saved)}\n")
			lf.write(f"Fields used in the header: {', '.join(header_fields)}\n")
			lf.write(f"Include GenBank ID in header: {'Yes' if include_id else 'No'}\n")
			if unique_species:
				lf.write("Flexible gene FASTA headers use organism-only species labels; selected header fields and the GenBank-ID option apply only to sequences.fasta.\n")
				lf.write("genes_matrix.xlsx contains one row per organism plus a gene_provenance sheet with accession, isolate, and specimen voucher.\n")
			lf.write("\n")

			if taxon_label_fallbacks:
				lf.write("[HEADER FALLBACKS] Organism was unavailable; GenBank ID was used as the Flexible supermatrix label for these records:\n")
				for fallback_id, fallback_fields in taxon_label_fallbacks:
					lf.write(f"{fallback_id}\n")
				lf.write("\n")

			if extract_all_genes:
				lf.write(f"[🌟 SEPARATE GENE EXTRACTION ENABLED]: {len(gene_dict)} genes grouped by synonym\n")

			if extract_all_genes and gene_dict:
				lf.write("[EXTRACT GENES SUMMARY]:\n")
				for gene in sorted(gene_dict):
					entries = gene_dict[gene]
					if not create_directory_safely(gene_folder, root):
						return
					safe_gene_name = sanitize_filename(gene)
					path = os.path.join(gene_folder, f"{safe_gene_name}.fasta")
					with open(path, "w", encoding='utf-8', newline='\n') as out:
						for head, seq, *_ in entries:
							out.write(f">{head}\n{seq}\n")
					lf.write(f"{gene}: {len(entries)} extracted sequences\n")

				grouper_combined = list(grouper_leu.items()) + list(grouper_ser.items())
				lf.write("\n[TRNA-LEU/SER CLUSTERINGS BY POSITIONAL PATTERN]: When using flexible mode, it's expected that some species will have more than two copies of tRNA-Leu and tRNA-Ser. This happens because GenBank entries often lack a clear specification about which copy of Leu/Ser a sequence refers to. Due to this, pynnotate strategy is that genes are extracted based on their genomic position. Therefore, it's up to the researcher to manually inspect and process these cases afterwards.\n")
				
				tags_set = sorted({tag for tag, _ in grouper_combined})
				
				for tag in tags_set:
					count = 0
					for name in gene_dict:
						if name == tag:
							for header, seq, oid in gene_dict[name]:
								count += 1
					lf.write(f"{tag}: {count} grouped sequences\n")
				lf.write("\n")

			if delete_overlap and org_type == "animal_mito":
				if not all_log_overlap:
					lf.write("[FIXING OVERLAP BETWEEN EXTRACTED GENES ENABLED] No overlap found!\n")
				else:
					lf.write("[FIXING OVERLAP BETWEEN EXTRACTED GENES ENABLED] Here are the bp removed:\n")
					for overlap_header, overlap_name, overlap_records in all_log_overlap:
						key = f"{overlap_name} ({overlap_header})"
						lf.write(f"[{key}]\n")
						for g1, s1, e1, g2, s2, e2, o in overlap_records:
							lf.write(f"⚠️ {g1} ({s1}..{e1}) overlap {g2} ({s2}..{e2}): {o} bp removed from the end of {g1}\n")
				lf.write("\n")

			if unknown_aliases:
				lf.write("[UNKNOWN ALIASES] Some aliases were not detected in the synonyms dictionary and were not extracted:\n")
				for alias in sorted(unknown_aliases):
					origins = ", ".join(sorted(unknown_aliases[alias]))
					lf.write(f"{alias} (from records: {origins})\n")
				lf.write("\n")

			if duplicated_aliases:
				lf.write("[REPEATED GENE ANNOTATIONS] These canonical gene names appeared more than once in the same GenBank record. Later occurrences were ignored to keep one sequence per gene. This report does not distinguish overlapping feature annotations from distinct biological gene copies:\n")
				for alias in sorted(duplicated_aliases):
					lf.write(alias + "\n")
				lf.write("\n")

			if ignored_ids:
				lf.write("[IGNORED RECORDS] Records with missing or undefined sequence content:\n")
				for ignored_id in sorted(set(ignored_ids)):
					lf.write(f"{ignored_id}\n")
				lf.write("\n")

			if log_missing and genes_selected_original:
				lf.write("[MISSING GENES] Some records did not present the genes of interest:\n")
				for record in filtered_records:
					found_genes = set()
					record_id = record.annotations.get("original_id", record.id)
					for gene, seqs in gene_dict.items():
						for _, _, seq_original_id in seqs:
							if seq_original_id == record_id:
								found_genes.add(gene)
								break
					missing = genes_selected_original - found_genes
					if missing:
						name = record.annotations.get("organism", "Unknown")
						chave = f"{name} ({record_id})"
						lf.write(f"[{chave}]\n")
						for gene in sorted(missing):
							lf.write(f"⚠️ {gene} missing\n")
						lf.write("\n")

			if removal_reasons:
				lf.write(f"[🔎 SEQUENCES REMOVED BY APPLIED FILTERS]: {len(all_records_orig) - len(filtered_records)}\n")
				all_reasons = defaultdict(list)
				for rec_id, reasons in removal_reasons.items():
					for reason in reasons:
						all_reasons[reason].append(rec_id)
						
				for reason, ids in all_reasons.items():
					ids_ordered = sorted(set(ids))
					lf.write(f"{reason}: {', '.join(ids_ordered)}\n")

				if unknown_aliases:
					lf.write("\n")
					lf.write(f"[💡 CUSTOM SYNONYMS TEMPLATE GENERATED]: A custom synonyms template with {len(unknown_aliases)} unrecognized names was generated at 'unrecognized_gene_synonyms.txt'. You can copy, edit, and paste its content directly into 'Add gene synonyms' for your next search.\n\n")

			else:
				lf.write("✅ No sequences were removed by filters.\n")
			lf.write("\n")

		if unknown_aliases:
			template_file = os.path.join(output_dir, "unrecognized_gene_synonyms.txt")
			with open(template_file, "w", encoding="utf-8", newline="\n") as tf:
				tf.write("================================================================================\n")
				tf.write("            PYNNOTATE - UNRECOGNIZED GENE SYNONYMS TEMPLATE                     \n")
				tf.write("================================================================================\n\n")
				tf.write("HOW TO USE THIS TEMPLATE:\n")
				tf.write("1. Copy the formatted YAML/dictionary structure at the bottom of this file.\n")
				tf.write("2. Open Pynnotate GUI and click 'Edit gene synonyms' (or paste into 'Add gene synonyms').\n")
				tf.write("3. Map each unrecognized alias to its correct CANONICAL GENE NAME.\n\n")
				tf.write("CRITICAL INSTRUCTION FOR GROUPING MULTIPLE ALIASES:\n")
				tf.write("- By default, each unrecognized name below has been given its own placeholder key.\n")
				tf.write("- IF MULTIPLE UNRECOGNIZED NAMES REFER TO THE SAME GENE (for example, 'CYTB',\n")
				tf.write("  'Cytb', and 'cytochrome b' all mean CYTB), YOU MUST GROUP THEM UNDER A SINGLE\n\n")
				tf.write("  KEY!\n\n")
				tf.write("EXAMPLE OF GROUPING MULTIPLE SYNONYMS:\n")
				tf.write("  CYTB:\n")
				tf.write("    - CYTB\n")
				tf.write("    - Cytb\n")
				tf.write("    - cytochrome b\n")
				tf.write("  16S:\n")
				tf.write("    - large subunit ribosomal RNA\n")
				tf.write("    - 16S rRNA\n\n")
				tf.write("================================================================================\n")
				tf.write("READY-TO-USE TEMPLATE (EDIT CANONICAL NAMES BEFORE PASTE):\n")
				tf.write("================================================================================\n\n")
			
				template_dict = {}
				for raw_alias in sorted(unknown_aliases.keys()):
					clean_key = raw_alias.upper().replace(" ", "_").replace("-", "_")
					template_dict[clean_key] = [raw_alias]
				
				json.dump(template_dict, tf, indent=4, ensure_ascii=False)
				tf.write("\n")

		metadata_by_gene(gene_dict, excel_data, fields_excel, output_dir, all_records, extract_all_genes)

		if extract_all_genes and len(gene_dict) > 1 and len(all_records) > 1:
			matrix_conflicts = metadata_gene_vs_sample(
				gene_dict, excel_data, header_fields, output_dir, config
			)
			if matrix_conflicts:
				with open(log_file, "a", encoding="utf-8", newline="\n") as lf:
					lf.write("\n[MATRIX CONFLICTS] Multiple GenBank accessions were found for the same matrix taxon and gene:\n")
					for conflict in matrix_conflicts:
						label = conflict["individual"]
						if conflict["organism"]:
							label = f'{conflict["organism"]} | {label}'
						accessions = "; ".join(conflict["accessions"])
						lf.write(f'{label} | {conflict["gene"]}: {accessions}\n')

		from pynnotate.graphic.utils import format_elapsed_time
		elapsed_text = format_elapsed_time(time.monotonic() - download_started_at)
		with open(log_file, "a", encoding="utf-8", newline="\n") as lf:
			lf.write(f"\nElapsed time: {elapsed_text}\n")
		completion_message = f"Download complete!\nElapsed time: {elapsed_text}"

	except Exception as e:
		traceback.print_exc()
		from pynnotate.graphic.utils import format_elapsed_time
		elapsed_text = format_elapsed_time(time.monotonic() - download_started_at)
		completion_message = f"Error processing sequences: {e}\nElapsed time: {elapsed_text}"
		completion_level = "error"
		completion_title = "Error"
	finally:
		if timer_state is not None:
			timer_state["running"] = False
		if root is not None and loading_win is not None:
			try:
				root.after(0, loading_win.destroy)
			except Exception:
				pass
		elif pbar is not None:
			pbar.close()

	if completion_message is not None:
		report_message(
			completion_message, completion_level, completion_title, root
		)
