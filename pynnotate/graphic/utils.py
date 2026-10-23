def center_modal(root, modal, width, height):

    root.update_idletasks()
    root_width = root.winfo_width()
    root_height = root.winfo_height()
    root_x = root.winfo_x()
    root_y = root.winfo_y()

    x = root_x + (root_width // 2) - (width // 2)
    y = root_y + (root_height // 2) - (height // 2)

    modal.geometry(f"{width}x{height}+{x}+{y}")


def copy_text_to_clipboard(root, text):
    """Copy text using Tk's cross-platform clipboard support."""

    root.clipboard_clear()
    root.clipboard_append(text)
    root.update_idletasks()


def format_elapsed_time(seconds):
    """Format an elapsed duration as HH:MM:SS."""

    total_seconds = max(0, int(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def confirm_query_download(root, record_count, query):
    """Show the generated Entrez query and return whether download was approved."""

    import tkinter as tk
    from tkinter import ttk
    from tkinter.scrolledtext import ScrolledText

    modal = tk.Toplevel(root)
    modal.title("Records found")
    modal.minsize(620, 360)
    modal.resizable(True, True)
    modal.transient(root)
    modal.grab_set()
    center_modal(root, modal, width=760, height=460)

    result = {"proceed": False}

    container = ttk.Frame(modal, padding=18)
    container.pack(fill="both", expand=True)

    ttk.Label(
        container,
        text=f"🔢 {record_count} records found",
        font=("Segoe UI", 15, "bold"),
    ).pack(anchor="w", pady=(0, 6))
    ttk.Label(
        container,
        text=(
            "Review the complete Entrez query below. You can copy it and paste it "
            "into GenBank for independent verification."
        ),
        wraplength=710,
        justify="left",
    ).pack(anchor="w", pady=(0, 10))

    query_text = ScrolledText(
        container,
        height=13,
        wrap="word",
        font=("TkFixedFont", 10),
        padx=8,
        pady=8,
    )
    query_text.pack(fill="both", expand=True)
    query_text.insert("1.0", query)
    query_text.configure(state="disabled")

    status_var = tk.StringVar(value="")
    ttk.Label(container, textvariable=status_var, foreground="#2e7d32").pack(
        anchor="w", pady=(5, 0)
    )

    def copy_query():
        copy_text_to_clipboard(root, query)
        status_var.set("Query copied to the clipboard.")

    def close_dialog(proceed=False):
        result["proceed"] = proceed
        try:
            modal.grab_release()
        except Exception:
            pass
        modal.destroy()

    button_frame = ttk.Frame(container)
    button_frame.pack(fill="x", pady=(10, 0))
    ttk.Button(button_frame, text="Copy query", command=copy_query).pack(side="left")
    cancel_button = ttk.Button(
        button_frame,
        text="Cancel",
        command=lambda: close_dialog(False),
    )
    download_button = ttk.Button(
        button_frame,
        text="Start download",
        command=lambda: close_dialog(True),
    )
    download_button.pack(side="right")
    cancel_button.pack(side="right", padx=(8, 8))

    modal.protocol("WM_DELETE_WINDOW", lambda: close_dialog(False))
    modal.bind("<Escape>", lambda _event: close_dialog(False))
    modal.bind("<Return>", lambda _event: close_dialog(True))
    modal.after_idle(download_button.focus_set)
    _activate_modal(root, modal)
    modal.wait_window()
    return result["proceed"]

def show_field_selector(entry_widget, root):

    import tkinter as tk

    valid_fields = ["organism", "isolate", "specimen_voucher", "geo_loc_name",
        "country", "collection_date", "host", "strain", "isolation_source",
        "lat_lon"]

    modal = tk.Toplevel(entry_widget.winfo_toplevel())
    modal.transient(entry_widget.winfo_toplevel())
    modal.title("Select header fields")
    center_modal(root, modal, width=300, height=450)
    _activate_modal(root, modal)
    tk.Label(modal, text="Select header fields:").pack(pady=5)
    tk.Label(modal, text="At least one required field (default: organism)", fg="gray", font=("Arial", 12, "italic")).pack(pady=(0, 10))
    check_vars = {field: tk.BooleanVar(value=(field.lower() == "organism")) for field in valid_fields}
    for field in valid_fields:
        cb = tk.Checkbutton(modal, text=field, variable=check_vars[field])
        cb.pack(anchor="w")

    def apply_selection():
        selected = [f for f, var in check_vars.items() if var.get()]
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, ", ".join(selected))
        modal.destroy()

    tk.Button(modal, text="Apply", command=apply_selection).pack(pady=10)

def show_gene_selector(entry_widget, org_type, root):

    import tkinter as tk

    try:
        from pynnotate.main.synonym_updates import get_official_alias_map
    except ModuleNotFoundError:
        from ..main.synonym_updates import get_official_alias_map

    GENES_MITO_ANIMAL_MAIN = [
        "12S", "16S", "ATP6", "ATP8", "COI", "COII", "COIII",
        "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"
    ]

    GENES_NUCLEAR_ANIMAL_MAIN = [
        "18S", "28S", "ITS1", "ITS2", "RAG1", "RAG2", "BDNF", "CMOS", "POMC", 
        "TYR", "RHO", "SIAH1", "CRYB", "TNS3", "S7", "TYRP1", "EIF3", "CXCR4", 
        "NT3", "FIB7", "MYH6", "CMYC", "CMYC2", "H3A", "ZEB2", "MC1R", "SNCAIP",
        "PTPRC", "OPA1"
    ]

    GENES_MITO_PLANT_MAIN = [
        "26S", "ATP1", "ATP4", "ATP6", "ATP8", "ATP9", "COI", "COII", "COIII",
        "CYTB", "ND1", "ND2", "ND3", "ND4", "NADL", "ND5", "ND6",
        "ND7", "ND9", "MATR", "RPS3", "L16"
    ]

    GENES_CHLOROPLASTO_MAIN = [
        "RBCL", "MATK", "PSBA", "PSBB", "PSBC", "PSBD", "ATPA", "ATPB",
        "NDHA", "NDHB", "NDHC", "NDHD", "NDHE", "NDHF", "NDHG", "NDHI", "NDHK",
        "YCF1", "YCF2", "RPOA", "RPOB", "RPOC1", "RPOC2", "CLPP", "ACCD"
    ]

    if org_type == "animal_mito":
        alias_map = get_official_alias_map(org_type)
        main_list = GENES_MITO_ANIMAL_MAIN
    elif org_type == "animal_nuclear":
        alias_map = get_official_alias_map(org_type)
        main_list = GENES_NUCLEAR_ANIMAL_MAIN
    elif org_type == "plant_mito":
        alias_map = get_official_alias_map(org_type)
        main_list = GENES_MITO_PLANT_MAIN
    elif org_type == "plant_chloro":
        alias_map = get_official_alias_map(org_type)
        main_list = GENES_CHLOROPLASTO_MAIN
    else:
        alias_map = {}
        main_list = []

    genes_available = []
    for g in main_list:
        for k in alias_map:
            if g.upper() == k.upper():
                genes_available.append(k)
                break

    modal = tk.Toplevel(entry_widget.winfo_toplevel())
    modal.transient(entry_widget.winfo_toplevel())
    modal.title("Select genes")
    center_modal(root, modal, width=600, height=480)
    _activate_modal(root, modal)

    label = tk.Label(modal, 
        text="Select the genes you want to search for:",
        font=("Segoe UI", 14, "bold")
    )
    label.pack(padx=0, pady=(30, 5))
    tk.Label(modal, text="💡 Only the main names are displayed.\nSynonyms will be automatically applied in the search.").pack(pady=2)

    frame_genes = tk.Frame(modal)
    frame_genes.pack(padx=10, pady=10)

    vars_dict = {}
    for i, gene in enumerate(sorted(genes_available)):
        var = tk.BooleanVar()
        cb = tk.Checkbutton(frame_genes, text=gene, variable=var)
        cb.grid(row=i % 8, column=i // 8, sticky="w", padx=5, pady=2)
        vars_dict[gene] = var

    def apply_selection():
        selected = [g for g, v in vars_dict.items() if v.get()]
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, ", ".join(selected))
        modal.destroy()

    tk.Button(modal, text="Apply", command=apply_selection).pack(pady=10)

def create_loading_window_with_progress(root):

    import time
    import tkinter as tk
    from tkinter import ttk

    win = tk.Toplevel(root)
    win.title("Downloading")
    root_width = root.winfo_width()
    root_height = root.winfo_height()
    root_x = root.winfo_x()
    root_y = root.winfo_y()

    x = root_x + (root_width // 2) - (440 // 2)
    y = root_y + (root_height // 2) - (155 // 2)
    win.geometry(f"440x155+{x}+{y}")
    win.transient(root)
    win.grab_set()
    tk.Label(win, text="🔄 Downloading and processing GenBank data...\nPlease wait.").pack(pady=(10, 5))
    progress_var = tk.DoubleVar()
    progress = ttk.Progressbar(win, orient="horizontal", mode="determinate", maximum=100, length=330, variable=progress_var)
    progress.pack(pady=(5, 7))

    timer_state = {"started_at": time.monotonic(), "running": True}
    win._pynnotate_timer_state = timer_state
    elapsed_var = tk.StringVar(value="Elapsed time: 00:00:00")
    ttk.Label(win, textvariable=elapsed_var).pack(pady=(0, 8))

    def update_elapsed_time():
        if not timer_state["running"]:
            return
        elapsed = time.monotonic() - timer_state["started_at"]
        elapsed_var.set(f"Elapsed time: {format_elapsed_time(elapsed)}")
        win.after(1000, update_elapsed_time)

    update_elapsed_time()
    win.update_idletasks()
    _activate_modal(root, win)
    return win, progress_var, progress

def show_help(root):
    
    import tkinter as tk
    from tkinter import ttk

    help_window = tk.Toplevel()
    help_window.transient(root)
    help_window.grab_set()
    help_window.title("Gene Synonyms")
    center_modal(root, help_window, width=500, height=400)
    help_window.resizable(False, False)

    emoji_label = tk.Label(help_window, text="🔧", font=("Segoe UI", 28))
    emoji_label.pack(pady=(15, 0))

    text_label = tk.Label(help_window, text="Gene synonyms", font=("Segoe UI", 18, "bold"))
    text_label.pack(pady=(10, 0))

    intro_text = (
        "Pynnotate already includes an internal dictionary of gene name synonyms to aid extraction. "
        "You can provide additional synonyms for genes not automatically recognized. "
        "We recommend running the program first to identify any unrecognized gene synonyms. "
        "Add any missing synonyms here to improve matching. Input example:"
    )
    label_intro = tk.Label(help_window, text=intro_text, font=("Segoe UI", 13), wraplength=460, justify="center")
    label_intro.pack(padx=10, pady=(10, 5))

    example_text = (
        "{\n"
        "    'CYTB': ['CYTOCHROME B', 'CYT B'],\n"
        "    'COI': ['CO1', 'COX1']\n"
        "}"
    )
    label_code = tk.Label(
        help_window, text=example_text, font=("Courier", 13), 
        justify="left", anchor="w", bg="#f4f4f4", padx=10, pady=5
    )
    label_code.pack(padx=20, pady=(0, 5), fill="x")

    attention_text = (
        "⚠️ ATTENTION: When selecting the genome type and adding synonyms, they will be incorporated into the internal dictionary for that specific genome type. "
        "However, if the genome type selected is 'Other', only the synonyms provided by the user will be used."
    )
    label_attention = tk.Label(help_window, text=attention_text, font=("Segoe UI", 13), wraplength=460, justify="center")
    label_attention.pack(padx=10, pady=(5, 10))

    btn_close = ttk.Button(help_window, text="Close", command=help_window.destroy)
    btn_close.pack(pady=(0, 5))

    help_window.after_idle(lambda: _activate_modal(root, help_window))

def open_genome_guide(root):

    import tkinter as tk
    from tkinter import ttk

    modal = tk.Toplevel(root)
    modal.transient(root)
    modal.grab_set()
    modal.title("Gene Synonyms Guide")
    center_modal(root, modal, width=800, height=720)
    modal.resizable(False, False)

    emoji_label = tk.Label(modal, text="🚨", font=("Segoe UI", 28))
    emoji_label.pack(pady=(15, 0))

    text_label = tk.Label(modal, text="How to choose the right genome type", font=("Segoe UI", 18, "bold"))
    text_label.pack(pady=(5, 10))

    guide_widget = tk.Text(
        modal,
        height=12,
        width=55,
        wrap="word",
        font=("Segoe UI", 13),
        bg=modal.cget("bg"),
        bd=0,
        highlightthickness=0,
        relief="flat"
    )
    guide_widget.pack(fill="both", expand=True, padx=20, pady=(0, 10))

    guide_widget.tag_configure("bold_question", font=("Segoe UI", 14, "bold"))
    guide_widget.tag_configure("center", justify="center")

    guide_widget.insert(
        "end",
        "Selecting a genome type loads its synonym dictionary to extract genes. If extraction is disabled, full GenBank entries are downloaded.\n\n",
        "center"
    )

    guide_widget.insert(
        "end",
        "🎯 Which option should you select?\n",
        ("bold_question", "center")
    )

    guide_widget.insert(
        "end",
        "Select 'Animal mitochondria' / 'Animal nuclear' / 'Plant mitochondria' / 'Chloroplast' if all your target genes belong to that specific genomic compartment. "
        "Pynnotate will automatically load its official built-in synonym dictionary for extraction.\n"
        "Select 'Other' if your target genes are NOT in Pynnotate's built-in dictionaries, "
        "or if you are working with organisms/genomes not covered above.\n\n",
        "center"
    )

    guide_widget.insert(
        "end",
        "🌟 What are the impacts of 'genome type' on filtering modes (Step 4)?\n",
        ("bold_question", "center")
    )
    guide_widget.insert(
        "end",
        "'Strict' and 'Flexible' modes rely entirely on your selected dictionary to count and prioritize genes. "
        "When 'Strict' mode selects the individual with the 'highest number of genes', or when 'Flexible' mode combines "
        "complementary genes at the species level, it only counts genes that are successfully recognized by your active "
        "synonym dictionary. If important target genes are missing from the dictionary, they won't be factored into "
        "these filtering decisions, so always make sure your custom synonyms are up to date for your specific dataset!\n\n",
        "center"
    )

    guide_widget.insert(
        "end",
        "✨ Need to download genes from different genome types at once (for example, nuclear and mitochondrial)?\n",
        ("bold_question", "center")
    )
    guide_widget.insert(
        "end",
        "Pynnotate's built-in dictionaries are genome-specific. Selecting 'Animal mitochondria' "
        "will automatically ignore nuclear genes (and vice versa). You can either run two separate "
        "downloads (one for 'Animal mitochondria', one for 'Animal nuclear', for example) or you "
        "can select 'Other' and click 'Edit gene synonyms' below to supply your own combined list of genes.\n\n",
        "center"
    )

    guide_widget.insert(
        "end",
        "🗣️ Genes NOT in Pynnotate's default lists?\n",
        ("bold_question", "center")
    )
    guide_widget.insert(
        "end",
        "Select 'Other' to use ONLY your custom synonyms, or click 'Edit gene synonyms' below "
        "to COMPLEMENT a built-in dictionary.\n\n",
        "center"
    )

    guide_widget.insert(
        "end",
        "❓ Unsure about missing synonyms?\n",
        ("bold_question", "center")
    )
    steps = [
        "1. Run a test search first with your target species or records.",
        "2. Open the 'log.txt' file located inside your output results folder.",
        "3. Look for the section '[UNKNOWN ALIASES]' to see exact names.",
        "4. Open 'unrecognized_gene_synonyms.txt': copy the JSON structure, group equivalent names, and paste it into the box above.",
        "5. Merging rule: Custom synonyms complement built-in dictionaries unless 'Other' is selected (which uses ONLY your custom list)."
    ]
    
    for step in steps:
        guide_widget.insert("end", f"{step}\n", "center")

    guide_widget.configure(state="disabled")

    btn_close = ttk.Button(modal, text="Got it!", command=modal.destroy)
    btn_close.pack(anchor="center", pady=(0, 15))

    modal.after_idle(lambda: _activate_modal(root, modal))

def _activate_modal(root, modal):
    try:
        # Ensure the dialog is truly modal and visible on macOS/Windows/Linux
        parent = root if hasattr(root, 'winfo_exists') else modal.master
        if parent:
            modal.transient(parent)
        try:
            modal.grab_set()        # grab focus
        except Exception:
            pass
        try:
            modal.lift()            # bring to front
        except Exception:
            pass
        modal.after_idle(modal.focus_set)
    except Exception:
        pass
