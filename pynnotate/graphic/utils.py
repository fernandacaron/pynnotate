def center_modal(root, modal, width, height):

    root.update_idletasks()
    root_width = root.winfo_width()
    root_height = root.winfo_height()
    root_x = root.winfo_x()
    root_y = root.winfo_y()

    x = root_x + (root_width // 2) - (width // 2)
    y = root_y + (root_height // 2) - (height // 2)

    modal.geometry(f"{width}x{height}+{x}+{y}")

def show_field_selector(entry_widget, root):

    import tkinter as tk

    valid_fields = ["organism", "isolate", "specimen_voucher", "country",
        "collection_date", "host", "strain", "isolation_source", "lat_lon"]

    modal = tk.Toplevel(entry_widget.winfo_toplevel())
    modal.transient(entry_widget.winfo_toplevel())
    modal.title("Select header fields")
    center_modal(root, modal, width=300, height=450)
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

    from main.alias_maps import alias_map_animal, alias_map_mito_plant, alias_map_chloroplast

    GENES_MITO_ANIMAL_MAIN = [
        "12S", "16S", "ATP6", "ATP8", "COI", "COII", "COIII",
        "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"
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
        alias_map = alias_map_animal
        main_list = GENES_MITO_ANIMAL_MAIN
    elif org_type == "plant_mito":
        alias_map = alias_map_mito_plant
        main_list = GENES_MITO_PLANT_MAIN
    elif org_type == "plant_chloro":
        alias_map = alias_map_chloroplast
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

    label = tk.Label(modal, 
        text="Select the genes you want to search for:",
        font=("Segoe UI", 14, "bold")
    )
    label.pack(padx=0, pady=(30, 5))
    tk.Label(modal, text="💡 Only the main names are displayed.\nSynonyms will be automatically applied in the search..").pack(pady=2)

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

    import tkinter as tk
    from tkinter import ttk

    win = tk.Toplevel(root)
    win.title("Downloading")
    root_width = root.winfo_width()
    root_height = root.winfo_height()
    root_x = root.winfo_x()
    root_y = root.winfo_y()

    x = root_x + (root_width // 2) - (400 // 2)
    y = root_y + (root_height // 2) - (120 // 2)
    win.geometry(f"400x120+{x}+{y}")
    win.transient(root)
    win.grab_set()
    tk.Label(win, text="🔄 Downloading data from GenBank...\nPlease wait a few minutes.").pack(pady=10)
    progress_var = tk.DoubleVar()
    progress = ttk.Progressbar(win, orient="horizontal", mode="determinate", maximum=100, length=300, variable=progress_var)
    progress.pack(pady=10)
    win.update_idletasks()
    return win, progress_var, progress

def check_errors_gui_loop(root, error_queue):

    import queue
    from tkinter import messagebox

    try:
        while True:
            error = error_queue.get_nowait()
            messagebox.showerror("Error", f"Error processing sequences: {error}")
    except queue.Empty:
        pass
    root.after(500, lambda: check_errors_gui_loop(root, error_queue))

def open_synonyms_editor(root):

    import tkinter as tk
    from tkinter import ttk

    modal = tk.Toplevel(root)
    modal.title("Edit Gene Synonyms")
    root_width = root.winfo_width()
    root_height = root.winfo_height()
    root_x = root.winfo_x()
    root_y = root.winfo_y()
    x = root_x + (root_width // 2) - (700 // 2)
    y = root_y + (root_height // 2) - (450 // 2)
    modal.geometry(f"700x450+{x}+{y}")
    modal.transient(root)
    modal.grab_set()
    modal.resizable(False, False)

    result = None

    container = ttk.Frame(modal, borderwidth=50)
    container.pack(fill="both", expand=True, padx=10, pady=10)

    label = ttk.Label(
        container,
        text="🔧 Provide additional gene name synonyms in JSON format:",
        font=("Segoe UI", 14, "bold")
    )
    label.pack(pady=(10, 10))

    text_synonyms = tk.Text(container, height=10, wrap="word", font=("Segoe UI", 10))
    text_synonyms.pack(fill="both", expand=True, pady=(0, 10))

    example = '''{
        "CYTB": ["CYTOCHROME B", "CYT B"],
        "COI": ["CO1", "COX1"]
    }'''

    text_synonyms.tag_configure("placeholder", foreground="gray")
    text_synonyms.insert("1.0", example, "placeholder")

    def clear_placeholder(event):
        if text_synonyms.get("1.0", "end-1c") == example:
            text_synonyms.delete("1.0", tk.END)
            text_synonyms.tag_remove("placeholder", "1.0", "end")

    def add_placeholder(event):
        if text_synonyms.get("1.0", "end-1c").strip() == "":
            text_synonyms.insert("1.0", example, "placeholder")

    def save_and_close():
        nonlocal result
        content = text_synonyms.get("1.0", "end-1c").strip()
        if content != "" and content != example:
            result = content
        modal.destroy()

    text_synonyms.bind("<FocusIn>", clear_placeholder)
    text_synonyms.bind("<FocusOut>", add_placeholder)

    button_frame = ttk.Frame(container)
    button_frame.pack(pady=5)

    btn_cancel = ttk.Button(button_frame, text="Cancel", command=modal.destroy)
    btn_cancel.pack(side="left", padx=5)

    btn_save = ttk.Button(button_frame, text="Save", command=save_and_close)
    btn_save.pack(side="left", padx=5)

    root.wait_window(modal)

    return result

def show_help(root):
    
    import tkinter as tk
    from tkinter import ttk

    help_window = tk.Toplevel()
    help_window.transient(root)
    help_window.grab_set()
    help_window.title("Gene Synonyms")
    center_modal(root, help_window, width=500, height=350)
    help_window.resizable(False, False)

    emoji_label = tk.Label(help_window, text="🔧", font=("Segoe UI", 28))
    emoji_label.pack(pady=(15, 0))

    text_label = tk.Label(help_window, text="Gene synonyms", font=("Segoe UI", 18, "bold"))
    text_label.pack(pady=(10, 0))

    help_text = (
        "Pynnotate already includes an internal dictionary of gene name synonyms to aid extraction. "
        "You can provide additional synonyms for genes not automatically recognized. "
        "We recommend running the program first to identify any unrecognized gene synonyms. "
        "Add any missing synonyms here to improve matching.\n\n"
        "⚠️ ATTENTION: When selecting the genome type and adding synonyms, they will be incorporated into the internal dictionary for that specific genome type. "
        "However, if the genome type selected is 'Other', only the synonyms provided by the user will be used."
    )
    label = tk.Label(help_window, text=help_text, font=("Segoe UI", 13), wraplength=460, justify="center")
    label.pack(padx=10, pady=10)

    btn_close = ttk.Button(help_window, text="Close", command=help_window.destroy)
    btn_close.pack(pady=(0, 5))
