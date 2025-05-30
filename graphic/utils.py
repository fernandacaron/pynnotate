def show_field_selector(entry_widget):

    import tkinter as tk

    valid_fields = ["organism", "isolate", "specimen_voucher", "country",
        "collection_date", "host", "strain", "isolation_source", "lat_lon"]

    modal = tk.Toplevel(entry_widget.winfo_toplevel())
    modal.transient(entry_widget.winfo_toplevel())
    modal.title("Select header fields")
    modal.geometry("300x400")
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

    tk.Button(modal, text="Aplicar", command=apply_selection).pack(pady=10)

def show_gene_selector(entry_widget, org_type):

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
        alias_map = alias_map_mito_planta
        main_list = GENES_MITO_PLANT_MAIN
    elif org_type == "plant_chloro":
        alias_map = alias_map_cloroplasto
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
    modal.geometry("600x480")

    tk.Label(modal, text="Select the genes you want to search for:").pack(pady=5)
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

def update_checkbox_overlay_state(*args):
    type = var_org_type.get()
    if type in ("plant_mito", "plant_chloro"):
        chk_delete_overlap.config(state="disabled")
        var_delete_overlap.set(False)
    else:
        chk_delete_overlap.config(state="normal")

def update_extraction_options():
    if var_extract_all_genes.get():
        chk_delete_overlap.pack(anchor="center", pady=0)
        chk_log_missing.pack(anchor="center", pady=0)
        update_checkbox_overlay_state()
    else:
        chk_delete_overlap.pack_forget()
        chk_log_missing.pack_forget()

def create_loading_window_with_progress(parent):

    import tkinter as tk
    from tkinter import ttk

    win = tk.Toplevel(parent)
    win.title("Downloading")
    win.geometry("400x120")
    win.transient(parent)
    win.grab_set()
    tk.Label(win, text="🔄 Downloading data from GenBank...\nPlease wait a few minutes.").pack(pady=10)
    progress_var = tk.DoubleVar()
    progress = ttk.Progressbar(win, orient="horizontal", mode="determinate", maximum=100, length=300, variable=progress_var)
    progress.pack(pady=10)
    win.update_idletasks()
    return win, progress_var, progress