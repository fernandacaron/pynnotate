# -*- coding: utf-8 -*-
# pynnotate_graphic.py (revisado)

import sys, os, traceback, pathlib, logging, queue, threading, webbrowser
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from PIL import Image, ImageTk

# --- Patch SSL para certificados (tolerante) ---
import ssl
try:
    import certifi
    ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
except Exception:
    # Se certifi não estiver disponível, segue com o contexto padrão
    pass

# ===== Helpers para rodar empacotado e logar erros =====
def _app_base():
    # Quando empacotado, PyInstaller expõe _MEIPASS com os arquivos do bundle
    return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))

BASE = _app_base()

def resource_path(*parts):
    """Resolve caminho tanto em dev quanto no bundle do PyInstaller."""
    return os.path.join(BASE, *parts)

# (Opcional, mas ajuda com código que ainda usa caminhos relativos)
if getattr(sys, "_MEIPASS", None):
    os.chdir(BASE)

# Diretório de logs do app
LOG_DIR = pathlib.Path.home() / ".pynnotate" / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / "app_startup.log"
logging.basicConfig(
    filename=str(LOG_FILE),
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# ===== Imports do projeto =====
from pynnotate import __version__
from pynnotate.graphic.utils import (
    show_gene_selector, show_field_selector, show_help, open_genome_guide
)
from pynnotate.main.main import begin_search
from pynnotate.main.synonym_updates import (
    check_for_synonym_updates,
    get_synonym_database_info,
    install_synonym_update,
    restore_bundled_synonym_database,
)


# ========= Lógica do botão principal (começar busca) =========
def start_download():
    """Coleta configurações da GUI e chama begin_search()."""
    selected_mode = filter_mode_var.get()

    if selected_mode == "strict":
        prioritize_more_genes = True
        unique_species = False
    elif selected_mode == "flexible":
        prioritize_more_genes = False
        unique_species = True
    else:
        prioritize_more_genes = False
        unique_species = False

    config = {
        "email": entry_email.get().strip(),
        "title": entry_title.get().strip(),
        "organisms": entry_org.get().strip(),
        "genes": entry_genes.get().strip(),
        "fields": entry_fields.get().strip(),
        "output_name": entry_path.get().strip(),
        "folder_name": folder_name_var.get().strip(),
        "include_id": var_id.get(),
        "mito": var_mito.get(),
        "mitogenome": var_mitogenome.get(),
        "ids_text": entry_ids.get("1.0", tk.END).strip(),
        "chloroplast": var_chloroplast.get(),
        "min_len": entry_min_len.get(),
        "max_len": entry_max_len.get(),
        "delete_unverified": var_delete_unverified.get(),
        "extraction": var_extract_all_genes.get(),
        "overlap": var_delete_overlap.get(),
        "unique": unique_species,
        "prioritize_more_genes": prioritize_more_genes,
        "logmissing": var_log_missing.get(),
        "org_type": var_org_type.get(),
        "api_key": entry_api_key.get().strip(),
        "add_synonyms": entry_syn.get("1.0", tk.END).strip(),
        "additional": entry_add.get().strip(),
    }
    begin_search(config, root=root)


# ========= Construção da GUI =========
def create_centered_frame():
    frame = ttk.Frame(root, padding=20)
    frame.grid(row=0, column=0, sticky="nsew")

    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)

    content = ttk.Frame(frame)
    content.grid(row=0, column=0)

    frames.append((frame, content))
    return frame, content


def set_synonyms_from_editor(root, entry_widget):
    current_value = entry_widget.get("1.0", tk.END).strip()

    def apply_result(result):
        entry_widget.delete("1.0", tk.END)
        if result:
            entry_widget.insert("1.0", result)

    open_synonyms_editor(
        root,
        on_save=apply_result,
        initial_value=current_value,
    )


def _synonym_status_text():
    info = get_synonym_database_info()
    source = "downloaded update" if info["source"] == "downloaded" else "included with app"
    return (
        f"Current base: {info['dictionary_version']} "
        f"(revision {info['database_revision']}, {source})"
    )


def refresh_synonym_status():
    try:
        synonym_status_var.set(_synonym_status_text())
    except Exception as exc:
        synonym_status_var.set(f"Synonym database unavailable: {exc}")


def check_synonym_updates_from_gui():
    """Check manually without blocking Tk's main thread."""
    btn_check_synonyms.config(state="disabled")
    btn_restore_synonyms.config(state="disabled")
    synonym_status_var.set("Checking the official synonym database...")
    result_queue = queue.Queue()

    def worker():
        try:
            result_queue.put(("ok", check_for_synonym_updates()))
        except Exception as exc:
            result_queue.put(("error", exc))

    def poll_result():
        try:
            status, payload = result_queue.get_nowait()
        except queue.Empty:
            root.after(100, poll_result)
            return

        btn_check_synonyms.config(state="normal")
        btn_restore_synonyms.config(state="normal")
        if status == "error":
            refresh_synonym_status()
            messagebox.showerror(
                "Synonym update",
                "The update could not be checked. The current local database was kept.\n\n"
                f"Details: {payload}",
            )
            return

        if not payload["update_available"]:
            refresh_synonym_status()
            messagebox.showinfo(
                "Synonym update",
                "The synonym database is already up to date!\n\n"
                f"Local: {payload['current_version']} (revision {payload['current_revision']})\n"
                f"Official: {payload['remote_version']} (revision {payload['remote_revision']})",
            )
            return

        changes = payload["changes"]
        install = messagebox.askyesno(
            "Synonym update available",
            f"A new complete database is available.\n\n"
            f"Current: {payload['current_version']} (revision {payload['current_revision']})\n"
            f"New: {payload['remote_version']} (revision {payload['remote_revision']})\n\n"
            f"Genes added/removed: {changes['genes_added']}/{changes['genes_removed']}\n"
            f"Aliases added/removed: {changes['aliases_added']}/{changes['aliases_removed']}\n\n"
            "Install it now?",
        )
        if not install:
            refresh_synonym_status()
            return
        try:
            info = install_synonym_update(payload)
            refresh_synonym_status()
            messagebox.showinfo(
                "Synonym update",
                f"Database {info['dictionary_version']} installed successfully.\n"
                "All four official genome dictionaries were updated.",
            )
        except Exception as exc:
            refresh_synonym_status()
            messagebox.showerror(
                "Synonym update",
                "The update was not installed; the previous database was preserved.\n\n"
                f"Details: {exc}",
            )

    threading.Thread(target=worker, name="synonym-update-check", daemon=True).start()
    root.after(100, poll_result)


def restore_bundled_synonyms_from_gui():
    try:
        info = get_synonym_database_info()
    except Exception as exc:
        messagebox.showerror("Synonym database", f"The database could not be read.\n\nDetails: {exc}")
        return
    if info["source"] == "bundled":
        messagebox.showinfo("Synonym database", "The app is already using its included database.")
        return
    if not messagebox.askyesno(
        "Restore included database",
        "Remove the downloaded update and return to the database included with this app?",
    ):
        return
    try:
        restore_bundled_synonym_database()
        refresh_synonym_status()
        messagebox.showinfo("Synonym database", "The included database was restored.")
    except Exception as exc:
        refresh_synonym_status()
        messagebox.showerror(
            "Synonym database",
            f"The downloaded database could not be removed.\n\nDetails: {exc}",
        )


def suggest_synonym_on_github():
    url = (
        "https://github.com/fernandacaron/pynnotate/issues/new"
        "?title=Synonym%20database%20suggestion"
        "&body=Genome%20type%3A%0ACanonical%20gene%3A%0ANew%20alias%3A%0AReference%20or%20example%20accession%3A"
    )
    webbrowser.open(url)


def update_checkbox_overlay_state(*args):
    org_type = var_org_type.get()
    if org_type in ("plant_mito", "plant_chloro"):
        chk_delete_overlap.config(state="disabled")
        var_delete_overlap.set(False)
    else:
        chk_delete_overlap.config(state="normal")


def update_extraction_warning(*args):
    if var_extract_all_genes.get():
        warning_extract_all_genes.pack(anchor="center", padx=5, pady=(10, 5))
    else:
        warning_extract_all_genes.pack_forget()


def update_extraction_options(*args):
    if filter_mode_var.get() == "flexible" and not var_extract_all_genes.get():
        var_extract_all_genes.set(True)
        return
    if var_extract_all_genes.get():
        chk_delete_overlap.pack(anchor="center", padx=5, pady=(0, 8))
        chk_delete_overlap.config(state="normal")
        chk_log_missing.pack(anchor="center", padx=5, pady=(0, 8))
        chk_log_missing.config(state="normal")
        update_checkbox_overlay_state()
    else:
        chk_delete_overlap.config(state="disabled")
        var_delete_overlap.set(False)
        chk_log_missing.config(state="disabled")
        var_log_missing.set(False)


def update_filtering_warning(*args):
    if filter_mode_var.get() == "strict":
        warning_label_strict.pack(anchor="center", padx=5, pady=(10, 5))
        warning_label_unconstrained.pack_forget()
        warning_label_flexible.pack_forget()
    elif filter_mode_var.get() == "unconstrained":
        warning_label_strict.pack_forget()
        warning_label_unconstrained.pack(anchor="center", padx=5, pady=(10, 5))
        warning_label_flexible.pack_forget()
    else:
        var_extract_all_genes.set(True)
        warning_label_strict.pack_forget()
        warning_label_unconstrained.pack_forget()
        warning_label_flexible.pack(anchor="center", padx=5, pady=(10, 5))


def choose_folder():
    folder = filedialog.askdirectory(title="Choose a folder to save the results")
    if folder:
        path_var.set(folder)


def show_frame(index):
    global current_index
    for f, _ in frames:
        f.grid_remove()
    frames[index][0].grid()
    current_index = index


def next_frame():
    if current_index + 1 < len(frames):
        show_frame(current_index + 1)


def previous_frame():
    if current_index - 1 >= 0:
        show_frame(current_index - 1)

def show_about():
    from tkinter import messagebox
    messagebox.showinfo(
        "About Pynnotate",
        f"Pynnotate v{__version__}\n"
        "© 2025 Fernanda Caron, Felipe Magalhães, Matheus Salles & Fabricius Domingos\n\n"
        "For more information, please check: https://github.com/fernandacaron/pynnotate\n"
    )


# ========= Inicialização segura (com logging) =========
try:
    # Janela principal
    root = tk.Tk()
    root.title(f"Pynnotate v{__version__}")

    # Ícone / logo — busca logo_transparente.png em locais prováveis e mantém referência viva
    try:
        # Preferimos PNG via iconphoto (macOS/Linux). Windows também aceita.
        icon_candidates = [
            resource_path("data", "pynnotate.png"),
            resource_path("logo_transparente.png"),
        ]
        icon_path = next((p for p in icon_candidates if os.path.exists(p)), None)
        if icon_path:
            _logo_img_small = Image.open(icon_path).resize((64, 64))
            _logo_photo_small = ImageTk.PhotoImage(_logo_img_small)
            root.iconphoto(True, _logo_photo_small)
            root._icon_ref = _logo_photo_small  # evita GC

    except Exception as e:
        logging.warning("Falha ao definir ícone/logo da janela: %s", e)

    # Posição/tamanho: tenta screeninfo, senão fallback
    try:
        from screeninfo import get_monitors
        monitors = get_monitors()
        if len(monitors) > 1:
            m = monitors[1]
            x = m.x + 520
            y = m.y + 100
        else:
            m = monitors[0]
            x = m.x + m.width // 2 - 300
            y = m.y + m.height // 2 - 200
        root.geometry(f"650x700+{x}+{y}")
    except Exception as e:
        logging.warning("screeninfo indisponível (%s). Usando fallback.", e)
        root.geometry("650x700+100+80")

    if sys.platform == "darwin":
        _base_scaling = root.tk.call('tk', 'scaling')
        _last_monitor_geom = {"value": None}

        def _get_current_monitor_geom():
            try:
                from screeninfo import get_monitors
                wx = root.winfo_x()
                wy = root.winfo_y()
                for m in get_monitors():
                    if m.x <= wx < m.x + m.width and m.y <= wy < m.y + m.height:
                        return (m.x, m.y, m.width, m.height)
            except Exception as e:
                logging.warning("screeninfo indisponível no monitor de reajuste: %s", e)
            return None

        def _reapply_scaling_and_geometry(event=None):
            geom = _get_current_monitor_geom()
            if geom is None or geom == _last_monitor_geom["value"]:
                return
            _last_monitor_geom["value"] = geom
            try:
                root.tk.call('tk', 'scaling', _base_scaling)
                w = root.winfo_width()
                h = root.winfo_height()
                x = root.winfo_x()
                y = root.winfo_y()
                if w > 1 and h > 1:
                    root.geometry(f"{w}x{h}+{x}+{y}")
                root.update_idletasks()
            except Exception as e:
                logging.warning("Falha ao reaplicar escala/geometria: %s", e)
        _configure_after_id = {"id": None}

        def _on_configure(event):
            if event.widget is not root:
                return
            if _configure_after_id["id"] is not None:
                root.after_cancel(_configure_after_id["id"])
            _configure_after_id["id"] = root.after(200, _reapply_scaling_and_geometry)

        root.bind("<Configure>", _on_configure)

    # Tema
    style = ttk.Style()
    style.configure("TLabel", font=("Segoe UI", 14))
    style.configure("TButton", font=("Segoe UI", 12), padding=6)
    style.configure("TRadiobutton", font=("Segoe UI", 12))
    style.configure("Bold.TCheckbutton", font=("Arial", 14, "bold"))
    style.configure('type.TRadiobutton', font=("Segoe UI", 13))

    # Frames
    current_index = 0
    frames = []

    # ---- Frame 0 ----
    frame0, content0 = create_centered_frame()

    try:
        # tenta nas pastas mais prováveis: data/, pynnotate/, pynnotate/main/
        candidates = [
            resource_path("data", "pynnotate.png"),
            resource_path("logo_transparente.png"),
        ]
        logo_path = next((p for p in candidates if os.path.exists(p)), None)
        if not logo_path:
            raise FileNotFoundError("Pynnotate logo not found in the application resources")

        logo_big = Image.open(logo_path).resize((200, 200))
        logo_photo = ImageTk.PhotoImage(logo_big)
        logo_label = tk.Label(content0, image=logo_photo)
        logo_label.image = logo_photo  # evita GC
        logo_label.pack(pady=(0, 10))

        # opcional: também define o ícone da janela com o mesmo PNG
        try:
            root.iconphoto(True, logo_photo)
            root._icon_ref = logo_photo  # evita GC
        except Exception:
            pass

    except Exception as e:
        logging.warning("Não foi possível carregar o logo grande: %s", e)

    msg_label = ttk.Label(content0, text="Welcome to Pynnotate!", font=("Segoe UI", 18, "bold"))
    msg_label.pack(pady=(0, 20))

    intro_text = (
        "Pynnotate is your go-to tool for efficient gene extraction and analysis.\n\n"
        "Explore tutorials and tips to get the most out of the app."
    )
    intro_label = ttk.Label(content0, text=intro_text, font=("Segoe UI", 12), justify="center")
    intro_label.pack(pady=(0, 20))

    next_btn = ttk.Button(content0, text="Start!", command=next_frame)
    next_btn.pack(pady=10)

    btn_about0 = ttk.Button(content0, text="Help / About", command=show_about)
    btn_about0.pack(pady=(0, 10))


    # ---- Frame 1 ----
    frame1, content1 = create_centered_frame()

    label_title = ttk.Label(content1, text="Step 1: Choose genome type", font=("Segoe UI", 18, "bold"))
    label_title.pack(pady=(0, 10))

    info_box = ttk.Frame(content1, padding=2)
    info_box.pack(fill="x", padx=10, pady=(4, 4))

    btn_guide_popup = tk.Button(
        info_box, text="⚠️ READ THIS BEFORE CONTINUING ⚠️", 
        font=("Segoe UI", 12, "bold"),
        command=lambda: open_genome_guide(root), padx=12, pady=5
    )
    btn_guide_popup.pack(anchor="center")

    var_org_type = tk.StringVar(value="animal_mito")
    var_org_type.trace_add("write", update_checkbox_overlay_state)

    options_frame = ttk.Frame(content1)
    options_frame.pack(pady=10)

    r1 = ttk.Radiobutton(options_frame, text="Animal mitochondria", variable=var_org_type, value="animal_mito", style="type.TRadiobutton")
    r1.pack(anchor="w", pady=6)
    r_nuc = ttk.Radiobutton(options_frame, text="Animal nuclear", variable=var_org_type, value="animal_nuclear", style="type.TRadiobutton")
    r_nuc.pack(anchor="w", pady=6)
    r2 = ttk.Radiobutton(options_frame, text="Plant mitochondria", variable=var_org_type, value="plant_mito", style="type.TRadiobutton")
    r2.pack(anchor="w", pady=6)
    r3 = ttk.Radiobutton(options_frame, text="Chloroplast", variable=var_org_type, value="plant_chloro", style="type.TRadiobutton")
    r3.pack(anchor="w", pady=6)
    r4 = ttk.Radiobutton(options_frame, text="Other", variable=var_org_type, value="other", style="type.TRadiobutton")
    r4.pack(anchor="w", pady=6)

    frame_syn_header = ttk.Frame(content1)
    frame_syn_header.pack(fill="x", pady=(15, 0))

    inner_header_frame = ttk.Frame(frame_syn_header)
    inner_header_frame.pack(anchor="center")

    info_button = tk.Button(
        inner_header_frame, text="?", font=("Segoe UI", 12, "bold"),
        fg="#555", bg='#f0f0f0', bd=0, relief="flat", highlightthickness=0,
        command=lambda: show_help(root), width=1
    )
    info_button.pack(side="left", padx=(0, 5))

    label_syn = ttk.Label(
        inner_header_frame, 
        text="🔧 Add gene synonyms [optional] (JSON format):", 
        font=("Segoe UI", 14, "bold")
    )
    label_syn.pack(side="left")
    
    entry_syn = tk.Text(content1, height=5, wrap="word", font=("Segoe UI", 11))
    entry_syn.pack(fill="x", padx=25, pady=(5, 10))

    style = ttk.Style()
    style.configure("Custom.TLabelframe.Label", font=("Segoe UI", 14, "bold"))

    synonym_update_frame = ttk.LabelFrame(
        content1,
        text="Official synonym database",
        padding=10,
        style="Custom.TLabelframe"
    )
    synonym_update_frame.pack(fill="x", padx=25, pady=(10, 5))

    synonym_status_var = tk.StringVar()
    synonym_status_label = ttk.Label(
        synonym_update_frame,
        textvariable=synonym_status_var,
        font=("Segoe UI", 12),
    )
    synonym_status_label.pack(anchor="center", pady=(0, 7))

    synonym_update_buttons = ttk.Frame(synonym_update_frame)
    synonym_update_buttons.pack(anchor="center")
    btn_check_synonyms = ttk.Button(
        synonym_update_buttons,
        text="Check for updates",
        command=check_synonym_updates_from_gui,
    )
    btn_check_synonyms.pack(side="left", padx=4)
    btn_restore_synonyms = ttk.Button(
        synonym_update_buttons,
        text="Restore included",
        command=restore_bundled_synonyms_from_gui,
    )
    btn_restore_synonyms.pack(side="left", padx=4)
    btn_suggest_synonym = ttk.Button(
        synonym_update_buttons,
        text="Suggest a synonym",
        command=suggest_synonym_on_github,
    )
    btn_suggest_synonym.pack(side="left", padx=4)
    ttk.Label(
        synonym_update_frame,
        text="Manual update of the complete 'Animal mitochondrial', 'Plant mitochondrial', and 'Chloroplast' base.",
        font=("Segoe UI", 12),
        foreground="#555555",
        wraplength=540,
        justify="center",
    ).pack(anchor="center", pady=(7, 0))
    refresh_synonym_status()

    button_frame = ttk.Frame(content1)
    btn_previous = ttk.Button(button_frame, text="◀ Previous", command=previous_frame)
    btn_previous.pack(side="left", padx=10)
    btn_next = ttk.Button(button_frame, text="Next ▶", command=next_frame)
    btn_next.pack(side="left", padx=10)
    button_frame.pack(pady=20)

    # ---- Frame 2 ----
    frame2, content2 = create_centered_frame()

    label_title = ttk.Label(content2, text="Step 2: Build custom search query", font=("Segoe UI", 18, "bold"))
    label_title.pack(pady=(0, 5))

    ttk.Label(content2, text="📋 GenBank IDs (one per line or comma separated) [optional]:").pack(anchor="center", pady=(10, 5))
    entry_ids = tk.Text(content2, height=3, width=100)
    entry_ids.pack(pady=(0, 10), anchor="w")

    ttk.Label(content2, text="🧬 Gene(s) [ex: 16S, COI] [optional]:").pack(anchor="center", pady=(10, 5))
    frame_genes = ttk.Frame(content2)
    frame_genes.pack(fill="x", pady=(0, 10), padx=0)
    entry_genes = ttk.Entry(frame_genes, width=50)
    entry_genes.grid(row=0, column=0, sticky="ew", padx=(0, 5))
    btn_genes = ttk.Button(frame_genes, text="Select genes", command=lambda: show_gene_selector(entry_genes, var_org_type.get(), root))
    btn_genes.grid(row=0, column=1)
    frame_genes.columnconfigure(0, weight=1)

    ttk.Label(content2, text="🧫 Organism(s) [optional]:").pack(anchor="center", pady=(10, 5))
    entry_org = ttk.Entry(content2, width=100)
    entry_org.pack(pady=(0, 10), anchor="w")

    ttk.Label(content2, text="🔎 Publication term (title, authors, year) [optional]:").pack(anchor="center", pady=(10, 5))
    entry_title = ttk.Entry(content2, width=100)
    entry_title.pack(pady=(0, 10), anchor="w")

    ttk.Label(content2, text="💭 Any additional query [ex: NOT sp.] [optional]:").pack(anchor="center", pady=(10, 5))
    entry_add = ttk.Entry(content2, width=100)
    entry_add.pack(pady=(0, 10), anchor="w")

    ttk.Label(content2, text="Refine your search terms to [optional]:").pack(anchor="center", pady=(10, 5))

    frame_checkboxes = ttk.Frame(content2)
    var_mito = tk.BooleanVar(value=False)
    chk_mito = ttk.Checkbutton(frame_checkboxes, text="Mitochondrial gene", variable=var_mito)
    chk_mito.pack(side="left", padx=5)
    var_mitogenome = tk.BooleanVar(value=False)
    chk_mitogenome = ttk.Checkbutton(frame_checkboxes, text="Mitogenome", variable=var_mitogenome)
    chk_mitogenome.pack(side="left", padx=5)
    var_chloroplast = tk.BooleanVar(value=False)
    chk_chloroplast = ttk.Checkbutton(frame_checkboxes, text="Chloroplast", variable=var_chloroplast)
    chk_chloroplast.pack(side="left", padx=5)
    var_delete_unverified = tk.BooleanVar(value=False)
    chk_delete_unverified = ttk.Checkbutton(frame_checkboxes, text="Delete unannotated", variable=var_delete_unverified)
    chk_delete_unverified.pack(side="left", padx=5)
    frame_checkboxes.pack(pady=(0, 10), anchor="center")

    button_frame2 = ttk.Frame(content2)
    btn_previous = ttk.Button(button_frame2, text="◀ Previous", command=previous_frame)
    btn_previous.pack(side="left", padx=10)
    btn_next_2 = ttk.Button(button_frame2, text="Next ▶", command=next_frame)
    btn_next_2.pack(side="left", padx=10)
    button_frame2.pack(pady=20)

    # ---- Frame 3 ----
    frame3, content3 = create_centered_frame()

    label_title = ttk.Label(content3, text="Step 3: Formatting and filtering options for download", font=("Segoe UI", 18, "bold"))
    label_title.pack(pady=(0, 20))

    label_flexible = ttk.Label(content3, text="🧪 Choose the filtering mode:", font=("Segoe UI", 15, "bold"))
    label_flexible.pack(pady=(25, 5))
    
    filter_mode_var = tk.StringVar(value="unconstrained")

    rbtn_unconstrained = ttk.Radiobutton(content3, text="🌐 Unconstrained: Retain all retrieved records", variable=filter_mode_var, value="unconstrained")
    rbtn_unconstrained.pack(anchor="center", pady=(5, 2))

    rbtn_strict = ttk.Radiobutton(
        content3,
        text="🔒 Strict: Select one representative individual per species",
        variable=filter_mode_var, value="strict"
    )
    rbtn_strict.pack(anchor="center", pady=2)

    rbtn_flexible = ttk.Radiobutton(
        content3,
        text="🌱 Flexible: Combine complementary genes at the species level",
        variable=filter_mode_var, value="flexible"
    )
    rbtn_flexible.pack(anchor="center", pady=(2, 10))

    warning_label_strict = ttk.Label(
        content3,
        text=("ℹ️ Strict Mode: This mode filters and selects only one representative individual per species, "
              "prioritizing specimens with the highest number of genes. Note that filtering is strictly based "
              "on the gene set defined in your selected genome type and/or custom user-provided synonym dictionary."),
        wraplength=550, justify="center", font=("Segoe UI", 11)
    )
    warning_label_unconstrained = ttk.Label(
        content3,
        text=("ℹ️ Unconstrained Mode: Retains all retrieved records without filtering out duplicates or multiple sequences per species. "
              "If separate gene extraction is enabled in Step 4, only genes recognised by the selected and/or "
              "user-provided synonym dictionary will be extracted."),
        wraplength=550, justify="center", font=("Segoe UI", 11)
    )
    warning_label_flexible = ttk.Label(
        content3,
        text=("ℹ️ Flexible Mode: Automatically builds species supermatrices by combining complementary genes from different "
              "identified individuals. Gene FASTA headers contain only the species name; selected header fields and the "
              "GenBank-ID option apply exclusively to sequences.fasta, while full accession, isolate, and voucher "
              "provenance details are securely retained in genes_matrix.xlsx."),
        wraplength=550, justify="center", font=("Segoe UI", 11)
    )

    filter_mode_var.trace_add("write", update_filtering_warning)

    ttk.Label(content3, text="📋 Complete sequences.fasta header fields (select):").pack(anchor="center", pady=(30, 5))
    fields_frame = ttk.Frame(content3)
    fields_frame.pack(fill="x", pady=(0, 10), padx=0)
    entry_fields = ttk.Entry(fields_frame, width=50)
    entry_fields.grid(row=0, column=0, sticky="ew", padx=(0, 5))
    entry_fields.insert(0, "organism")
    btn_fields = ttk.Button(fields_frame, text="Select", command=lambda: show_field_selector(entry_fields, root))
    btn_fields.grid(row=0, column=1)
    fields_frame.columnconfigure(0, weight=1)

    center_frame = ttk.Frame(content3)
    center_frame.pack(pady=2)

    var_id = tk.BooleanVar(value=True)
    chk_id = ttk.Checkbutton(center_frame, text="Include GenBank ID in fasta header [optional]", variable=var_id)
    chk_id.pack(anchor="center", pady=10)

    frame_min_len = ttk.Frame(content3)
    frame_min_len.pack(anchor="w", pady=(10, 10), fill="x")
    ttk.Label(frame_min_len, text="🔢 Minimum sequence size (bp) [optional]:").grid(row=0, column=0, sticky="w", padx=(0, 10))
    entry_min_len = ttk.Entry(frame_min_len, width=20)
    entry_min_len.grid(row=0, column=1, sticky="ew")
    frame_min_len.columnconfigure(1, weight=1)

    frame_max_len = ttk.Frame(content3)
    frame_max_len.pack(anchor="w", pady=(0, 10), fill="x")
    ttk.Label(frame_max_len, text="🔢 Maximum sequence size (bp) [optional]:").grid(row=0, column=0, sticky="w", padx=(0, 10))
    entry_max_len = ttk.Entry(frame_max_len, width=20)
    entry_max_len.grid(row=0, column=1, sticky="ew")
    frame_max_len.columnconfigure(1, weight=1)

    button_frame3 = ttk.Frame(content3)
    button_frame3.pack(pady=20)
    btn_previous3 = ttk.Button(button_frame3, text="◀ Previous", command=previous_frame)
    btn_previous3.pack(side="left", padx=10)
    btn_next3 = ttk.Button(button_frame3, text="Next ▶", command=next_frame)
    btn_next3.pack(side="left", padx=10)

    warning_label_unconstrained.pack(anchor="center", padx=5, pady=(10, 5))

    # ---- Frame 4 ----
    frame4, content4 = create_centered_frame()

    label_title = ttk.Label(content4, text="Step 4: Additional options", font=("Segoe UI", 18, "bold"))
    label_title.pack(pady=(0, 20))

    var_extract_all_genes = tk.BooleanVar(value=True)
    var_delete_overlap = tk.BooleanVar(value=False)
    var_log_missing = tk.BooleanVar(value=False)

    chk_extract_all_genes = ttk.Checkbutton(
        content4,
        text="🌟 Extract all annotated genes separately [optional]",
        variable=var_extract_all_genes,
        command=lambda: update_extraction_options(),
        style="Bold.TCheckbutton"
    )
    chk_extract_all_genes.pack(anchor="center", padx=5, pady=(0, 8))

    warning_extract_all_genes = ttk.Label(
        content4,
        text=("⚠️ WARNING: Gene extraction will be limited to the selected synonym dictionary. "
              "For example, when selecting 'plant_chloro', only chloroplast genes will be extracted."
              "This is automatically required in 'Flexible' Mode, but optional in 'Strict' and 'Unconstrained' modes."),
        wraplength=550, justify="center", font=("Segoe UI", 12)
    )
    var_extract_all_genes.trace_add("write", update_extraction_warning)

    chk_delete_overlap = ttk.Checkbutton(
        content4,
        text="Fix overlap between extracted genes (only for animal mitochondria) [optional]",
        variable=var_delete_overlap
    )
    chk_delete_overlap.pack(anchor="center", padx=5, pady=(0, 8))

    chk_log_missing = ttk.Checkbutton(
        content4,
        text="Generate log of missing genes per sample (useful for mitogenomes) [optional]",
        variable=var_log_missing
    )
    chk_log_missing.pack(anchor="center", padx=5, pady=(0, 8))

    ttk.Label(
        content4,
        text="🔑 NCBI API key [optional]:",
        font=("Segoe UI", 14),
    ).pack(anchor="center", pady=(12, 5))
    entry_api_key = ttk.Entry(content4, width=52)
    entry_api_key.pack(anchor="center", padx=5, pady=(0, 8))
    ttk.Label(
        content4,
        text="Recommended for larger searches. The key is not written to the output log.",
        wraplength=520,
        justify="center",
        font=("Segoe UI", 12),
    ).pack(anchor="center", padx=5, pady=(0, 8))

    var_org_type.trace_add("write", update_checkbox_overlay_state)

    button_frame4 = ttk.Frame(content4)
    button_frame4.pack(pady=20)
    btn_previous4 = ttk.Button(button_frame4, text="◀ Previous", command=previous_frame)
    btn_previous4.pack(side="left", padx=10)
    btn_next4 = ttk.Button(button_frame4, text="Next ▶", command=next_frame)
    btn_next4.pack(side="left", padx=10)

    warning_extract_all_genes.pack(anchor="center", padx=5, pady=(10, 5))

    # ---- Frame 5 ----
    frame5, content5 = create_centered_frame()

    label_title = ttk.Label(content5, text="Step 5: Final options", font=("Segoe UI", 18, "bold"))
    label_title.pack(pady=(0, 20))

    ttk.Label(content5, text="📧 User e-mail to access NCBI database:").pack(anchor="center", pady=(10, 5))
    entry_email = ttk.Entry(content5, width=100)
    entry_email.pack(anchor="w", pady=(0, 10))

    ttk.Label(content5, text="🗂️ Set output folder:").pack(anchor="center", pady=(10, 5))
    frame_path = ttk.Frame(content5)
    frame_path.pack(fill="x", pady=(0, 10), padx=0)
    path_var = tk.StringVar()
    entry_path = ttk.Entry(frame_path, textvariable=path_var, width=50)
    entry_path.grid(row=0, column=0, sticky="ew", padx=(0, 5))
    btn_path = ttk.Button(frame_path, text="Select folder", command=choose_folder)
    btn_path.grid(row=0, column=1)
    frame_path.columnconfigure(0, weight=1)

    folder_name_var = tk.StringVar()
    ttk.Label(content5, text="📁 Folder name to create inside selected location [optional]:").pack(anchor="center", pady=(5, 5))
    frame_subfolder = ttk.Frame(content5)
    frame_subfolder.pack(fill="x", padx=0, pady=(0, 10))
    entry_folder_name = ttk.Entry(frame_subfolder, textvariable=folder_name_var, width=50)
    entry_folder_name.pack(side="left", fill="x", expand=True, padx=(0, 5))

    btn_previous5 = ttk.Button(content5, text="◀ Previous", command=previous_frame)
    btn_previous5.pack(pady=25)

    btn = tk.Button(
        content5,
        text="💾 Search and download sequences",
        command=start_download,
        bg="gray", fg="black",
        font=("Segoe UI", 14, "bold"),
        padx=10, pady=6
    )
    btn.pack(pady=(10, 0))

    # Inicia com o primeiro frame
    show_frame(0)

    # Mainloop
    root.mainloop()

except Exception:
    # Crash amigável ao abrir pelo Finder: mostra diálogo e grava log
    try:
        with open(LOG_FILE, "w", encoding="utf-8") as f:
            f.write("=== Pynnotate crash ===\n")
            traceback.print_exc(file=f)
    except Exception:
        pass
    try:
        _r = tk.Tk(); _r.withdraw()
        messagebox.showerror("Erro ao iniciar", f"Ocorreu um erro.\nVeja o log em:\n{LOG_FILE}")
    except Exception:
        pass
    raise
