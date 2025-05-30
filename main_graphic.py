from datetime import datetime
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from PIL import Image, ImageTk
import traceback
import os
import threading
import logging

from graphic.utils import show_gene_selector, show_field_selector, update_checkbox_overlay_state, update_extraction_options

from main.pynnotate import begin_search

def main():

    config = {
        "email": entry_email.get().strip(),
        "title": entry_title.get().strip(),
        "organisms": entry_org.get().strip(),
        "genes": entry_genes.get().strip(),
        "fields": entry_fields.get().strip(),
        "output_name": entry_path.get().strip(),
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
        "unique": var_unique_species.get(),
        "logmissing": var_log_missing.get(),
        "org_type": var_org_type.get(),
        "add_synonyms": entry_syn.get("1.0", tk.END).strip(),
        "additional": entry_add.get().strip()
    }
    print(entry_syn.get("1.0", tk.END).strip())
    begin_search(config, root=root)

root = tk.Tk()
if os.name == "nt":  # Windows
    root.iconbitmap("logo.ico")
else:
    img = tk.PhotoImage(file="@Pynnotate_v9/logo_transparente.png")
    root.iconphoto(True, img)
root.title("Pynnotate v.0.1")
root.geometry("600x700")

style = ttk.Style()
style.configure("TLabel", font=("Segoe UI", 14))
style.configure("TButton", font=("Segoe UI", 12), padding=6)
style.configure("TRadiobutton", font=("Segoe UI", 12))
style.configure("Bold.TCheckbutton", font=("Arial", 14, "bold"))

## Setting frames up

current_index = 0
frames = []

def show_frame(index):
    global current_index
    for f, _ in frames:
        f.grid_remove()
    frames[index][0].grid()
    current_index = index

def next_frame():
    global current_index
    if current_index + 1 < len(frames):
        show_frame(current_index + 1)

def previous_frame():
    global current_index
    if current_index - 1 >= 0:
        show_frame(current_index - 1)

def create_centered_frame():
    from tkinter import ttk

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

# Frame 0
frame0, content0 = create_centered_frame()

logo_img = Image.open("@Pynnotate_v9/logo_transparente.png")
logo_img = logo_img.resize((200, 200))
logo_photo = ImageTk.PhotoImage(logo_img)

logo_label = tk.Label(content0, image=logo_photo)
logo_label.image = logo_photo
logo_label.pack(pady=(0, 10))

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

## Frame 1
frame1, content1 = create_centered_frame()

label_title = ttk.Label(content1, text="Step 1: Choose genome type", font=("Segoe UI", 18, "bold"))
label_title.pack(pady=(0, 20))

var_org_type = tk.StringVar(value="animal_mito")

options_frame = ttk.Frame(content1)
options_frame.pack(pady=10)

r1 = ttk.Radiobutton(options_frame, text="Animal mitochondria", variable=var_org_type, value="animal_mito")
r1.pack(anchor="w", pady=6)

r2 = ttk.Radiobutton(options_frame, text="Plant mitochondria", variable=var_org_type, value="plant_mito")
r2.pack(anchor="w", pady=6)

r3 = ttk.Radiobutton(options_frame, text="Chloroplast", variable=var_org_type, value="plant_chloro")
r3.pack(anchor="w", pady=6)

r4 = ttk.Radiobutton(options_frame, text="Other", variable=var_org_type, value="other")
r4.pack(anchor="w", pady=6)

button_frame = ttk.Frame(content1)
btn_previous = ttk.Button(button_frame, text="◀ Previous", command=previous_frame)
btn_previous.pack(side="left", padx=10)

btn_next = ttk.Button(button_frame, text="Next ▶", command=next_frame)
btn_next.pack(side="left", padx=10)

button_frame.pack(pady=20)

## Frame 2
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
btn_genes = ttk.Button(frame_genes, text="Select genes", command=lambda: show_gene_selector(entry_genes, var_org_type.get()))
btn_genes.grid(row=0, column=1)
frame_genes.columnconfigure(0, weight=1)

ttk.Label(content2, text="🧫 Organism(s) [optional]:").pack(anchor="center", pady=(10, 5))
entry_org = ttk.Entry(content2, width=60)
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

button_frame = ttk.Frame(content2)
btn_previous = ttk.Button(button_frame, text="◀ Previous", command=previous_frame)
btn_previous.pack(side="left", padx=10)

btn_next_2 = ttk.Button(button_frame, text="Next ▶", command=next_frame)
btn_next_2.pack(side="left", padx=10)

button_frame.pack(pady=20)

## Frame 3
frame3, content3 = create_centered_frame()

label_title = ttk.Label(content3, text="Step 3: Formatting and filtering options for download", font=("Segoe UI", 18, "bold"))
label_title.pack(pady=(0, 20))

# Header fields selection
ttk.Label(content3, text="📋 Header fields (select):").pack(anchor="center", pady=(10, 5))
fields_frame = ttk.Frame(content3)
fields_frame.pack(fill="x", pady=(0, 10), padx=0)
entry_fields = ttk.Entry(fields_frame, width=50)
entry_fields.grid(row=0, column=0, sticky="ew", padx=(0, 5))
entry_fields.insert(0, "organism")
btn_fields = ttk.Button(fields_frame, text="Select", command=lambda: show_field_selector(entry_fields))
btn_fields.grid(row=0, column=1)
fields_frame.columnconfigure(0, weight=1) 

center_frame = ttk.Frame(content3)
center_frame.pack(pady=2)

# Include GenBank ID checkbox
var_id = tk.BooleanVar(value=True)
chk_id = ttk.Checkbutton(center_frame, text="Include GenBank ID in fasta header", variable=var_id)
chk_id.pack(anchor="w", pady=10)

# Unique species checkbox
var_unique_species = tk.BooleanVar(value=False)
chk_unique_species = ttk.Checkbutton(center_frame, text="Include only 1 individual per species", variable=var_unique_species)
chk_unique_species.pack(anchor="w", pady=4)

ttk.Label(content3, text='🔧 Add gene synonyms (YAML-style) [optional]\n[ex: { "CYTB": ["CYTOCHROME B", "CYT B"], "COI": ["CO1", "COX1"] }]:', wraplength=500, justify="center").pack(anchor="center", pady=(10, 5))
frame_syn = ttk.Frame(content3)
frame_syn.pack(fill="both", expand=True, pady=(0, 10), padx=0)
entry_syn = tk.Text(frame_syn, height=5, wrap="word")
entry_syn.pack(fill="both", expand=True)

# Minimum sequence size
frame_min_len = ttk.Frame(content3)
frame_min_len.pack(anchor="w", pady=(10, 10), fill="x")

ttk.Label(frame_min_len, text="🔢 Minimum sequence size (bp) [optional]:").grid(row=0, column=0, sticky="w", padx=(0, 10))

entry_min_len = ttk.Entry(frame_min_len, width=20)
entry_min_len.grid(row=0, column=1, sticky="ew")
frame_min_len.columnconfigure(1, weight=1)

# Maximum sequence size
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

# Frame 4
frame4, content4 = create_centered_frame()

label_title = ttk.Label(content4, text="Step 4: Additional options", font=("Segoe UI", 18, "bold"))
label_title.pack(pady=(0, 20))

# Variables
var_extract_all_genes = tk.BooleanVar(value=True)
var_delete_overlap = tk.BooleanVar(value=False)
var_log_missing = tk.BooleanVar(value=False)

# Main checkbutton
chk_extract_all_genes = ttk.Checkbutton(
    content4,
    text="🌟 Extract all annotated genes separately [optional]",
    variable=var_extract_all_genes,
    command=lambda: update_extraction_options(),
    style="Bold.TCheckbutton"
)
chk_extract_all_genes.pack(anchor="center", padx=5, pady=(0, 8))

# Additional checkbuttons (initially hidden)
chk_delete_overlap = ttk.Checkbutton(
    content4,
    text="Fix overlap between extracted genes [only for animal mitochondria]",
    variable=var_delete_overlap
)
chk_delete_overlap.pack(anchor="center", padx=5, pady=(0, 8))

chk_log_missing = ttk.Checkbutton(
    content4,
    text="Generate log of missing genes per sample (useful for mitogenomes)",
    variable=var_log_missing
)
chk_log_missing.pack(anchor="center", padx=5, pady=(0, 8))

var_org_type.trace_add("write", update_checkbox_overlay_state)

button_frame4 = ttk.Frame(content4)
button_frame4.pack(pady=20)

btn_previous4 = ttk.Button(button_frame4, text="◀ Previous", command=previous_frame)
btn_previous4.pack(side="left", padx=10)

btn_next4 = ttk.Button(button_frame4, text="Next ▶", command=next_frame)
btn_next4.pack(side="left", padx=10)

# Frame 5
frame5, content5 = create_centered_frame()

label_title = ttk.Label(content5, text="Step 5: Final options", font=("Segoe UI", 18, "bold"))
label_title.pack(pady=(0, 20))

ttk.Label(content5, text="📧 User e-mail to access NCBI database:").pack(anchor="center", pady=(10, 5))
entry_email = ttk.Entry(content5, width=100)
entry_email.pack(anchor="w", pady=(0, 10))

def choose_folder():
    folder = filedialog.askdirectory(title="Choose a folder to save the results")
    if folder:
        path_var.set(folder)

ttk.Label(content5, text="🗂️ Set output folder:").pack(anchor="center", pady=(10, 5))
frame_path = ttk.Frame(content5)
frame_path.pack(fill="x", pady=(0, 10), padx=0)
path_var = tk.StringVar()
entry_path = ttk.Entry(frame_path, textvariable=path_var, width=50)
entry_path.grid(row=0, column=0, sticky="ew", padx=(0, 5))
btn_path = ttk.Button(frame_path, text="Select folder", command=choose_folder)
btn_path.grid(row=0, column=1)
frame_path.columnconfigure(0, weight=1) 

btn_previous5 = ttk.Button(content5, text="◀ Previous", command=previous_frame)
btn_previous5.pack(pady=25)

btn = tk.Button(
    content5,
    text="💾 Search and download sequences",
    command=main,
    bg="gray",
    fg="black",
    font=("Arial", 14, "bold"),
    padx=10,
    pady=6
)

btn.pack(pady=(10, 0))

show_frame(0)

root.mainloop()