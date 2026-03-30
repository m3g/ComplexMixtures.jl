module GUI

using WGLMakie
using Bonito
using Bonito: DOM
using ComplexMixtures
using ComplexMixtures: Result, Options, AtomSelection, mddf,
    SoluteGroup, SolventGroup, ResidueContributions,
    contributions, overview, Overview,
    _set_clims_and_colorscale!
import PDBTools
using Statistics: mean, std

# ─────────────────────────────────────────────────────────────────────────
# Native file dialog helper (runs server-side)
# ─────────────────────────────────────────────────────────────────────────

function _pick_file(; title="Select file")
    path = ""
    try
        if Sys.islinux()
            if !isnothing(Sys.which("zenity"))
                path = strip(read(`zenity --file-selection --title=$title`, String))
            elseif !isnothing(Sys.which("kdialog"))
                path = strip(read(`kdialog --getopenfilename . --title $title`, String))
            end
        elseif Sys.isapple()
            path = strip(read(`osascript -e 'POSIX path of (choose file with prompt "Select file")'`, String))
        elseif Sys.iswindows()
            ps = """
            Add-Type -AssemblyName System.Windows.Forms
            \$d = New-Object System.Windows.Forms.OpenFileDialog
            \$d.Title = "$title"
            if (\$d.ShowDialog() -eq 'OK') { \$d.FileName }
            """
            path = strip(read(`powershell -NoProfile -Command $ps`, String))
        end
    catch
    end
    return path
end

function _pick_save_file(; title="Save file")
    path = ""
    try
        if Sys.islinux()
            if !isnothing(Sys.which("zenity"))
                path = strip(read(`zenity --file-selection --save --confirm-overwrite --title=$title`, String))
            elseif !isnothing(Sys.which("kdialog"))
                path = strip(read(`kdialog --getsavefilename . --title $title`, String))
            end
        elseif Sys.isapple()
            path = strip(read(`osascript -e 'POSIX path of (choose file name with prompt "Save file")'`, String))
        elseif Sys.iswindows()
            ps = """
            Add-Type -AssemblyName System.Windows.Forms
            \$d = New-Object System.Windows.Forms.SaveFileDialog
            \$d.Title = "$title"
            \$d.Filter = "JSON files (*.json)|*.json"
            if (\$d.ShowDialog() -eq 'OK') { \$d.FileName }
            """
            path = strip(read(`powershell -NoProfile -Command $ps`, String))
        end
    catch
    end
    return path
end

# ─────────────────────────────────────────────────────────────────────────
# CSS styling
# ─────────────────────────────────────────────────────────────────────────

const _CSS = DOM.style("""
body { font-family: 'Segoe UI', Arial, sans-serif; font-size: 13px; margin: 8px; }
.cm-title { font-size: 16px; font-weight: bold; color: #333; margin: 4px 0 8px 0; text-align: center; letter-spacing: 0.5px; }
.cm-main { display: grid; grid-template-columns: 340px 1fr; gap: 12px; height: calc(100vh - 60px); }
.cm-sidebar { overflow-y: auto; padding: 10px; border: 1px solid #ddd; border-radius: 6px; background: #fafafa; }
.cm-plots { overflow-y: auto; overflow-x: auto; padding: 8px; width: fit-content; }
.cm-section-title { font-weight: bold; font-size: 14px; margin: 10px 0 4px 0; color: #333; }
.cm-row { display: flex; align-items: center; gap: 6px; margin: 3px 0; }
.cm-row label { min-width: 100px; text-align: right; font-size: 12px; color: #555; }
.cm-row-compact { display: flex; align-items: center; gap: 6px; margin: 1px 0 3px 0; padding-left: 6px; }
.cm-row-compact select { max-width: 130px; font-size: 11px; }
.cm-row-compact input { width: 50px; font-size: 11px; }
.cm-row input, .cm-row select { font-size: 12px; padding: 3px 6px; }
.cm-row input[type=text], .cm-row input[type=textfield] { flex: 1; }
.cm-row input[type=number] { width: 70px; }
.cm-required label { color: #b00; font-weight: bold; }
.cm-browse-btn button { padding: 2px 6px; font-size: 11px; min-width: 28px; }
button { font-size: 12px !important; padding: 3px 8px !important; }
input, select, textarea { font-size: 12px !important; }
.cm-grp-list { margin: 4px 0; border: 1px solid #e0e0e0; border-radius: 4px; padding: 4px 8px; overflow-y: auto; background: #fff; }
.cm-grp-item { display: flex; align-items: center; gap: 6px; padding: 2px 0; font-size: 12px; }
.cm-tab2-body { display: flex; gap: 10px; align-items: stretch; }
.cm-grp-panel { width: 220px; flex-shrink: 0; border: 1px solid #ddd; border-radius: 4px; padding: 8px; background: #fafafa; display: flex; flex-direction: column; }
.cm-export-section { margin-top: auto; padding-top: 8px; border-top: 1px solid #ddd; }
.cm-export-field { display: flex; flex-direction: column; gap: 2px; margin: 3px 0; }
.cm-export-field label { color: #555; font-size: 11px; }
.cm-export-field input[type=text], .cm-export-field input[type=textfield] { width: 100%; font-size: 11px !important; box-sizing: border-box; }
.cm-export-field-inline { display: flex; align-items: center; gap: 6px; margin: 3px 0; }
.cm-export-field-inline label { white-space: nowrap; color: #555; font-size: 11px; min-width: 56px; }
.cm-export-field-inline select { flex: 1; font-size: 11px !important; }
.cm-export-btns { display: flex; gap: 4px; margin: 4px 0; flex-wrap: wrap; }
.cm-export-btns button { flex: 1; font-size: 10px !important; padding: 2px 4px !important; }
.cm-grp-panel-title { font-weight: bold; font-size: 12px; color: #333; margin-bottom: 4px; }
.cm-comp-tabs { display: flex; gap: 0; margin-bottom: 6px; border-bottom: 2px solid #ccc; width: fit-content; }
.cm-comp-tab button { font-size: 11px !important; padding: 3px 8px !important; border: 1px solid #ccc; border-bottom: none; border-radius: 3px 3px 0 0; background: #f0f0f0; cursor: pointer; margin-right: 2px; width: fit-content !important; }
.cm-comp-tab.active button { background: #3970d1 !important; color: white !important; font-weight: bold; border-color: #3970d1; }
.cm-fig-wrap { display: flex; flex-direction: column; }
.cm-fig-row { display: flex; align-items: stretch; }
.cm-ylabel { writing-mode: vertical-rl; transform: rotate(180deg); font-size: 12px; color: #444; text-align: center; padding: 0 3px; white-space: nowrap; }
.cm-xlabel { text-align: center; font-size: 12px; color: #444; margin: 1px 0 3px 0; }
.cm-plot-title { text-align: center; font-size: 14px; font-weight: bold; color: #222; margin: 6px 0 1px 0; }
.cm-row-left label { min-width: auto !important; text-align: left !important; font-weight: normal; }
.cm-tabs { display: flex; gap: 4px; margin-bottom: 0px; }
.cm-tab-btn { padding: 6px 16px; border: 1px solid #ccc; border-radius: 4px 4px 0 0;
              cursor: pointer; background: #eee; font-size: 12px; border-bottom: none; user-select: none; }
.cm-tab-btn.active { background: #3970d1; color: white; font-weight: bold; }
.cm-tab-content { display: none; border: 1px solid #ccc; border-radius: 0 4px 4px 4px; padding: 8px; background: #fff; }
.cm-tab-content.active { display: block; }
.cm-status { font-size: 11px; color: #888; padding: 4px 8px; border-top: 1px solid #ddd; }
.cm-toggle-row { display: flex; align-items: center; justify-content: center; gap: 10px; margin: 3px 0; font-size: 11px; color: #777; }
.cm-advanced-sep { border: none; border-top: 1px solid #e0e0e0; margin: 6px 0 4px 0; }
.cm-lims-row { display: flex; align-items: center; gap: 3px; margin: 3px 0; font-size: 10px; flex-wrap: wrap; }
.cm-lims-row span { color: #666; font-size: 10px; }
.cm-lims-row input[type=number] { width: 50px; font-size: 10px; padding: 1px 3px; }
.cm-lims-row button { font-size: 10px; padding: 2px 6px; }
.cm-overview { font-family: 'Cascadia Mono', 'Consolas', monospace; font-size: 11px; white-space: pre-wrap; padding: 6px;
               border: 1px solid #ddd; border-radius: 4px; background: #fff; max-height: 400px; overflow-y: auto; }
.cm-overview table { border-collapse: collapse; width: 100%; font-size: 11px; }
.cm-overview td { padding: 2px 6px; }
.cm-overview td:first-child { color: #555; text-align: right; white-space: nowrap; }
.cm-overview td:last-child { font-weight: 500; }
.cm-overview .cm-ov-section { font-weight: bold; font-size: 12px; color: #333; padding: 6px 0 2px 0; }
.cm-sidebar-tabs { display: flex; gap: 2px; margin-bottom: 4px; }
.cm-sidebar-tab { padding: 4px 12px; border: 1px solid #ccc; border-radius: 4px 4px 0 0;
                  cursor: pointer; background: #eee; font-size: 11px; border-bottom: none; user-select: none; }
.cm-sidebar-tab.active { background: #555; color: white; font-weight: bold; }
.cm-sidebar-content { display: none; }
.cm-sidebar-content.active { display: block; }
""")

# ─────────────────────────────────────────────────────────────────────────
# Format overview as HTML
# ─────────────────────────────────────────────────────────────────────────

function _overview_text(R::Result)
    ov = overview(R)
    fmt(x) = string(round(x; sigdigits=6))
    ifar = trunc(Int, ov.R.nbins - 1.0 / ov.R.files[1].options.binstep)
    lr_mddf_mean = fmt(mean(ov.R.mddf[ifar:ov.R.nbins]))
    lr_mddf_std = fmt(std(ov.R.mddf[ifar:ov.R.nbins]))
    lr_rdf_mean = fmt(mean(ov.R.rdf[ifar:ov.R.nbins]))
    lr_rdf_std = fmt(std(ov.R.rdf[ifar:ov.R.nbins]))
    bulk_str = ComplexMixtures._bulk_range_from_R(ov.R)
    files_str = join(
        ["  $(normpath(ov.R.files[i].filename)) (w=$(ov.R.weights[i]))" for i in eachindex(ov.R.files)],
        "\n"
    )
    return """
Solvent
  Concentration:      $(fmt(ov.density.solvent)) mol/L
  Molar volume:       $(fmt(ov.solvent_molar_volume)) cm³/mol
  Bulk concentration: $(fmt(ov.density.solvent_bulk)) mol/L
  Bulk molar volume:  $(fmt(ov.solvent_molar_volume_bulk)) cm³/mol

Solute
  Concentration:      $(fmt(ov.density.solute)) mol/L
  Partial molar vol.: $(fmt(ov.solute_molar_volume)) cm³/mol

System
  Bulk range:         $(bulk_str)
  Domain molar vol.:  $(fmt(ov.domain_molar_volume)) cm³/mol
  Auto-correlation:   $(ov.R.autocorrelation)

Convergence
  Long-range MDDF:    $(lr_mddf_mean) ± $(lr_mddf_std)
  Long-range RDF:     $(lr_rdf_mean) ± $(lr_rdf_std)

Trajectory files
$(files_str)
"""
end

# ─────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────

"""
    ComplexMixtures.gui(; pdbfile=nothing, result=nothing, port=9384)

Launch the Bonito/WGLMakie web-based interface for ComplexMixtures analysis.
Opens a browser at `http://localhost:\$port`.
"""
function ComplexMixtures.gui(;
    pdbfile::Union{<:AbstractString,Nothing}=nothing,
    result::Union{<:AbstractString,Nothing}=nothing,
    port::Int=9384,
)
    app = App() do session::Bonito.Session
        # ── State ──────────────────────────────────────────────────────
        result_obs = Observable{Union{Nothing,Result}}(nothing)
        atoms_obs = Observable{Union{Nothing,Vector{PDBTools.Atom}}}(nothing)
        status_obs = Observable("Ready")
        overview_obs = Observable("No results loaded.")
        computing = Observable(false)

        # ── Preload if provided ────────────────────────────────────────
        if !isnothing(result)
            result_obs[] = ComplexMixtures.load(result)
        end
        if !isnothing(pdbfile)
            atoms_obs[] = PDBTools.read_pdb(pdbfile)
        end

        # ══════════════════════════════════════════════════════════════
        # LEFT PANEL — with two subtabs: "Setup / Run" and "Load / Results"
        # ══════════════════════════════════════════════════════════════

        # ── Shared PDB widget (appears in both tabs) ─────────────────
        tf_pdb = Bonito.TextField(isnothing(pdbfile) ? "" : pdbfile)
        btn_browse_pdb = Bonito.Button("📂")

        # ── Setup / Run subtab widgets ─────────────────────────────────
        tf_traj = Bonito.TextField("")
        btn_browse_traj = Bonito.Button("📂")
        tf_solute = Bonito.TextField("protein")
        dd_solute_opt = Bonito.Dropdown(["nmols", "natomspermol"])
        tf_solute_n = Bonito.TextField("1")
        tf_solvent = Bonito.TextField("resname TMAO")
        dd_solvent_opt = Bonito.Dropdown(["natomspermol", "nmols"])
        tf_solvent_n = Bonito.TextField("14")

        # Options
        tf_bulk_min = Bonito.NumberInput(8.0)
        tf_bulk_max = Bonito.NumberInput(12.0)
        tf_binstep = Bonito.NumberInput(0.02)
        tf_stride = Bonito.NumberInput(1.0)
        tf_firstframe = Bonito.NumberInput(1.0)
        tf_lastframe = Bonito.NumberInput(-1.0)
        tf_nrandom = Bonito.NumberInput(10.0)
        tf_seed = Bonito.NumberInput(321.0)

        cb_silent = Bonito.Checkbox(false)
        cb_gc = Bonito.Checkbox(true)
        cb_stablerng = Bonito.Checkbox(false)

        btn_run = Bonito.Button("Run MDDF")

        setup_content = DOM.div(
            DOM.div(class="cm-row", DOM.label("Trajectory:"), tf_traj, DOM.div(class="cm-browse-btn", btn_browse_traj)),
            DOM.div(class="cm-section-title", "Atom Selections"),
            DOM.div(class="cm-row cm-row-left", DOM.label("Solute:"), tf_solute),
            DOM.div(class="cm-row-compact", dd_solute_opt, tf_solute_n),
            DOM.div(class="cm-row cm-row-left", DOM.label("Solvent:"), tf_solvent),
            DOM.div(class="cm-row-compact", dd_solvent_opt, tf_solvent_n),
            DOM.div(class="cm-section-title", "Options"),
            DOM.div(class="cm-row", DOM.label("bulk from (Å):"), tf_bulk_min),
            DOM.div(class="cm-row", DOM.label("bulk to (Å):"), tf_bulk_max),
            DOM.div(class="cm-row", DOM.label("binstep (Å):"), tf_binstep),
            DOM.div(class="cm-row", DOM.label("stride:"), tf_stride),
            DOM.div(class="cm-row", DOM.label("first frame:"), tf_firstframe),
            DOM.div(class="cm-row", DOM.label("last frame:"), tf_lastframe),
            DOM.div(class="cm-row", DOM.label("n_random:"), tf_nrandom),
            DOM.div(class="cm-row", DOM.label("seed:"), tf_seed),
            DOM.div(style="text-align:center; margin: 8px 0;", btn_run),
            DOM.hr(class="cm-advanced-sep"),
            DOM.div(class="cm-toggle-row",
                cb_silent, "Silent", cb_gc, "GC", cb_stablerng, "StableRNG"),
        )

        # ── Load / Results subtab widgets ──────────────────────────────
        tf_json = Bonito.TextField(isnothing(result) ? "" : result)
        btn_browse_json = Bonito.Button("📂")
        btn_load = Bonito.Button("Load JSON")

        tf_save_path = Bonito.TextField("")
        btn_browse_save = Bonito.Button("📂")
        btn_save = Bonito.Button("Save Result to JSON")

        results_content = DOM.div(
            DOM.div(class="cm-section-title", "Load Result"),
            DOM.div(class="cm-row cm-required", DOM.label("JSON file:"), tf_json, DOM.div(class="cm-browse-btn", btn_browse_json)),
            DOM.div(style="text-align:center; margin: 6px 0;", btn_load),
            DOM.div(class="cm-section-title", "Save Result"),
            DOM.div(class="cm-row", DOM.label("Save as:"), tf_save_path, DOM.div(class="cm-browse-btn", btn_browse_save)),
            DOM.div(style="text-align:center; margin: 6px 0;", btn_save),
            DOM.div(class="cm-section-title", "Result Overview"),
            DOM.div(class="cm-overview", overview_obs),
        )

        # ── Sidebar tabs JS switching ──────────────────────────────────
        sidebar_tabs_js = DOM.script("""
        document.addEventListener('DOMContentLoaded', function() {
            setTimeout(function() {
                var sbtns = document.querySelectorAll('.cm-sidebar-tab');
                var stabs = document.querySelectorAll('.cm-sidebar-content');
                sbtns.forEach(function(b, i) {
                    b.addEventListener('click', function() {
                        sbtns.forEach(function(bb) { bb.classList.remove('active'); });
                        stabs.forEach(function(tt) { tt.classList.remove('active'); });
                        b.classList.add('active');
                        stabs[i].classList.add('active');
                    });
                });
            }, 500);
        });
        """)

        sidebar = DOM.div(class="cm-sidebar",
            DOM.div(class="cm-section-title", "Structure File"),
            DOM.div(class="cm-row cm-required", DOM.label("PDB / mmCIF:"), tf_pdb, DOM.div(class="cm-browse-btn", btn_browse_pdb)),
            DOM.div(class="cm-sidebar-tabs",
                DOM.div(class="cm-sidebar-tab active", "Setup / Run"),
                DOM.div(class="cm-sidebar-tab", "Load / Results"),
            ),
            DOM.div(class="cm-sidebar-content active", setup_content),
            DOM.div(class="cm-sidebar-content", results_content),
            sidebar_tabs_js,
        )

        # ══════════════════════════════════════════════════════════════
        # RIGHT PANEL — plots with tabs
        # ══════════════════════════════════════════════════════════════

        # ── Tab 1: MDDF & KB ──────────────────────────────────────────
        fig1_mddf = Figure(; size=(900, 295))
        ax_mddf = Axis(fig1_mddf[1, 1]; xticklabelsize=11, yticklabelsize=11)
        fig1_kb = Figure(; size=(900, 295))
        ax_kb = Axis(fig1_kb[1, 1]; xticklabelsize=11, yticklabelsize=11)

        # Tab 1 limits (below the plots)
        tf_mddf_xmin = Bonito.NumberInput(0.0)
        tf_mddf_xmax = Bonito.NumberInput(10.0)
        tf_mddf_ymin = Bonito.NumberInput(0.0)
        tf_mddf_ymax = Bonito.NumberInput(5.0)
        btn_mddf_lims = Bonito.Button("Apply")
        tf_kb_xmin = Bonito.NumberInput(0.0)
        tf_kb_xmax = Bonito.NumberInput(10.0)
        tf_kb_ymin = Bonito.NumberInput(-1.0)
        tf_kb_ymax = Bonito.NumberInput(1.0)
        btn_kb_lims = Bonito.Button("Apply")

        tab1_limits = DOM.div(
            DOM.div(class="cm-lims-row",
                DOM.span("MDDF:"),
                DOM.span("x:"), tf_mddf_xmin, DOM.span("–"), tf_mddf_xmax,
                DOM.span("y:"), tf_mddf_ymin, DOM.span("–"), tf_mddf_ymax,
                btn_mddf_lims),
            DOM.div(class="cm-lims-row",
                DOM.span("KB:"),
                DOM.span("x:"), tf_kb_xmin, DOM.span("–"), tf_kb_xmax,
                DOM.span("y:"), tf_kb_ymin, DOM.span("–"), tf_kb_ymax,
                btn_kb_lims),
        )

        # ── Tab 2 & 3: Group Contributions (Solute / Solvent) ────────────
        MAX_GROUPS = 10

        # Solute group tab figures
        fig_sol_mddf = Figure(; size=(900, 295))
        ax_sol_mddf = Axis(fig_sol_mddf[1, 1]; xticklabelsize=11, yticklabelsize=11)
        fig_sol_cn = Figure(; size=(900, 295))
        ax_sol_cn = Axis(fig_sol_cn[1, 1]; xticklabelsize=11, yticklabelsize=11)

        # Solvent group tab figures
        fig_slv_mddf = Figure(; size=(900, 295))
        ax_slv_mddf = Axis(fig_slv_mddf[1, 1]; xticklabelsize=11, yticklabelsize=11)
        fig_slv_cn = Figure(; size=(900, 295))
        ax_slv_cn = Axis(fig_slv_cn[1, 1]; xticklabelsize=11, yticklabelsize=11)

        # Solute group state
        tf_comp_sol = Bonito.TextField("protein")
        tf_newgrp_sol = Bonito.TextField("resname ARG")
        btn_addgrp_sol = Bonito.Button("Add")
        btn_rmgrp_sol = Bonito.Button("Remove unchecked")
        grp_names_sol  = [Observable("") for _ in 1:MAX_GROUPS]
        grp_labels_sol = [Observable("") for _ in 1:MAX_GROUPS]
        grp_checks_sol = [Bonito.Checkbox(true) for _ in 1:MAX_GROUPS]
        grp_active_sol = Observable(0)
        grp_total_sol = Bonito.Checkbox(true)

        # Solvent group state
        tf_comp_slv = Bonito.TextField("resname TMAO")
        tf_newgrp_slv = Bonito.TextField("element O")
        btn_addgrp_slv = Bonito.Button("Add")
        btn_rmgrp_slv = Bonito.Button("Remove unchecked")
        grp_names_slv  = [Observable("") for _ in 1:MAX_GROUPS]
        grp_labels_slv = [Observable("") for _ in 1:MAX_GROUPS]
        grp_checks_slv = [Bonito.Checkbox(true) for _ in 1:MAX_GROUPS]
        grp_active_slv = Observable(0)
        grp_total_slv = Bonito.Checkbox(true)

        checklist_sol_dom = map(grp_active_sol) do n
            rows = Any[DOM.div(class="cm-grp-item", grp_total_sol, DOM.span("Total MDDF"))]
            for i in 1:n
                push!(rows, DOM.div(class="cm-grp-item", grp_checks_sol[i], DOM.span(grp_labels_sol[i][])))
            end
            DOM.div(class="cm-grp-list", rows...)
        end

        checklist_slv_dom = map(grp_active_slv) do n
            rows = Any[DOM.div(class="cm-grp-item", grp_total_slv, DOM.span("Total MDDF"))]
            for i in 1:n
                push!(rows, DOM.div(class="cm-grp-item", grp_checks_slv[i], DOM.span(grp_labels_slv[i][])))
            end
            DOM.div(class="cm-grp-list", rows...)
        end

        # Group limits (shared)
        tf_grp_xmin = Bonito.NumberInput(0.0)
        tf_grp_xmax = Bonito.NumberInput(10.0)
        tf_grp_ymin = Bonito.NumberInput(0.0)
        tf_grp_ymax = Bonito.NumberInput(5.0)
        btn_grp_lims_sol = Bonito.Button("Apply")
        btn_grp_lims_slv = Bonito.Button("Apply")

        # Export controls - Solute
        tf_export_sol      = Bonito.TextField("solute_group_contributions")
        dd_export_fmt_sol  = Bonito.Dropdown(["svg", "png", "pdf"])
        btn_export_mddf_sol = Bonito.Button("Export MDDF plot")
        btn_export_cn_sol   = Bonito.Button("Export CN plot")
        btn_export_csv_sol  = Bonito.Button("Export data (CSV)")
        last_sol_data = Observable{Any}(nothing)

        # Export controls - Solvent
        tf_export_slv      = Bonito.TextField("solvent_group_contributions")
        dd_export_fmt_slv  = Bonito.Dropdown(["svg", "png", "pdf"])
        btn_export_mddf_slv = Bonito.Button("Export MDDF plot")
        btn_export_cn_slv   = Bonito.Button("Export CN plot")
        btn_export_csv_slv  = Bonito.Button("Export data (CSV)")
        last_slv_data = Observable{Any}(nothing)

        tab_sol_body = DOM.div(class="cm-tab2-body",
            DOM.div(
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Solute MDDF Group Contributions"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "MDDF(r)"),
                        fig_sol_mddf,
                    ),
                ),
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Solute Coord. Number Group Contributions"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "Coordination number"),
                        fig_sol_cn,
                    ),
                ),
                DOM.div(class="cm-xlabel", "r (Angstrom)"),
                DOM.div(class="cm-lims-row",
                    DOM.span("Limits:"),
                    DOM.span("x:"), tf_grp_xmin, DOM.span("–"), tf_grp_xmax,
                    DOM.span("y:"), tf_grp_ymin, DOM.span("–"), tf_grp_ymax,
                    btn_grp_lims_sol),
            ),
            DOM.div(class="cm-grp-panel",
                DOM.div(class="cm-row cm-row-left", DOM.label("Solute:"), tf_comp_sol),
                DOM.div(class="cm-grp-panel-title", "Groups"),
                DOM.div(class="cm-lims-row", tf_newgrp_sol, btn_addgrp_sol),
                checklist_sol_dom,
                DOM.div(style="margin: 4px 0;", btn_rmgrp_sol),
                DOM.div(class="cm-export-section",
                    DOM.div(class="cm-grp-panel-title", "Export"),
                    DOM.div(class="cm-export-field", DOM.label("File name:"), tf_export_sol),
                    DOM.div(class="cm-export-field-inline", DOM.label("Format:"), dd_export_fmt_sol),
                    DOM.div(class="cm-export-btns", btn_export_mddf_sol, btn_export_cn_sol),
                    DOM.div(class="cm-export-btns", btn_export_csv_sol),
                ),
            ),
        )

        tab_slv_body = DOM.div(class="cm-tab2-body",
            DOM.div(
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Solvent MDDF Group Contributions"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "MDDF(r)"),
                        fig_slv_mddf,
                    ),
                ),
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Solvent Coord. Number Group Contributions"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "Coordination number"),
                        fig_slv_cn,
                    ),
                ),
                DOM.div(class="cm-xlabel", "r (Angstrom)"),
                DOM.div(class="cm-lims-row",
                    DOM.span("Limits:"),
                    DOM.span("x:"), tf_grp_xmin, DOM.span("–"), tf_grp_xmax,
                    DOM.span("y:"), tf_grp_ymin, DOM.span("–"), tf_grp_ymax,
                    btn_grp_lims_slv),
            ),
            DOM.div(class="cm-grp-panel",
                DOM.div(class="cm-row cm-row-left", DOM.label("Solvent:"), tf_comp_slv),
                DOM.div(class="cm-grp-panel-title", "Groups"),
                DOM.div(class="cm-lims-row", tf_newgrp_slv, btn_addgrp_slv),
                checklist_slv_dom,
                DOM.div(style="margin: 4px 0;", btn_rmgrp_slv),
                DOM.div(class="cm-export-section",
                    DOM.div(class="cm-grp-panel-title", "Export"),
                    DOM.div(class="cm-export-field", DOM.label("File name:"), tf_export_slv),
                    DOM.div(class="cm-export-field-inline", DOM.label("Format:"), dd_export_fmt_slv),
                    DOM.div(class="cm-export-btns", btn_export_mddf_slv, btn_export_cn_slv),
                    DOM.div(class="cm-export-btns", btn_export_csv_slv),
                ),
            ),
        )

        # ── Tab 3: Residue Contributions ──────────────────────────────
        fig3_mddf = Figure(; size=(900, 295))
        ax_rc_mddf = Axis(fig3_mddf[1, 1]; xticklabelsize=9, yticklabelsize=11)
        fig3_cn = Figure(; size=(900, 295))
        ax_rc_cn = Axis(fig3_cn[1, 1]; xticklabelsize=9, yticklabelsize=11)
        tf_rc_sel = Bonito.TextField("protein")
        tf_dmin = Bonito.NumberInput(1.5)
        tf_dmax = Bonito.NumberInput(3.5)
        btn_rc_update = Bonito.Button("Update")

        tab3_controls = DOM.div(
            DOM.div(class="cm-row", DOM.label("Selection:"), tf_rc_sel),
            DOM.div(class="cm-lims-row", style="justify-content: center;",
                DOM.span("dmin:"), tf_dmin,
                DOM.span("dmax:"), tf_dmax,
                btn_rc_update,
            ),
        )

        # ── Plot tabs JS switching ────────────────────────────────────
        tab_switch_js = DOM.script("""
        document.addEventListener('DOMContentLoaded', function() {
            setTimeout(function() {
                var btns = document.querySelectorAll('.cm-tab-btn');
                var tabs = document.querySelectorAll('.cm-tab-content');
                btns.forEach(function(b, i) {
                    b.addEventListener('click', function() {
                        btns.forEach(function(bb) { bb.classList.remove('active'); });
                        tabs.forEach(function(tt) { tt.classList.remove('active'); });
                        b.classList.add('active');
                        tabs[i].classList.add('active');
                    });
                });
            }, 500);
        });
        """)

        plots_panel = DOM.div(class="cm-plots",
            DOM.div(class="cm-tabs",
                DOM.div(class="cm-tab-btn active", "MDDF & KB"),
                DOM.div(class="cm-tab-btn", "Solute group contributions"),
                DOM.div(class="cm-tab-btn", "Solvent group contributions"),
                DOM.div(class="cm-tab-btn", "Residue Contributions"),
            ),
            DOM.div(class="cm-tab-content active",
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "MDDF"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "MDDF(r)"),
                        fig1_mddf,
                    ),
                ),
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Kirkwood-Buff Integral"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "KB (L/mol)"),
                        fig1_kb,
                    ),
                ),
                DOM.div(class="cm-xlabel", "r (Angstrom)"),
                tab1_limits,
            ),
            DOM.div(class="cm-tab-content", tab_sol_body),
            DOM.div(class="cm-tab-content", tab_slv_body),
            DOM.div(class="cm-tab-content", tab3_controls,
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Residue Contributions to MDDF"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "r (Angstrom)"),
                        fig3_mddf,
                    ),
                ),
                DOM.div(class="cm-fig-wrap",
                    DOM.div(class="cm-plot-title", "Residue Contributions to Coordination Number"),
                    DOM.div(class="cm-fig-row",
                        DOM.span(class="cm-ylabel", "r (Angstrom)"),
                        fig3_cn,
                    ),
                ),
                DOM.div(class="cm-xlabel", "Residue"),
            ),
            tab_switch_js,
        )

        status_bar = DOM.div(class="cm-status", status_obs)

        # ══════════════════════════════════════════════════════════════
        # CALLBACKS
        # ══════════════════════════════════════════════════════════════

        # ── Browse buttons ─────────────────────────────────────────────
        on(btn_browse_pdb.value) do _
            path = _pick_file(; title="Select PDB / mmCIF file")
            if !isempty(path)
                tf_pdb.value[] = path
            end
        end
        on(btn_browse_traj.value) do _
            path = _pick_file(; title="Select trajectory file")
            isempty(path) || (tf_traj.value[] = path)
        end
        on(btn_browse_json.value) do _
            path = _pick_file(; title="Select JSON result file")
            isempty(path) || (tf_json.value[] = path)
        end
        on(btn_browse_save.value) do _
            path = _pick_save_file(; title="Save result JSON")
            isempty(path) || (tf_save_path.value[] = path)
        end

        # ── Group management ───────────────────────────────────────────
        function _grp_label(sel, comp_sel, at)
            if at === nothing
                return sel
            end
            combined = isempty(strip(comp_sel)) ? sel : "($comp_sel) and ($sel)"
            n_atoms = try length(PDBTools.select(at, combined)) catch; -1 end
            n_atoms < 0 ? sel : "$sel ($n_atoms atoms)"
        end

        on(btn_addgrp_sol.value) do _
            sel = strip(String(tf_newgrp_sol.value[]))
            isempty(sel) && return
            n = grp_active_sol[]
            n >= MAX_GROUPS && (status_obs[] = "Maximum $MAX_GROUPS groups reached"; return)
            for i in 1:n; grp_names_sol[i][] == sel && return; end
            grp_names_sol[n + 1][] = sel
            grp_labels_sol[n + 1][] = _grp_label(sel, String(tf_comp_sol.value[]), atoms_obs[])
            grp_checks_sol[n + 1].value[] = true
            grp_active_sol[] = n + 1
        end
        on(btn_rmgrp_sol.value) do _
            n = grp_active_sol[]
            kept_names = String[]; kept_labels = String[]
            for i in 1:n
                if grp_checks_sol[i].value[]
                    push!(kept_names, grp_names_sol[i][])
                    push!(kept_labels, grp_labels_sol[i][])
                end
            end
            for i in eachindex(kept_names)
                grp_names_sol[i][] = kept_names[i]
                grp_labels_sol[i][] = kept_labels[i]
                grp_checks_sol[i].value[] = true
            end
            grp_active_sol[] = length(kept_names)
        end
        on(btn_addgrp_slv.value) do _
            sel = strip(String(tf_newgrp_slv.value[]))
            isempty(sel) && return
            n = grp_active_slv[]
            n >= MAX_GROUPS && (status_obs[] = "Maximum $MAX_GROUPS groups reached"; return)
            for i in 1:n; grp_names_slv[i][] == sel && return; end
            grp_names_slv[n + 1][] = sel
            grp_labels_slv[n + 1][] = _grp_label(sel, String(tf_comp_slv.value[]), atoms_obs[])
            grp_checks_slv[n + 1].value[] = true
            grp_active_slv[] = n + 1
        end
        on(btn_rmgrp_slv.value) do _
            n = grp_active_slv[]
            kept_names = String[]; kept_labels = String[]
            for i in 1:n
                if grp_checks_slv[i].value[]
                    push!(kept_names, grp_names_slv[i][])
                    push!(kept_labels, grp_labels_slv[i][])
                end
            end
            for i in eachindex(kept_names)
                grp_names_slv[i][] = kept_names[i]
                grp_labels_slv[i][] = kept_labels[i]
                grp_checks_slv[i].value[] = true
            end
            grp_active_slv[] = length(kept_names)
        end

        # ── Load JSON ──────────────────────────────────────────────────
        on(btn_load.value) do _
            path = String(tf_json.value[])
            if isempty(path) || !isfile(path)
                status_obs[] = "File not found: $path"; return
            end
            try
                result_obs[] = ComplexMixtures.load(path)
                pdb_path = String(tf_pdb.value[])
                if !isempty(pdb_path) && isfile(pdb_path)
                    atoms_obs[] = PDBTools.read_pdb(pdb_path)
                end
                status_obs[] = "Loaded: $path"
            catch e
                status_obs[] = "Error: $(sprint(showerror, e))"
            end
        end

        # ── Save Result ────────────────────────────────────────────────
        on(btn_save.value) do _
            R = result_obs[]
            if R === nothing
                status_obs[] = "No results to save"; return
            end
            path = String(tf_save_path.value[])
            if isempty(path)
                status_obs[] = "Enter a file path for saving"; return
            end
            try
                ComplexMixtures.save(path, R)
                status_obs[] = "Saved: $path"
            catch e
                status_obs[] = "Error saving: $(sprint(showerror, e))"
            end
        end

        # ── Run MDDF ───────────────────────────────────────────────────
        on(btn_run.value) do _
            computing[] && return
            pdb_path = String(tf_pdb.value[])
            traj_path = String(tf_traj.value[])
            if isempty(pdb_path) || !isfile(pdb_path)
                status_obs[] = "PDB file not found: $pdb_path"; return
            end
            if isempty(traj_path) || !isfile(traj_path)
                status_obs[] = "Trajectory not found: $traj_path"; return
            end
            computing[] = true
            status_obs[] = "Computing MDDF…"
            Threads.@spawn begin
                try
                    all_atoms = PDBTools.read_pdb(pdb_path)
                    atoms_obs[] = all_atoms
                    solute_sel = String(tf_solute.value[])
                    solvent_sel = String(tf_solvent.value[])
                    solute_atoms = PDBTools.select(all_atoms, solute_sel)
                    solvent_atoms = PDBTools.select(all_atoms, solvent_sel)
                    solute_opt = String(dd_solute_opt.value[])
                    solvent_opt = String(dd_solvent_opt.value[])
                    solute_n = parse(Int, strip(String(tf_solute_n.value[])))
                    solvent_n = parse(Int, strip(String(tf_solvent_n.value[])))
                    solute_kw = solute_opt == "nmols" ? (; nmols=solute_n) : (; natomspermol=solute_n)
                    solvent_kw = solvent_opt == "nmols" ? (; nmols=solvent_n) : (; natomspermol=solvent_n)
                    solute = AtomSelection(solute_atoms; solute_kw...)
                    solvent = AtomSelection(solvent_atoms; solvent_kw...)
                    opts = Options(;
                        binstep=Float64(tf_binstep.value[]),
                        bulk_range=(Float64(tf_bulk_min.value[]), Float64(tf_bulk_max.value[])),
                        firstframe=round(Int, tf_firstframe.value[]),
                        lastframe=round(Int, tf_lastframe.value[]),
                        stride=round(Int, tf_stride.value[]),
                        n_random_samples=round(Int, tf_nrandom.value[]),
                        seed=round(Int, tf_seed.value[]),
                        silent=false,
                        GC=true,
                        StableRNG=false,
                    )
                    R = mddf(traj_path, solute, solvent, opts)
                    result_obs[] = R
                    status_obs[] = "Computation complete!"
                catch e
                    status_obs[] = "Error: $(sprint(showerror, e))"
                finally
                    computing[] = false
                end
            end
        end

        # ── Tab 1: auto-update on result change ───────────────────────
        on(result_obs) do R
            R === nothing && return
            empty!(ax_mddf)
            empty!(ax_kb)
            for c in copy(fig1_mddf.content)
                c isa Legend && delete!(c)
            end
            for c in copy(fig1_kb.content)
                c isa Legend && delete!(c)
            end
            lines!(ax_mddf, R.d, R.mddf; color=:dodgerblue, linewidth=1.5, label="MDDF")
            hlines!(ax_mddf, [1.0]; color=:gray60, linestyle=:dash)
            axislegend(ax_mddf; position=:rt, labelsize=10)
            lines!(ax_kb, R.d, R.kb ./ 1000; color=:orangered, linewidth=1.5, label="KB integral")
            axislegend(ax_kb; position=:rt, labelsize=10)
            # Update overview text
            try
                overview_obs[] = _overview_text(R)
            catch
                overview_obs[] = "Overview unavailable."
            end
            status_obs[] = "Plots updated"
        end

        # ── Tab 1: MDDF limits ────────────────────────────────────────
        on(btn_mddf_lims.value) do _
            xlims!(ax_mddf, Float64(tf_mddf_xmin.value[]), Float64(tf_mddf_xmax.value[]))
            ylims!(ax_mddf, Float64(tf_mddf_ymin.value[]), Float64(tf_mddf_ymax.value[]))
        end

        # ── Tab 1: KB limits ──────────────────────────────────────────
        on(btn_kb_lims.value) do _
            xlims!(ax_kb, Float64(tf_kb_xmin.value[]), Float64(tf_kb_xmax.value[]))
            ylims!(ax_kb, Float64(tf_kb_ymin.value[]), Float64(tf_kb_ymax.value[]))
        end

        # ── Tab 2/3: group contributions ──────────────────────────────
        palette = [:orangered, :green3, :purple, :goldenrod, :deeppink,
                   :teal, :slateblue, :sienna, :cyan4, :olive]

        function _update_sol_grp_plots!()
            R = result_obs[]; R === nothing && return
            at = atoms_obs[]; at === nothing && return
            comp_sel = String(tf_comp_sol.value[])
            n = grp_active_sol[]
            active_sels = String[]
            for i in 1:n; grp_checks_sol[i].value[] && push!(active_sels, grp_names_sol[i][]); end
            group_labels = String[]; mddf_curves = Vector{Float64}[]; cn_curves = Vector{Float64}[]
            for sel in active_sels
                combined_sel = isempty(strip(comp_sel)) ? sel : "($comp_sel) and ($sel)"
                local sel_atoms
                try; sel_atoms = PDBTools.select(at, combined_sel)
                catch e; status_obs[] = "Error in selection '$combined_sel': $(sprint(showerror, e))"; return; end
                isempty(sel_atoms) && (status_obs[] = "Selection '$combined_sel' matched no atoms"; return)
                grp = SoluteGroup(sel_atoms)
                push!(group_labels, sel)
                push!(mddf_curves, contributions(R, grp; type=:mddf))
                push!(cn_curves, contributions(R, grp; type=:coordination_number))
            end
            empty!(ax_sol_mddf)
            for c in copy(fig_sol_mddf.content); c isa Legend && delete!(c); end
            if grp_total_sol.value[]
                lines!(ax_sol_mddf, R.d, R.mddf; color=:dodgerblue, linewidth=2, label="Total MDDF")
                hlines!(ax_sol_mddf, [1.0]; color=:gray60, linestyle=:dash)
            end
            for (k, lab) in enumerate(group_labels)
                lines!(ax_sol_mddf, R.d, mddf_curves[k]; color=palette[mod1(k, length(palette))], linewidth=1.5, label=lab)
            end
            axislegend(ax_sol_mddf; position=:rt, labelsize=9)
            empty!(ax_sol_cn)
            for c in copy(fig_sol_cn.content); c isa Legend && delete!(c); end
            if grp_total_sol.value[]
                lines!(ax_sol_cn, R.d, R.coordination_number; color=:dodgerblue, linewidth=2, label="Total")
            end
            for (k, lab) in enumerate(group_labels)
                lines!(ax_sol_cn, R.d, cn_curves[k]; color=palette[mod1(k, length(palette))], linewidth=1.5, label=lab)
            end
            axislegend(ax_sol_cn; position=:lt, labelsize=9)
            last_sol_data[] = (
                d            = copy(R.d),
                total_mddf   = grp_total_sol.value[] ? copy(R.mddf) : nothing,
                total_cn     = grp_total_sol.value[] ? copy(R.coordination_number) : nothing,
                group_labels = copy(group_labels),
                mddf_curves  = copy(mddf_curves),
                cn_curves    = copy(cn_curves),
            )
            status_obs[] = "Solute group contributions updated ($(length(group_labels)) groups)"
        end

        function _update_slv_grp_plots!()
            R = result_obs[]; R === nothing && return
            at = atoms_obs[]; at === nothing && return
            comp_sel = String(tf_comp_slv.value[])
            n = grp_active_slv[]
            active_sels = String[]
            for i in 1:n; grp_checks_slv[i].value[] && push!(active_sels, grp_names_slv[i][]); end
            group_labels = String[]; mddf_curves = Vector{Float64}[]; cn_curves = Vector{Float64}[]
            for sel in active_sels
                combined_sel = isempty(strip(comp_sel)) ? sel : "($comp_sel) and ($sel)"
                local sel_atoms
                try; sel_atoms = PDBTools.select(at, combined_sel)
                catch e; status_obs[] = "Error in selection '$combined_sel': $(sprint(showerror, e))"; return; end
                isempty(sel_atoms) && (status_obs[] = "Selection '$combined_sel' matched no atoms"; return)
                grp = SolventGroup(sel_atoms)
                push!(group_labels, sel)
                push!(mddf_curves, contributions(R, grp; type=:mddf))
                push!(cn_curves, contributions(R, grp; type=:coordination_number))
            end
            empty!(ax_slv_mddf)
            for c in copy(fig_slv_mddf.content); c isa Legend && delete!(c); end
            if grp_total_slv.value[]
                lines!(ax_slv_mddf, R.d, R.mddf; color=:dodgerblue, linewidth=2, label="Total MDDF")
                hlines!(ax_slv_mddf, [1.0]; color=:gray60, linestyle=:dash)
            end
            for (k, lab) in enumerate(group_labels)
                lines!(ax_slv_mddf, R.d, mddf_curves[k]; color=palette[mod1(k, length(palette))], linewidth=1.5, label=lab)
            end
            axislegend(ax_slv_mddf; position=:rt, labelsize=9)
            empty!(ax_slv_cn)
            for c in copy(fig_slv_cn.content); c isa Legend && delete!(c); end
            if grp_total_slv.value[]
                lines!(ax_slv_cn, R.d, R.coordination_number; color=:dodgerblue, linewidth=2, label="Total")
            end
            for (k, lab) in enumerate(group_labels)
                lines!(ax_slv_cn, R.d, cn_curves[k]; color=palette[mod1(k, length(palette))], linewidth=1.5, label=lab)
            end
            axislegend(ax_slv_cn; position=:lt, labelsize=9)
            last_slv_data[] = (
                d            = copy(R.d),
                total_mddf   = grp_total_slv.value[] ? copy(R.mddf) : nothing,
                total_cn     = grp_total_slv.value[] ? copy(R.coordination_number) : nothing,
                group_labels = copy(group_labels),
                mddf_curves  = copy(mddf_curves),
                cn_curves    = copy(cn_curves),
            )
            status_obs[] = "Solvent group contributions updated ($(length(group_labels)) groups)"
        end

        on(grp_total_sol.value) do _; _update_sol_grp_plots!(); end
        on(grp_total_slv.value) do _; _update_slv_grp_plots!(); end
        on(grp_active_sol) do n; n > 0 && _update_sol_grp_plots!(); end
        on(grp_active_slv) do n; n > 0 && _update_slv_grp_plots!(); end
        for i in 1:MAX_GROUPS
            on(grp_checks_sol[i].value) do _; i <= grp_active_sol[] && _update_sol_grp_plots!(); end
            on(grp_checks_slv[i].value) do _; i <= grp_active_slv[] && _update_slv_grp_plots!(); end
        end
        on(tf_comp_sol.value) do _; grp_active_sol[] > 0 && _update_sol_grp_plots!(); end
        on(tf_comp_slv.value) do _; grp_active_slv[] > 0 && _update_slv_grp_plots!(); end

        # ── Tab 2/3: export helpers ───────────────────────────────────
        function _export_fig(fig, tf_name, dd_fmt)
            fmt  = String(dd_fmt.value[])
            base = strip(String(tf_name.value[]))
            isempty(base) && (base = "export")
            path = endswith(base, ".$fmt") ? base : "$base.$fmt"
            try
                WGLMakie.save(path, fig)
                status_obs[] = "Saved: $path"
            catch e
                status_obs[] = "Export error: $(sprint(showerror, e))"
            end
        end

        function _export_csv(data, tf_name)
            data === nothing && (status_obs[] = "No data to export — run a plot first"; return)
            base = strip(String(tf_name.value[]))
            isempty(base) && (base = "export")
            path = endswith(base, ".csv") ? base : "$base.csv"
            try
                open(path, "w") do io
                    # Header
                    cols = ["d"]
                    data.total_mddf !== nothing && push!(cols, "Total_MDDF")
                    for lab in data.group_labels; push!(cols, "$(lab)_mddf"); end
                    data.total_cn !== nothing && push!(cols, "Total_CN")
                    for lab in data.group_labels; push!(cols, "$(lab)_cn"); end
                    println(io, join(cols, ","))
                    # Rows
                    for j in eachindex(data.d)
                        vals = [string(data.d[j])]
                        data.total_mddf !== nothing && push!(vals, string(data.total_mddf[j]))
                        for k in eachindex(data.group_labels); push!(vals, string(data.mddf_curves[k][j])); end
                        data.total_cn !== nothing && push!(vals, string(data.total_cn[j]))
                        for k in eachindex(data.group_labels); push!(vals, string(data.cn_curves[k][j])); end
                        println(io, join(vals, ","))
                    end
                end
                status_obs[] = "Saved: $path"
            catch e
                status_obs[] = "Export error: $(sprint(showerror, e))"
            end
        end

        on(btn_export_mddf_sol.value) do _; _export_fig(fig_sol_mddf, tf_export_sol, dd_export_fmt_sol); end
        on(btn_export_cn_sol.value)   do _; _export_fig(fig_sol_cn,   tf_export_sol, dd_export_fmt_sol); end
        on(btn_export_csv_sol.value)  do _; _export_csv(last_sol_data[], tf_export_sol); end
        on(btn_export_mddf_slv.value) do _; _export_fig(fig_slv_mddf, tf_export_slv, dd_export_fmt_slv); end
        on(btn_export_cn_slv.value)   do _; _export_fig(fig_slv_cn,   tf_export_slv, dd_export_fmt_slv); end
        on(btn_export_csv_slv.value)  do _; _export_csv(last_slv_data[], tf_export_slv); end

        # ── Tab 2/3: apply group limits ───────────────────────────────
        on(btn_grp_lims_sol.value) do _
            xlo = Float64(tf_grp_xmin.value[]); xhi = Float64(tf_grp_xmax.value[])
            ylo = Float64(tf_grp_ymin.value[]); yhi = Float64(tf_grp_ymax.value[])
            xlims!(ax_sol_mddf, xlo, xhi); xlims!(ax_sol_cn, xlo, xhi)
            ylims!(ax_sol_mddf, ylo, yhi); ylims!(ax_sol_cn, 0, nothing)
        end
        on(btn_grp_lims_slv.value) do _
            xlo = Float64(tf_grp_xmin.value[]); xhi = Float64(tf_grp_xmax.value[])
            ylo = Float64(tf_grp_ymin.value[]); yhi = Float64(tf_grp_ymax.value[])
            xlims!(ax_slv_mddf, xlo, xhi); xlims!(ax_slv_cn, xlo, xhi)
            ylims!(ax_slv_mddf, ylo, yhi); ylims!(ax_slv_cn, 0, nothing)
        end

        # ── Tab 3: residue contributions ──────────────────────────────
        on(btn_rc_update.value) do _
            R = result_obs[]
            R === nothing && (status_obs[] = "No results loaded"; return)
            at = atoms_obs[]
            at === nothing && (status_obs[] = "Load a PDB file first"; return)

            sel_str = String(tf_rc_sel.value[])
            dmin = Float64(tf_dmin.value[])
            dmax = Float64(tf_dmax.value[])

            status_obs[] = "Computing residue contributions…"
            try
                sel_atoms = PDBTools.select(at, sel_str)
                if isempty(sel_atoms)
                    status_obs[] = "Selection '$sel_str' matched no atoms"; return
                end
                rc_mddf = ResidueContributions(R, sel_atoms; dmin=dmin, dmax=dmax, type=:mddf, silent=true)
                rc_cn   = ResidueContributions(R, sel_atoms; dmin=dmin, dmax=dmax, type=:coordination_number, silent=true)

                nres = length(rc_mddf.resnums)
                x_pos = rc_mddf.xticks[1]
                y_d   = rc_mddf.d

                rc_range = 1:nres
                if nres > 2000
                    rc_step = nres ÷ 2000
                    rc_range = 1:rc_step:nres
                end
                x_plot = x_pos[rc_range]

                step = max(1, length(rc_range) ÷ 50)
                tick_idx = 1:step:length(x_plot)
                orig_idx = collect(rc_range)[collect(tick_idx)]
                xtick_pos = rc_mddf.xticks[1][orig_idx]
                xtick_lab = rc_mddf.xticks[2][orig_idx]

                # MDDF plot
                zmat_mddf = collect(hcat(rc_mddf.residue_contributions...)')[rc_range, :]
                clims_mddf, cscale_mddf = _set_clims_and_colorscale!(rc_mddf)
                cmap_mddf = cscale_mddf == :bwr ? :RdBu : :tempo
                nlevels_mddf = cscale_mddf == :tempo ? 5 : 12
                empty!(ax_rc_mddf)
                contourf!(ax_rc_mddf, x_plot, y_d, zmat_mddf; colormap=cmap_mddf, levels=nlevels_mddf)
                ax_rc_mddf.xticks = (xtick_pos, xtick_lab)
                ax_rc_mddf.xticklabelrotation = π / 3
                if length(fig3_mddf.layout.content) > 1
                    try delete!(fig3_mddf.layout.content[end].content) catch end
                end
                Colorbar(fig3_mddf[1, 2]; colormap=cmap_mddf, limits=clims_mddf,
                    label="Contribution", labelsize=10, ticklabelsize=9)

                # Coordination number plot
                zmat_cn = collect(hcat(rc_cn.residue_contributions...)')[rc_range, :]
                clims_cn, cscale_cn = _set_clims_and_colorscale!(rc_cn)
                cmap_cn = cscale_cn == :bwr ? :RdBu : :tempo
                nlevels_cn = cscale_cn == :tempo ? 5 : 12
                empty!(ax_rc_cn)
                contourf!(ax_rc_cn, x_plot, y_d, zmat_cn; colormap=cmap_cn, levels=nlevels_cn)
                ax_rc_cn.xticks = (xtick_pos, xtick_lab)
                ax_rc_cn.xticklabelrotation = π / 3
                if length(fig3_cn.layout.content) > 1
                    try delete!(fig3_cn.layout.content[end].content) catch end
                end
                Colorbar(fig3_cn[1, 2]; colormap=cmap_cn, limits=clims_cn,
                    label="Contribution", labelsize=10, ticklabelsize=9)

                status_obs[] = "Residue contributions updated ($nres residues)"
            catch e
                status_obs[] = "Error: $(sprint(showerror, e))"
            end
        end

        # Trigger initial plot if data was preloaded
        if !isnothing(result_obs[])
            notify(result_obs)
        end

        # ── Assemble page ─────────────────────────────────────────────
        _version = string(pkgversion(ComplexMixtures))
        return DOM.div(
            _CSS,
            DOM.div(class="cm-title", "ComplexMixtures v$(_version)"),
            DOM.div(class="cm-main", sidebar, plots_panel),
            status_bar,
        )
    end

    server = Bonito.Server(app, "0.0.0.0", port)
    url = "http://localhost:$port"
    @info "ComplexMixtures GUI running at $url"
    try
        Bonito.open_browser(url)
    catch
        @info "Open $url in your browser"
    end
    return server
end

end # module GUI
