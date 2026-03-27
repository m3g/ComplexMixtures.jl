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
.cm-tab2-body { display: flex; gap: 10px; align-items: flex-start; }
.cm-grp-panel { width: 200px; flex-shrink: 0; border: 1px solid #ddd; border-radius: 4px; padding: 8px; background: #fafafa; }
.cm-grp-panel-title { font-weight: bold; font-size: 12px; color: #333; margin-bottom: 6px; }
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
        fig1 = Figure(; size=(900, 620))
        ax_mddf = Axis(fig1[1, 1]; ylabel="g(r)", ylabelsize=14, xticklabelsize=11, yticklabelsize=11)
        Label(fig1[1, 1, Top()], "MDDF"; fontsize=16, padding=(0, 0, 4, 0))
        ax_kb = Axis(fig1[2, 1]; ylabel="KB (cm3/mol)", ylabelsize=14, xticklabelsize=11, yticklabelsize=11)
        Label(fig1[2, 1, Top()], "Kirkwood-Buff Integral"; fontsize=16, padding=(0, 0, 4, 0))

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

        # ── Tab 2: Group Contributions ────────────────────────────────
        fig2 = Figure(; size=(900, 620))
        ax_grp_mddf = Axis(fig2[1, 1]; ylabel="g(r)", ylabelsize=14, xticklabelsize=11, yticklabelsize=11)
        Label(fig2[1, 1, Top()], "Group MDDF Contributions"; fontsize=16, padding=(0, 0, 4, 0))
        ax_grp_cn = Axis(fig2[2, 1]; ylabel="Coordination number", ylabelsize=14, xticklabelsize=11, yticklabelsize=11)
        Label(fig2[2, 1, Top()], "Group Coord. Number Contributions"; fontsize=16, padding=(0, 0, 4, 0))

        MAX_GROUPS = 10
        tf_newgrp = Bonito.TextField("resname ARG")
        btn_addgrp = Bonito.Button("Add")
        grp_names = [Observable("") for _ in 1:MAX_GROUPS]
        grp_checks = [Bonito.Checkbox(true) for _ in 1:MAX_GROUPS]
        grp_active = Observable(0)
        grp_total_check = Bonito.Checkbox(true)
        btn_rmgrp = Bonito.Button("Remove unchecked")
        dd_side = Bonito.Dropdown(["Solute", "Solvent"])

        checklist_dom = map(grp_active) do n
            rows = Any[DOM.div(class="cm-grp-item",
                grp_total_check,
                DOM.span("Total MDDF"),
            )]
            for i in 1:n
                push!(rows, DOM.div(class="cm-grp-item",
                    grp_checks[i],
                    DOM.span(grp_names[i][]),
                ))
            end
            DOM.div(class="cm-grp-list", rows...)
        end

        # Group limits (below plots)
        tf_grp_xmin = Bonito.NumberInput(0.0)
        tf_grp_xmax = Bonito.NumberInput(10.0)
        tf_grp_ymin = Bonito.NumberInput(0.0)
        tf_grp_ymax = Bonito.NumberInput(5.0)
        btn_grp_lims = Bonito.Button("Apply")

        tab2_body = DOM.div(class="cm-tab2-body",
            DOM.div(
                fig2,
                DOM.div(style="text-align:center; font-size:13px; color:#444; margin: 0 0 2px 0;", "r (Angstrom)"),
            ),
            DOM.div(class="cm-grp-panel",
                DOM.div(class="cm-grp-panel-title", "Groups"),
                DOM.div(class="cm-lims-row", tf_newgrp, btn_addgrp),
                DOM.div(class="cm-row cm-row-left", DOM.label("Side:"), dd_side),
                checklist_dom,
                DOM.div(style="margin: 4px 0;", btn_rmgrp),
            ),
        )

        tab2_limits = DOM.div(class="cm-lims-row",
            DOM.span("Limits:"),
            DOM.span("x:"), tf_grp_xmin, DOM.span("–"), tf_grp_xmax,
            DOM.span("y:"), tf_grp_ymin, DOM.span("–"), tf_grp_ymax,
            btn_grp_lims,
        )

        # ── Tab 3: Residue Contributions ──────────────────────────────
        fig3 = Figure(; size=(900, 500))
        ax_rc = Axis(fig3[1, 1]; ylabel="r (Angstrom)",
            ylabelsize=14, xticklabelsize=9, yticklabelsize=11)
        tf_rc_sel = Bonito.TextField("protein")
        tf_dmin = Bonito.NumberInput(1.5)
        tf_dmax = Bonito.NumberInput(3.5)
        btn_rc_update = Bonito.Button("Update")

        tab3_controls = DOM.div(
            DOM.div(class="cm-row", DOM.label("Selection:"), tf_rc_sel),
            DOM.div(class="cm-lims-row",
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
                DOM.div(class="cm-tab-btn", "Group Contributions"),
                DOM.div(class="cm-tab-btn", "Residue Contributions"),
            ),
            DOM.div(class="cm-tab-content active",
                fig1,
                DOM.div(style="text-align:center; font-size:13px; color:#444; margin: 0 0 2px 0;", "r (Angstrom)"),
                tab1_limits,
            ),
            DOM.div(class="cm-tab-content", tab2_body, tab2_limits),
            DOM.div(class="cm-tab-content", tab3_controls,
                DOM.div(style="text-align:center; font-size:15px; font-weight:bold; color:#222; margin:4px 0 0 0;", "Residue Contributions"),
                fig3,
                DOM.div(style="text-align:center; font-size:13px; color:#444; margin-top:0;", "Residue"),
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
        on(btn_addgrp.value) do _
            sel = strip(String(tf_newgrp.value[]))
            isempty(sel) && return
            n = grp_active[]
            n >= MAX_GROUPS && (status_obs[] = "Maximum $MAX_GROUPS groups reached"; return)
            for i in 1:n
                grp_names[i][] == sel && return  # duplicate
            end
            grp_names[n + 1][] = sel
            grp_checks[n + 1].value[] = true
            grp_active[] = n + 1
        end
        on(btn_rmgrp.value) do _
            n = grp_active[]
            kept = String[]
            for i in 1:n
                grp_checks[i].value[] && push!(kept, grp_names[i][])
            end
            for i in eachindex(kept)
                grp_names[i][] = kept[i]
                grp_checks[i].value[] = true
            end
            grp_active[] = length(kept)
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
            ax_mddf.ylabel = "g(r)"
            empty!(ax_kb)
            ax_kb.ylabel = "KB (cm3/mol)"
            # Remove old legends from fig1
            for c in copy(fig1.content)
                c isa Legend && delete!(c)
            end
            lines!(ax_mddf, R.d, R.mddf; color=:dodgerblue, linewidth=1.5, label="MDDF")
            hlines!(ax_mddf, [1.0]; color=:gray60, linestyle=:dash)
            axislegend(ax_mddf; position=:rt, labelsize=10)
            lines!(ax_kb, R.d, R.kb; color=:orangered, linewidth=1.5, label="KB integral")
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

        # ── Tab 2: group contributions ────────────────────────────────
        function _update_grp_plots!()
            R = result_obs[]
            R === nothing && return
            at = atoms_obs[]
            at === nothing && return

            side = String(dd_side.value[])
            GrpT = side == "Solute" ? SoluteGroup : SolventGroup

            n = grp_active[]
            active_sels = String[]
            for i in 1:n
                grp_checks[i].value[] && push!(active_sels, grp_names[i][])
            end

            group_labels = String[]
            mddf_curves = Vector{Float64}[]
            cn_curves = Vector{Float64}[]
            for sel in active_sels
                local sel_atoms
                try
                    sel_atoms = PDBTools.select(at, sel)
                catch e
                    status_obs[] = "Error in selection '$sel': $(sprint(showerror, e))"; return
                end
                if isempty(sel_atoms)
                    status_obs[] = "Selection '$sel' matched no atoms"; return
                end
                grp = GrpT(sel_atoms)
                push!(group_labels, sel)
                push!(mddf_curves, contributions(R, grp; type=:mddf))
                push!(cn_curves, contributions(R, grp; type=:coordination_number))
            end

            palette = [:orangered, :green3, :purple, :goldenrod, :deeppink,
                       :teal, :slateblue, :sienna, :cyan4, :olive]
            empty!(ax_grp_mddf)
            ax_grp_mddf.ylabel = "g(r)"
            for c in copy(fig2.content)
                c isa Legend && delete!(c)
            end
            if grp_total_check.value[]
                lines!(ax_grp_mddf, R.d, R.mddf; color=:dodgerblue, linewidth=2, label="Total MDDF")
                hlines!(ax_grp_mddf, [1.0]; color=:gray60, linestyle=:dash)
            end
            for (k, lab) in enumerate(group_labels)
                c = palette[mod1(k, length(palette))]
                lines!(ax_grp_mddf, R.d, mddf_curves[k]; color=c, linewidth=1.5, label=lab)
            end
            axislegend(ax_grp_mddf; position=:rt, labelsize=9)
            empty!(ax_grp_cn)
            ax_grp_cn.ylabel = "Coordination number"
            if grp_total_check.value[]
                lines!(ax_grp_cn, R.d, R.coordination_number; color=:dodgerblue, linewidth=2, label="Total")
            end
            for (k, lab) in enumerate(group_labels)
                c = palette[mod1(k, length(palette))]
                lines!(ax_grp_cn, R.d, cn_curves[k]; color=c, linewidth=1.5, label=lab)
            end
            axislegend(ax_grp_cn; position=:lt, labelsize=9)
            status_obs[] = "Group contributions updated ($(length(group_labels)) checked, $side)"
        end

        on(grp_total_check.value) do _
            _update_grp_plots!()
        end
        on(dd_side.value) do _
            _update_grp_plots!()
        end
        on(grp_active) do n
            n > 0 && _update_grp_plots!()
        end
        for i in 1:MAX_GROUPS
            on(grp_checks[i].value) do _
                i <= grp_active[] && _update_grp_plots!()
            end
        end

        # ── Tab 2: apply group limits ─────────────────────────────────
        on(btn_grp_lims.value) do _
            xlo = Float64(tf_grp_xmin.value[])
            xhi = Float64(tf_grp_xmax.value[])
            ylo = Float64(tf_grp_ymin.value[])
            yhi = Float64(tf_grp_ymax.value[])
            xlims!(ax_grp_mddf, xlo, xhi)
            xlims!(ax_grp_cn, xlo, xhi)
            ylims!(ax_grp_mddf, ylo, yhi)
            ylims!(ax_grp_cn, 0, nothing)
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
                rc = ResidueContributions(R, sel_atoms; dmin=dmin, dmax=dmax, silent=true)

                nres = length(rc.resnums)
                x_pos = rc.xticks[1]
                y_d = rc.d
                zmat = collect(hcat(rc.residue_contributions...)')  # (nres × nd)

                rc_range = 1:nres
                if nres > 2000
                    rc_step = nres ÷ 2000
                    rc_range = 1:rc_step:nres
                end
                x_plot = x_pos[rc_range]
                z_plot = zmat[rc_range, :]

                clims, cscale = _set_clims_and_colorscale!(rc)
                cmap = cscale == :bwr ? :RdBu : :tempo
                nlevels = cscale == :tempo ? 5 : 12

                step = max(1, length(rc_range) ÷ 50)
                tick_idx = 1:step:length(x_plot)
                orig_idx = collect(rc_range)[collect(tick_idx)]
                xtick_pos = rc.xticks[1][orig_idx]
                xtick_lab = rc.xticks[2][orig_idx]

                empty!(ax_rc)
                ax_rc.ylabel = "r (Angstrom)"
                contourf!(ax_rc, x_plot, y_d, z_plot; colormap=cmap, levels=nlevels)
                ax_rc.xticks = (xtick_pos, xtick_lab)
                ax_rc.xticklabelrotation = π / 3

                # Add colorbar (remove previous if exists)
                if length(fig3.layout.content) > 1
                    try delete!(fig3.layout.content[end].content) catch end
                end
                Colorbar(fig3[1, 2]; colormap=cmap, limits=clims,
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
