import csv
import html
import os
import shutil
import zipfile
from pathlib import Path

OUT_DIR = Path("summary_reports")
OUT_DIR.mkdir(exist_ok=True)
PPTX_PATH = OUT_DIR / "meta_analysis_results_presentation.pptx"

SLIDE_W = 13_333_500
SLIDE_H = 7_500_000
EMU_PER_IN = 914400

BG = "F8FAFC"
INK = "172033"
MUTED = "64748B"
BLUE = "2563EB"
TEAL = "0F766E"
GREEN = "16A34A"
RED = "DC2626"
AMBER = "D97706"
LIGHT_BLUE = "DBEAFE"
LIGHT_TEAL = "CCFBF1"
LIGHT_GREEN = "DCFCE7"
LIGHT_RED = "FEE2E2"
LIGHT_AMBER = "FEF3C7"
WHITE = "FFFFFF"


def emu(inches):
    return int(inches * EMU_PER_IN)


def esc(text):
    return html.escape(str(text), quote=False)


def read_csv(path):
    if not Path(path).exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


filter_counts = read_csv("meta_filter_sensitivity_output/filter_isolation_counts.csv")
strict_genes = read_csv("meta_filter_sensitivity_output/strict_high_confidence_adj0.1_i2lt50_abs_effect_gt0.25.csv")
relaxed_summary = read_csv("meta_filter_sensitivity_output/strict_vs_relaxed_filter_summary.csv")
loo_summary = read_csv("leave_one_study_out_meta_output/leave_one_study_out_summary.csv")
gprofiler = read_csv("gprofiler_enrichment_output/gprofiler_ora_results.csv")
gprofiler_relaxed = read_csv("gprofiler_enrichment_output/gprofiler_ora_results_relaxed_no_i2_abs0.1.csv")


def metric(summary, name, default="NA"):
    for row in summary:
        if row.get("metric") == name:
            return row.get("value", default)
    return default


def filter_count(name):
    for row in filter_counts:
        if row.get("filter") == name:
            return row.get("count", "NA")
    return "NA"


def fmt_float(value, digits=3):
    try:
        x = float(value)
    except Exception:
        return "NA"
    if abs(x) < 0.001:
        return f"{x:.2e}"
    return f"{x:.{digits}g}"


def top_terms(rows, query_name=None, n=5):
    filtered = rows
    if query_name is not None:
        filtered = [r for r in rows if r.get("query_name") == query_name]
    return sorted(filtered, key=lambda x: float(x.get("p_value", "inf")))[:n]


def p_xml(text, size=24, color=INK, bold=False, align="l"):
    b = '<a:b/>' if bold else ''
    return (
        f'<a:p><a:pPr algn="{align}"/>'
        f'<a:r><a:rPr lang="en-US" sz="{int(size*100)}" dirty="0">{b}'
        f'<a:solidFill><a:srgbClr val="{color}"/></a:solidFill></a:rPr>'
        f'<a:t>{esc(text)}</a:t></a:r></a:p>'
    )


def text_box(shape_id, name, x, y, w, h, paragraphs, fill=None, line=None, radius=False):
    fill_xml = '<a:noFill/>' if fill is None else f'<a:solidFill><a:srgbClr val="{fill}"/></a:solidFill>'
    line_xml = '<a:ln><a:noFill/></a:ln>' if line is None else f'<a:ln w="12700"><a:solidFill><a:srgbClr val="{line}"/></a:solidFill></a:ln>'
    prst = "roundRect" if radius else "rect"
    body = ''.join(paragraphs)
    return f"""
      <p:sp>
        <p:nvSpPr><p:cNvPr id="{shape_id}" name="{esc(name)}"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>
        <p:spPr>
          <a:xfrm><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm>
          <a:prstGeom prst="{prst}"><a:avLst/></a:prstGeom>
          {fill_xml}
          {line_xml}
        </p:spPr>
        <p:txBody>
          <a:bodyPr wrap="square" anchor="mid"><a:spAutoFit/></a:bodyPr>
          <a:lstStyle/>
          {body}
        </p:txBody>
      </p:sp>
    """


def rect(shape_id, name, x, y, w, h, fill, line=None, radius=False):
    return text_box(shape_id, name, x, y, w, h, [], fill=fill, line=line, radius=radius)


def image_xml(shape_id, name, r_id, x, y, w, h):
    return f"""
      <p:pic>
        <p:nvPicPr><p:cNvPr id="{shape_id}" name="{esc(name)}"/><p:cNvPicPr/><p:nvPr/></p:nvPicPr>
        <p:blipFill>
          <a:blip r:embed="{r_id}"/>
          <a:stretch><a:fillRect/></a:stretch>
        </p:blipFill>
        <p:spPr>
          <a:xfrm><a:off x="{x}" y="{y}"/><a:ext cx="{w}" cy="{h}"/></a:xfrm>
          <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
        </p:spPr>
      </p:pic>
    """


def slide_xml(elements, bg=BG):
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
       xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
       xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
  <p:cSld>
    <p:bg><p:bgPr><a:solidFill><a:srgbClr val="{bg}"/></a:solidFill></p:bgPr></p:bg>
    <p:spTree>
      <p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>
      <p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>
      {''.join(elements)}
    </p:spTree>
  </p:cSld>
  <p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr>
</p:sld>"""


def title_elements(title, subtitle=None, kicker=None):
    els = []
    if kicker:
        els.append(text_box(10, "kicker", emu(0.7), emu(0.45), emu(6.8), emu(0.35), [p_xml(kicker.upper(), 12, BLUE, True)]))
    els.append(text_box(11, "title", emu(0.7), emu(0.78), emu(11.4), emu(0.9), [p_xml(title, 31, INK, True)]))
    if subtitle:
        els.append(text_box(12, "subtitle", emu(0.72), emu(1.48), emu(10.8), emu(0.55), [p_xml(subtitle, 15, MUTED)]))
    return els


slides = []
image_assets = []


def add_slide(elements, images=None, bg=BG):
    slides.append({"xml": slide_xml(elements, bg), "images": images or []})


# 1. Cover
cover = [
    rect(2, "left-band", 0, 0, emu(4.15), SLIDE_H, INK),
    rect(3, "accent", emu(0.55), emu(5.72), emu(2.15), emu(0.08), BLUE),
    text_box(4, "cover-kicker", emu(0.65), emu(0.85), emu(3.2), emu(0.35), [p_xml("ADIPOSE TRANSCRIPTOMICS", 12, "93C5FD", True)]),
    text_box(5, "cover-title", emu(0.65), emu(1.35), emu(3.2), emu(2.2), [p_xml("Meta-analysis results", 34, WHITE, True)]),
    text_box(6, "cover-subtitle", emu(0.68), emu(3.55), emu(2.9), emu(0.9), [p_xml("Exercise and lifestyle intervention signals across GEO cohorts", 15, "CBD5E1")]),
    text_box(7, "cover-date", emu(0.68), emu(6.0), emu(3.2), emu(0.35), [p_xml("Generated from current pipeline outputs", 11, "94A3B8")]),
]
for i, (label, value, color) in enumerate([
    ("studies", "5", BLUE),
    ("meta-analyzed genes", "29,193", TEAL),
    ("strict genes", "12", AMBER),
    ("relaxed genes", "105", GREEN),
]):
    x = emu(4.75 + (i % 2) * 3.85)
    y = emu(1.25 + (i // 2) * 2.35)
    cover.append(text_box(20 + i, f"metric-{label}", x, y, emu(3.0), emu(1.35), [
        p_xml(value, 39, color, True),
        p_xml(label, 13, MUTED)
    ]))
cover.append(text_box(30, "cover-thesis", emu(4.75), emu(5.35), emu(7.55), emu(0.8), [
    p_xml("A small conservative gene set emerges from broad adipose intervention data, while sensitivity analyses show where the signal is fragile.", 17, INK, True)
]))
add_slide(cover)


# 2. Study question and data shape
els = title_elements("What this analysis is asking", "A meta-analysis of adipose tissue expression changes across public intervention datasets.", "Framing")
for i, (label, text, color) in enumerate([
    ("Unit of analysis", "Each GEO study is analyzed separately with limma, then combined gene by gene.", BLUE),
    ("Effect exported", "Study-level limma logFC; variance estimated from logFC / t.", TEAL),
    ("Model", "Random-effects meta-analysis, REML estimator, DL fallback only if needed.", AMBER),
]):
    y = emu(2.3 + i * 1.28)
    els.append(rect(40 + i, f"rule-{i}", emu(0.85), y + emu(0.12), emu(0.09), emu(0.72), color))
    els.append(text_box(50 + i, f"point-{i}", emu(1.05), y, emu(10.8), emu(0.9), [
        p_xml(label, 18, INK, True),
        p_xml(text, 15, MUTED)
    ]))
els.append(text_box(60, "bottom-line", emu(0.85), emu(6.25), emu(11.6), emu(0.55), [
    p_xml("Best wording: exercise/lifestyle-associated adipose transcriptional response, not a pure endurance-training signature.", 17, RED, True)
]))
add_slide(els)


# 3. Random-effects model logic
els = title_elements("How REMA combines studies", "For every gene, each study contributes an effect estimate and its uncertainty.", "Model logic")
model_steps = [
    ("Study effect", "yi = study-level limma logFC", BLUE),
    ("Study variance", "vi = estimated variance from limma logFC and t", TEAL),
    ("Random-effects weight", "wi* = 1 / (vi + tau2)", AMBER),
    ("Pooled effect", "mu-hat = weighted mean of yi using wi*", GREEN),
]
for i, (head, body, color) in enumerate(model_steps):
    x = emu(0.85 + (i % 2) * 6.15)
    y = emu(2.05 + (i // 2) * 1.65)
    els.append(rect(65 + i, f"model-rule-{i}", x, y + emu(0.1), emu(0.09), emu(0.9), color))
    els.append(text_box(70 + i, f"model-step-{i}", x + emu(0.22), y, emu(5.25), emu(1.05), [
        p_xml(head, 17, INK, True),
        p_xml(body, 15, MUTED)
    ]))
els.append(text_box(80, "model-note", emu(1.15), emu(5.65), emu(11.0), emu(0.85), [
    p_xml("tau2 is between-study heterogeneity. We estimate tau2 with REML, then report the pooled logFC, confidence interval, heterogeneity, and prediction interval.", 17, INK, True, "ctr")
]))
add_slide(els)


# 4. Pipeline figure
img = Path("figures/prisma_like_analysis_flow.png")
els = title_elements("From GEO data to meta-analysis", "The same pipeline turns study-level limma outputs into pooled gene-level evidence.", "Workflow")
if img.exists():
    els.append(image_xml(70, "pipeline-flow", "rId2", emu(0.95), emu(1.9), emu(11.4), emu(4.95)))
    add_slide(els, images=[img])
else:
    els.append(text_box(70, "missing", emu(1), emu(2), emu(10), emu(1), [p_xml("Pipeline figure not found.", 22, RED, True)]))
    add_slide(els)


# 5. Filter sensitivity
els = title_elements("Filter sensitivity: what changes?", "The strict 12 genes mainly come from the effect-size threshold, not the I2 filter.", "Results")
steps = [
    ("All meta-analyzed", filter_count("all_meta_analyzed_genes"), LIGHT_BLUE, BLUE),
    ("FDR < 0.1", filter_count("FDR_lt_0.1"), LIGHT_TEAL, TEAL),
    ("FDR + I2 < 50", filter_count("FDR_lt_0.1_and_I2_lt_50"), LIGHT_GREEN, GREEN),
    ("FDR + |beta| > 0.1", filter_count("FDR_lt_0.1_and_abs_effect_gt_0.1"), LIGHT_AMBER, AMBER),
    ("Strict: FDR + I2 + |beta| > 0.25", filter_count("FDR_lt_0.1_and_I2_lt_50_and_abs_effect_gt_0.25"), LIGHT_RED, RED),
]
for i, (label, value, fill, color) in enumerate(steps):
    x = emu(0.75 + i * 2.52)
    els.append(text_box(80 + i, f"filter-{i}", x, emu(2.25), emu(2.15), emu(1.75), [
        p_xml(value, 28, color, True, "ctr"),
        p_xml(label, 11, INK, True, "ctr")
    ], fill=fill, line=color, radius=True))
els.append(text_box(90, "sensitivity-takeaway", emu(1.25), emu(5.0), emu(10.6), emu(0.9), [
    p_xml("Relaxed no-I2 sensitivity gives 105 genes. The strict set is intentionally conservative: 5 upregulated, 7 downregulated.", 20, INK, True, "ctr")
]))
add_slide(els)


# 6. Final 12 genes table
els = title_elements("The strict high-confidence 12", "Genes passing FDR < 0.1, I2 < 50%, and absolute pooled beta > 0.25.", "Core result")
headers = ["Gene", "Direction", "beta", "FDR", "I2", "Studies"]
xs = [0.8, 2.45, 4.05, 5.35, 6.75, 8.0]
ws = [1.4, 1.25, 1.0, 1.1, 0.85, 4.4]
y0 = 1.95
for j, h in enumerate(headers):
    els.append(text_box(100 + j, f"header-{h}", emu(xs[j]), emu(y0), emu(ws[j]), emu(0.32), [p_xml(h, 10.5, WHITE, True)], fill=INK))
for i, row in enumerate(strict_genes[:12]):
    y = y0 + 0.38 + i * 0.36
    direction = "up" if float(row["meta_logFC"]) > 0 else "down"
    dir_color = GREEN if direction == "up" else RED
    vals = [
        row["gene"],
        direction,
        fmt_float(row["meta_logFC"], 2),
        fmt_float(row["adj_p_value"], 2),
        fmt_float(row["I2"], 2),
        row.get("studies", ""),
    ]
    for j, val in enumerate(vals):
        color = dir_color if j == 1 else INK
        els.append(text_box(120 + i * 10 + j, f"gene-{i}-{j}", emu(xs[j]), emu(y), emu(ws[j]), emu(0.29), [p_xml(val, 9.2, color, j in [0, 1])], fill="FFFFFF" if i % 2 == 0 else "F1F5F9"))
add_slide(els)


# 7. Forest plot example
img = Path("meta_filter_sensitivity_output/forest_plots_strict_12/APOE_forest_plot.png")
els = title_elements("Forest plot example: APOE", "Study-level logFC estimates point in the same positive direction across all five study contexts.", "Evidence")
if img.exists():
    els.append(image_xml(250, "apoe-forest", "rId2", emu(0.72), emu(1.8), emu(7.65), emu(4.9)))
els.append(text_box(251, "apoe-note", emu(8.65), emu(2.05), emu(3.55), emu(2.4), [
    p_xml("What to read", 18, INK, True),
    p_xml("Each point is a study-level limma effect.", 13.5, MUTED),
    p_xml("The pooled diamond is the random-effects estimate.", 13.5, MUTED),
    p_xml("Prediction intervals are now exported in the meta-analysis CSV.", 13.5, MUTED),
    p_xml("APOE links this signal to lipid handling in adipose tissue.", 13.5, TEAL, True)
]))
add_slide(els, images=[img] if img.exists() else [])


# 8. Leave-one-out
img = Path("figures/leave_one_study_out_sensitivity.png")
els = title_elements("Sensitivity: leave one study out", "A robustness check repeats the meta-analysis after dropping each GEO dataset.", "Robustness")
if img.exists():
    els.append(image_xml(270, "loo", "rId2", emu(0.72), emu(1.8), emu(7.5), emu(4.85)))
els.append(text_box(271, "loo-note", emu(8.55), emu(1.95), emu(3.7), emu(3.05), [
    p_xml("Key read", 18, INK, True),
    p_xml("No gene survived all six runs.", 14, RED, True),
    p_xml("TMEM170B was strongest by survival count.", 14, MUTED),
    p_xml("GSE58559 strongly influences the broader FDR-only signal.", 14, MUTED)
]))
add_slide(els, images=[img] if img.exists() else [])


# 9. Up/down network and enrichment
img = Path("figures/consistent_up_down_gene_network.png")
els = title_elements("Biological interpretation is suggestive, not definitive", "The strict gene list is small, so pathway enrichment has limited power.", "Biology")
if img.exists():
    els.append(image_xml(290, "network", "rId2", emu(0.65), emu(1.75), emu(6.9), emu(4.9)))
term = "No enriched ORA terms found"
if gprofiler:
    top = gprofiler[0]
    term = f"{top.get('term_name', 'term')} ({top.get('source', '')}); p={fmt_float(top.get('p_value', 'NA'), 2)}; genes={top.get('intersection', '')}"
els.append(text_box(291, "bio-note", emu(7.85), emu(2.05), emu(4.15), emu(2.7), [
    p_xml("gProfiler ORA", 18, INK, True),
    p_xml(term, 13.5, MUTED),
    p_xml("Current enrichment is APOE-driven, so it is better treated as a clue than a pathway conclusion.", 14, AMBER, True)
]))
add_slide(els, images=[img] if img.exists() else [])


# 10. Relaxed ORA results
els = title_elements(
    "Relaxed ORA reveals broader pathway signal",
    "Using FDR < 0.1, no I2 filter, and absolute effect > 0.1 gives a 105-gene list for pathway exploration.",
    "Relaxed enrichment"
)
for i, (label, value, color) in enumerate([
    ("genes", "105", BLUE),
    ("upregulated", "45", GREEN),
    ("downregulated", "60", RED),
    ("ORA terms", str(len(gprofiler_relaxed)), AMBER),
]):
    x = emu(0.85 + i * 3.05)
    els.append(text_box(350 + i, f"relaxed-metric-{i}", x, emu(1.95), emu(2.45), emu(1.05), [
        p_xml(value, 27, color, True, "ctr"),
        p_xml(label, 11.5, INK, True, "ctr")
    ], fill="FFFFFF", line=color, radius=True))

relaxed_query_counts = {}
for row in gprofiler_relaxed:
    relaxed_query_counts[row.get("query_name", "")] = relaxed_query_counts.get(row.get("query_name", ""), 0) + 1

els.append(text_box(360, "relaxed-breakdown", emu(0.95), emu(3.35), emu(4.2), emu(1.2), [
    p_xml("Where the 66 ORA terms came from", 15, INK, True),
    p_xml(f"All relaxed genes: {relaxed_query_counts.get('relaxed_no_i2_abs0.1_all_significant', 0)} terms", 12.5, MUTED),
    p_xml(f"Upregulated genes: {relaxed_query_counts.get('relaxed_no_i2_abs0.1_upregulated', 0)} terms", 12.5, MUTED),
    p_xml(f"Downregulated genes: {relaxed_query_counts.get('relaxed_no_i2_abs0.1_downregulated', 0)} terms", 12.5, MUTED),
]))

terms = top_terms(gprofiler_relaxed, query_name="relaxed_no_i2_abs0.1_all_significant", n=6)
els.append(text_box(370, "relaxed-top-title", emu(5.65), emu(3.18), emu(6.7), emu(0.35), [
    p_xml("Top relaxed all-gene terms", 15, INK, True)
]))
for i, row in enumerate(terms):
    y = emu(3.65 + i * 0.42)
    term_text = f"{row.get('term_name', '')} ({row.get('source', '')}), p={fmt_float(row.get('p_value', 'NA'), 2)}"
    els.append(text_box(380 + i, f"relaxed-term-{i}", emu(5.65), y, emu(6.7), emu(0.32), [
        p_xml(term_text, 10.8, INK if i < 3 else MUTED)
    ], fill="F8FAFC" if i % 2 == 0 else "FFFFFF"))

els.append(text_box(390, "relaxed-caution", emu(1.0), emu(6.15), emu(11.2), emu(0.58), [
    p_xml("Interpretation: the relaxed set gives better pathway discovery power, but it is exploratory because the I2 filter is removed and the effect-size threshold is lower.", 13.5, AMBER, True, "ctr")
]))
add_slide(els)


# 11. Limitations and next steps
els = title_elements("What this supports now", "A strong exploratory result, with clear manual-review steps before publication.", "Takeaway")
for i, (head, body, color) in enumerate([
    ("Defensible claim", "A small core set of adipose genes is meta-analytically associated with exercise/lifestyle intervention contexts.", GREEN),
    ("Main caveat", "The studies are related but not identical; this is not yet a pure endurance-training signature.", RED),
    ("Next step", "Use the audit tables, forest plots, and diet-arm sensitivity to decide what belongs in a manuscript.", BLUE),
]):
    y = emu(2.0 + i * 1.35)
    els.append(rect(310 + i, f"take-rule-{i}", emu(0.9), y + emu(0.1), emu(0.09), emu(0.75), color))
    els.append(text_box(320 + i, f"take-{i}", emu(1.1), y, emu(10.9), emu(0.9), [
        p_xml(head, 18, INK, True),
        p_xml(body, 14.5, MUTED)
    ]))
els.append(text_box(340, "outputs", emu(1.1), emu(6.15), emu(10.6), emu(0.45), [
    p_xml("Companion files: summary PDF, strict/relaxed CSVs, forest plots, leave-one-out figures, manuscript tables.", 12.5, MUTED, False, "ctr")
]))
add_slide(els)


def rels_xml(slide_count):
    rels = []
    for i in range(slide_count):
        rels.append(
            f'<Relationship Id="rId{i+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide{i+1}.xml"/>'
        )
    rels.append(
        f'<Relationship Id="rId{slide_count+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="slideMasters/slideMaster1.xml"/>'
    )
    rels.append(
        f'<Relationship Id="rId{slide_count+2}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/presProps" Target="presProps.xml"/>'
    )
    rels.append(
        f'<Relationship Id="rId{slide_count+3}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/viewProps" Target="viewProps.xml"/>'
    )
    rels.append(
        f'<Relationship Id="rId{slide_count+4}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/tableStyles" Target="tableStyles.xml"/>'
    )
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{''.join(rels)}</Relationships>"""


def presentation_xml(slide_count):
    ids = []
    for i in range(slide_count):
        ids.append(f'<p:sldId id="{256+i}" r:id="rId{i+1}"/>')
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
 <p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId{slide_count+1}"/></p:sldMasterIdLst>
 <p:sldIdLst>{''.join(ids)}</p:sldIdLst>
 <p:sldSz cx="{SLIDE_W}" cy="{SLIDE_H}" type="wide"/>
 <p:notesSz cx="6858000" cy="9144000"/>
 <p:defaultTextStyle>
  <a:defPPr><a:defRPr lang="en-US"/></a:defPPr>
  <a:lvl1pPr marL="0" algn="l" defTabSz="914400" rtl="0" eaLnBrk="1" fontAlgn="base" hangingPunct="1"><a:defRPr sz="1800" kern="1200"><a:solidFill><a:schemeClr val="tx1"/></a:solidFill><a:latin typeface="Aptos"/></a:defRPr></a:lvl1pPr>
 </p:defaultTextStyle>
</p:presentation>"""


def content_types(slide_count, images):
    overrides = [
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Default Extension="png" ContentType="image/png"/>',
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
        '<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>',
        '<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>',
        '<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>',
        '<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>',
        '<Override PartName="/ppt/presProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presProps+xml"/>',
        '<Override PartName="/ppt/viewProps.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.viewProps+xml"/>',
        '<Override PartName="/ppt/tableStyles.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.tableStyles+xml"/>',
    ]
    for i in range(slide_count):
        overrides.append(f'<Override PartName="/ppt/slides/slide{i+1}.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>')
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">{''.join(overrides)}</Types>"""


ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="ppt/presentation.xml"/>
 <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
 <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>"""


THEME = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="MetaAnalysisTheme">
 <a:themeElements>
  <a:clrScheme name="Custom"><a:dk1><a:srgbClr val="172033"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="334155"/></a:dk2><a:lt2><a:srgbClr val="F8FAFC"/></a:lt2><a:accent1><a:srgbClr val="2563EB"/></a:accent1><a:accent2><a:srgbClr val="0F766E"/></a:accent2><a:accent3><a:srgbClr val="D97706"/></a:accent3><a:accent4><a:srgbClr val="16A34A"/></a:accent4><a:accent5><a:srgbClr val="DC2626"/></a:accent5><a:accent6><a:srgbClr val="64748B"/></a:accent6><a:hlink><a:srgbClr val="2563EB"/></a:hlink><a:folHlink><a:srgbClr val="7C3AED"/></a:folHlink></a:clrScheme>
  <a:fontScheme name="Custom"><a:majorFont><a:latin typeface="Aptos Display"/></a:majorFont><a:minorFont><a:latin typeface="Aptos"/></a:minorFont></a:fontScheme>
  <a:fmtScheme name="Custom"><a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst><a:lnStyleLst><a:ln w="9525"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst><a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst><a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme>
 </a:themeElements>
</a:theme>"""


def slide_rels(images):
    rels = [
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>'
    ]
    for idx, target in enumerate(images, start=2):
        rels.append(
            f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/{target.name}"/>'
        )
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{''.join(rels)}</Relationships>"""


SLIDE_MASTER = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
 <p:cSld>
  <p:spTree>
   <p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>
   <p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>
  </p:spTree>
 </p:cSld>
 <p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>
 <p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>
 <p:txStyles>
  <p:titleStyle><a:lvl1pPr algn="l"><a:defRPr sz="4400" b="1"><a:latin typeface="Aptos Display"/></a:defRPr></a:lvl1pPr></p:titleStyle>
  <p:bodyStyle><a:lvl1pPr algn="l"><a:defRPr sz="1800"><a:latin typeface="Aptos"/></a:defRPr></a:lvl1pPr></p:bodyStyle>
  <p:otherStyle><a:lvl1pPr algn="l"><a:defRPr sz="1800"><a:latin typeface="Aptos"/></a:defRPr></a:lvl1pPr></p:otherStyle>
 </p:txStyles>
</p:sldMaster>"""


SLIDE_MASTER_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>
 <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme" Target="../theme/theme1.xml"/>
</Relationships>"""


SLIDE_LAYOUT = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" type="blank" preserve="1">
 <p:cSld name="Blank">
  <p:spTree>
   <p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>
   <p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>
  </p:spTree>
 </p:cSld>
 <p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr>
</p:sldLayout>"""


SLIDE_LAYOUT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster" Target="../slideMasters/slideMaster1.xml"/>
</Relationships>"""


PRES_PROPS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentationPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"/>"""


VIEW_PROPS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:viewPr xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
 <p:normalViewPr><p:restoredLeft sz="15620"/><p:restoredTop sz="94660"/></p:normalViewPr>
 <p:slideViewPr><p:cSldViewPr><p:cViewPr varScale="1"><p:scale><a:sx n="100" d="100"/><a:sy n="100" d="100"/></p:scale><p:origin x="0" y="0"/></p:cViewPr><p:guideLst/></p:cSldViewPr></p:slideViewPr>
 <p:notesTextViewPr><p:cViewPr><p:scale><a:sx n="100" d="100"/><a:sy n="100" d="100"/></p:scale><p:origin x="0" y="0"/></p:cViewPr></p:notesTextViewPr>
 <p:gridSpacing cx="72008" cy="72008"/>
</p:viewPr>"""


TABLE_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:tblStyleLst xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" def="{5C22544A-7EE6-4342-B048-85BDC9FD1C3A}"/>"""


CORE_PROPS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/"
 xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/"
 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
 <dc:title>Meta-analysis results presentation</dc:title>
 <dc:creator>Codex</dc:creator>
 <cp:lastModifiedBy>Codex</cp:lastModifiedBy>
 <dcterms:created xsi:type="dcterms:W3CDTF">2026-06-11T00:00:00Z</dcterms:created>
 <dcterms:modified xsi:type="dcterms:W3CDTF">2026-06-11T00:00:00Z</dcterms:modified>
</cp:coreProperties>"""


def app_props(slide_count):
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
 <Application>Microsoft PowerPoint</Application>
 <PresentationFormat>On-screen Show (16:9)</PresentationFormat>
 <Slides>{slide_count}</Slides>
 <ScaleCrop>false</ScaleCrop>
 <HeadingPairs><vt:vector size="2" baseType="variant"><vt:variant><vt:lpstr>Theme</vt:lpstr></vt:variant><vt:variant><vt:i4>1</vt:i4></vt:variant></vt:vector></HeadingPairs>
 <TitlesOfParts><vt:vector size="1" baseType="lpstr"><vt:lpstr>Office Theme</vt:lpstr></vt:vector></TitlesOfParts>
</Properties>"""


media_files = {}
for slide in slides:
    copied = []
    for img in slide["images"]:
        if not img.exists():
            continue
        name = img.name
        if name in media_files and media_files[name] != img:
            stem, suffix = img.stem, img.suffix
            name = f"{stem}_{len(media_files)+1}{suffix}"
        media_files[name] = img
        copied.append(Path(name))
    slide["rel_images"] = copied

tmp_dir = OUT_DIR / "_pptx_build"
if tmp_dir.exists():
    shutil.rmtree(tmp_dir)
(tmp_dir / "_rels").mkdir(parents=True)
(tmp_dir / "ppt" / "_rels").mkdir(parents=True)
(tmp_dir / "ppt" / "slides" / "_rels").mkdir(parents=True)
(tmp_dir / "ppt" / "slideMasters" / "_rels").mkdir(parents=True)
(tmp_dir / "ppt" / "slideLayouts" / "_rels").mkdir(parents=True)
(tmp_dir / "ppt" / "media").mkdir(parents=True)
(tmp_dir / "ppt" / "theme").mkdir(parents=True)
(tmp_dir / "docProps").mkdir(parents=True)

(tmp_dir / "[Content_Types].xml").write_text(content_types(len(slides), media_files), encoding="utf-8")
(tmp_dir / "_rels" / ".rels").write_text(ROOT_RELS, encoding="utf-8")
(tmp_dir / "ppt" / "presentation.xml").write_text(presentation_xml(len(slides)), encoding="utf-8")
(tmp_dir / "ppt" / "_rels" / "presentation.xml.rels").write_text(rels_xml(len(slides)), encoding="utf-8")
(tmp_dir / "ppt" / "theme" / "theme1.xml").write_text(THEME, encoding="utf-8")
(tmp_dir / "ppt" / "slideMasters" / "slideMaster1.xml").write_text(SLIDE_MASTER, encoding="utf-8")
(tmp_dir / "ppt" / "slideMasters" / "_rels" / "slideMaster1.xml.rels").write_text(SLIDE_MASTER_RELS, encoding="utf-8")
(tmp_dir / "ppt" / "slideLayouts" / "slideLayout1.xml").write_text(SLIDE_LAYOUT, encoding="utf-8")
(tmp_dir / "ppt" / "slideLayouts" / "_rels" / "slideLayout1.xml.rels").write_text(SLIDE_LAYOUT_RELS, encoding="utf-8")
(tmp_dir / "ppt" / "presProps.xml").write_text(PRES_PROPS, encoding="utf-8")
(tmp_dir / "ppt" / "viewProps.xml").write_text(VIEW_PROPS, encoding="utf-8")
(tmp_dir / "ppt" / "tableStyles.xml").write_text(TABLE_STYLES, encoding="utf-8")
(tmp_dir / "docProps" / "core.xml").write_text(CORE_PROPS, encoding="utf-8")
(tmp_dir / "docProps" / "app.xml").write_text(app_props(len(slides)), encoding="utf-8")

for idx, slide in enumerate(slides, start=1):
    (tmp_dir / "ppt" / "slides" / f"slide{idx}.xml").write_text(slide["xml"], encoding="utf-8")
    rel = slide_rels(slide["rel_images"])
    if rel:
        (tmp_dir / "ppt" / "slides" / "_rels" / f"slide{idx}.xml.rels").write_text(rel, encoding="utf-8")

for name, src in media_files.items():
    shutil.copyfile(src, tmp_dir / "ppt" / "media" / name)

if PPTX_PATH.exists():
    PPTX_PATH.unlink()
with zipfile.ZipFile(PPTX_PATH, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for path in tmp_dir.rglob("*"):
        if path.is_file():
            z.write(path, path.relative_to(tmp_dir).as_posix())

shutil.rmtree(tmp_dir)
print(PPTX_PATH)
