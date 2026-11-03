"""Generate and embed a Python/matplotlib root-locus diagram in the HTML guide."""
from __future__ import annotations

import base64
import io
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HTML = ROOT / "HTML" / "control_systems_design_reference.html"

# Unity-feedback plant G(s)=1/[s(s+1)]: s^2+s+K=0.
gain = np.linspace(0.0, 2.0, 600)
disc = 1.0 - 4.0 * gain
root_a = (-1.0 + np.sqrt(disc.astype(complex))) / 2.0
root_b = (-1.0 - np.sqrt(disc.astype(complex))) / 2.0

fig, ax = plt.subplots(figsize=(8.2, 3.4), constrained_layout=True)
ax.plot(root_a.real, root_a.imag, color="#2563eb", lw=2.5)
ax.plot(root_b.real, root_b.imag, color="#2563eb", lw=2.5)
ax.scatter([0, -1], [0, 0], marker="x", s=80, linewidths=2, color="#b91c1c", label="open-loop poles")
ax.scatter([-0.5], [0], marker="o", s=55, facecolors="none", edgecolors="#15803d", label="breakaway")
ax.axhline(0, color="#64748b", lw=0.8)
ax.axvline(0, color="#64748b", lw=0.8)
ax.set(xlabel="Real axis, Re{s}", ylabel="Imaginary axis, Im{s}",
       title="Root locus of G(s) = 1 / [s(s + 1)] as K increases")
ax.grid(True, alpha=0.22)
ax.legend(loc="upper left", frameon=False, fontsize=8)
ax.set_xlim(-2.2, 0.25)
ax.set_ylim(-1.8, 1.8)

buffer = io.BytesIO()
fig.savefig(buffer, format="svg", metadata={"Date": None})
plt.close(fig)
svg = buffer.getvalue()
encoded = base64.b64encode(svg).decode("ascii")
root_image = f'<div class="fig"><img alt="Python-generated root locus for G(s)=1/[s(s+1)]" style="width:100%;height:auto;border:1px solid #d8e0ea;border-radius:10px" src="data:image/svg+xml;base64,{encoded}"></div>'

# A matching Bode diagram for the same plant, also generated with Python.
w = np.logspace(-2, 2, 500)
mag_db = 20 * np.log10(1 / (w * np.sqrt(w**2 + 1)))
phase = -90 - np.degrees(np.arctan(w))
fig, (axm, axp) = plt.subplots(2, 1, figsize=(8.2, 4.4), sharex=True, constrained_layout=True)
axm.semilogx(w, mag_db, color="#2563eb", lw=2.2)
axm.set_ylabel("Magnitude (dB)")
axm.grid(True, which="both", alpha=0.22)
axm.set_title("Bode response of G(s) = 1 / [s(s + 1)]")
axp.semilogx(w, phase, color="#ea580c", lw=2.2)
axp.set(xlabel="Frequency (rad/s)", ylabel="Phase (degrees)", ylim=(-190, -75))
axp.grid(True, which="both", alpha=0.22)
buffer = io.BytesIO()
fig.savefig(buffer, format="svg", metadata={"Date": None})
plt.close(fig)
bode_encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
bode_image = f'<div class="fig"><img alt="Python-generated Bode plot for G(s)=1/[s(s+1)]" style="width:100%;height:auto;border:1px solid #d8e0ea;border-radius:10px" src="data:image/svg+xml;base64,{bode_encoded}"></div>'

document = HTML.read_text(encoding="utf-8")
updated = document
pattern = r'<div class="fig"><svg viewBox="0 0 700 240".*?</svg></div>'
if 'alt="Python-generated root locus' not in updated:
    updated, count = re.subn(pattern, lambda _: root_image, updated, count=1, flags=re.DOTALL)
    if count != 1:
        raise RuntimeError("Could not find the root-locus placeholder in the HTML document")
if 'alt="Python-generated Bode plot' not in updated:
    root_pattern = r'<div class="fig"><img alt="Python-generated root locus[^>]*></div>'
    roots = list(re.finditer(root_pattern, updated))
    if len(roots) >= 2:
        match = roots[1]
        updated = updated[:match.start()] + bode_image + updated[match.end():]
    else:
        updated, count = re.subn(pattern, lambda _: bode_image, updated, count=1, flags=re.DOTALL)
        if count != 1:
            raise RuntimeError("Could not find the Bode placeholder in the HTML document")

# Restore explicit inline math delimiters around equation-sized parenthetical
# expressions so MathJax renders the formulas rather than browser text.
def normalize_math(paragraph: str) -> str:
    result: list[str] = []
    inside = False
    nested_parens = 0
    i = 0
    while i < len(paragraph):
        if paragraph.startswith(r"\(", i):
            if inside:
                result.append("(")
                nested_parens += 1
            else:
                result.append(r"\(")
                inside = True
                nested_parens = 0
            i += 2
        elif paragraph.startswith(r"\)", i):
            if inside and nested_parens:
                result.append(")")
                nested_parens -= 1
            elif inside:
                result.append(r"\)")
                inside = False
            else:
                result.append(")")
            i += 2
        else:
            char = paragraph[i]
            result.append(char)
            if inside and char == "(":
                nested_parens += 1
            elif inside and char == ")" and nested_parens:
                nested_parens -= 1
            i += 1
    return "".join(result)

def delimit_equations(paragraph: str) -> str:
    result: list[str] = []
    i = 0
    while i < len(paragraph):
        if paragraph[i] != "(" or (i > 0 and paragraph[i - 1] == "\\"):
            result.append(paragraph[i])
            i += 1
            continue
        depth = 0
        end = i
        for end in range(i, len(paragraph)):
            if paragraph[end] == "(":
                depth += 1
            elif paragraph[end] == ")":
                depth -= 1
                if depth == 0:
                    break
        if depth != 0:
            result.append(paragraph[i])
            i += 1
            continue
        expression = paragraph[i + 1:end]
        if re.search(r"(?:=|\\[A-Za-z]+|[_^]|[<>]|≈|±|\d\s*/\s*\d)", expression):
            result.append(r"\(" + expression + r"\)")
        else:
            result.append(paragraph[i:end + 1])
        i = end + 1
    return "".join(result)

updated = re.sub(r"<p>(.*?)</p>", lambda match: "<p>" + normalize_math(match.group(1)) + "</p>", updated, flags=re.DOTALL)
updated = updated.replace(r"\int e,dt", r"\int e(t)\,dt")
updated = updated.replace("ẋs", r"\dot{s}")
HTML.write_text(updated, encoding="utf-8")
print(f"Updated Python-generated diagrams and math in {HTML}")
