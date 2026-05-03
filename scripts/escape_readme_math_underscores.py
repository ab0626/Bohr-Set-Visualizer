"""
Escape Markdown '_' inside GitHub $...$ and $$...$$ math as '\\_'
so CommonMark does not strip subscript underscores before KaTeX.

Usage (repo root): python scripts/escape_readme_math_underscores.py README.md
"""
from __future__ import annotations

import re
import sys


def escape_underscores_in_tex(tex: str) -> str:
    """Turn subscript '_' into '\\_' when not already escaped (for Markdown → math)."""
    s = tex
    # Braced subscripts: \mathbb{Z}_{N}, \Lambda_{\Theta,\varepsilon}
    s = re.sub(r"(?<!\\)_(?=\{)", r"\\_", s)
    # Single-token subscripts: \Theta_i, \varepsilon_0
    s = re.sub(r"(?<!\\)_(?=[a-zA-Z0-9])", r"\\_", s)
    return s


def process_markdown(text: str) -> str:
    def repl_dd(m: re.Match[str]) -> str:
        return "$$" + escape_underscores_in_tex(m.group(1)) + "$$"

    text = re.sub(r"\$\$([\s\S]*?)\$\$", repl_dd, text)

    def repl_inline(m: re.Match[str]) -> str:
        return "$" + escape_underscores_in_tex(m.group(1)) + "$"

    text = re.sub(r"\$([^$\n]+)\$", repl_inline, text)
    return text


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "README.md"
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    out = process_markdown(raw)
    if out != raw:
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(out)
        print("Updated:", path)
    else:
        print("No changes:", path)


if __name__ == "__main__":
    main()
