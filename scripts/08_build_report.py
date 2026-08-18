#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Render RESULTS.md into a standalone results/report.html.

RESULTS.md is the source of truth; this only styles it and inlines the figure as a data
URI, so the page can be shared without carrying the folder. Run it after editing the
report or regenerating the figure.
"""

import argparse
import base64
import re
from pathlib import Path

import markdown
from PIL import Image

BENCH_ROOT = Path(__file__).resolve().parent.parent

STYLE = """
  :root {
    color-scheme: light;
    --ground:  #fbfcfb;
    --surface: #f1f5f3;
    --ink:     #131c19;
    --muted:   #5a6b64;
    --rule:    #dce4e0;
    --accent:  #0f7d61;
    --shadow:  0 1px 2px rgba(19, 28, 25, .06), 0 8px 24px rgba(19, 28, 25, .05);
    --serif: ui-serif, "Iowan Old Style", "Source Serif 4", Palatino, Georgia, serif;
    --sans: ui-sans-serif, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
    --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      color-scheme: dark;
      --ground:  #0e1412;
      --surface: #161f1b;
      --ink:     #edf3f0;
      --muted:   #9aaba3;
      --rule:    #26312d;
      --accent:  #3cbe93;
      --shadow:  0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px rgba(0, 0, 0, .3);
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --ground:  #0e1412;
    --surface: #161f1b;
    --ink:     #edf3f0;
    --muted:   #9aaba3;
    --rule:    #26312d;
    --accent:  #3cbe93;
    --shadow:  0 1px 2px rgba(0, 0, 0, .4), 0 8px 24px rgba(0, 0, 0, .3);
  }
  body {
    background: var(--ground);
    color: var(--ink);
    font-family: var(--sans);
    font-size: 16.5px;
    line-height: 1.62;
    margin: 0;
    padding: clamp(2rem, 5vw, 4.5rem) 1.25rem 6rem;
  }
  main { max-width: 74ch; margin: 0 auto; }
  main > * { margin: 0 0 1.15rem; }
  h1 {
    font-family: var(--serif);
    font-size: clamp(2rem, 4.4vw, 2.8rem);
    line-height: 1.12;
    letter-spacing: -.015em;
    font-weight: 600;
    text-wrap: balance;
    margin-bottom: 1.6rem;
  }
  h2 {
    font-family: var(--serif);
    font-size: 1.5rem;
    font-weight: 600;
    letter-spacing: -.01em;
    border-top: 1px solid var(--rule);
    padding-top: 1.5rem;
    margin-top: 2.4rem;
    text-wrap: balance;
  }
  h3 {
    font-family: var(--mono);
    font-size: .82rem;
    text-transform: uppercase;
    letter-spacing: .07em;
    color: var(--accent);
    font-weight: 600;
    margin-top: 2rem;
  }
  strong { font-weight: 630; }
  a { color: var(--accent); }
  code {
    font-family: var(--mono);
    font-size: .88em;
    background: var(--surface);
    padding: .1em .34em;
    border-radius: 3px;
  }
  pre {
    font-family: var(--mono);
    font-size: .84rem;
    background: var(--surface);
    border: 1px solid var(--rule);
    border-radius: 6px;
    padding: .9rem 1rem;
    overflow-x: auto;
    line-height: 1.55;
  }
  pre code { background: none; padding: 0; font-size: 1em; }
  ul, ol { padding-left: 1.3rem; }
  li { margin-bottom: .45rem; }
  li::marker { color: var(--muted); }
  figure { margin: 2.6rem 0; }
  img {
    display: block;
    width: 100%;
    max-width: 1140px;
    height: auto;
    margin: 0 auto;
    border: 1px solid var(--rule);
    border-radius: 6px;
    background: #fff;
    box-shadow: var(--shadow);
  }
  .tablewrap {
    overflow-x: auto;
    border: 1px solid var(--rule);
    border-radius: 6px;
    background: var(--surface);
    margin-bottom: 1.15rem;
  }
  table { border-collapse: collapse; width: 100%; font-size: .92rem; }
  th, td { padding: .48rem .9rem; text-align: right; white-space: nowrap; }
  th:first-child, td:first-child { text-align: left; white-space: normal; }
  thead th {
    font-family: var(--mono);
    font-size: .74rem;
    letter-spacing: .05em;
    text-transform: uppercase;
    color: var(--muted);
    font-weight: 600;
    border-bottom: 1px solid var(--rule);
  }
  tbody td { font-variant-numeric: tabular-nums; border-bottom: 1px solid var(--rule); }
  tbody tr:last-child td { border-bottom: 0; }
  :focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=str(BENCH_ROOT / "RESULTS.md"))
    parser.add_argument("--out", default=str(BENCH_ROOT / "results" / "report.html"))
    parser.add_argument("--title", default="SATIVA on EPA-ng")
    args = parser.parse_args()

    text = Path(args.source).read_text()

    # The figure travels inside the page: half size, base64, as a data URI.
    def inline_image(match):
        alt, rel = match.group(1), match.group(2)
        path = (BENCH_ROOT / rel).resolve()
        web = BENCH_ROOT / "results" / "figures" / "_web.png"
        image = Image.open(path)
        image.resize((image.width // 2, image.height // 2), Image.LANCZOS).save(web, optimize=True)
        data = base64.b64encode(web.read_bytes()).decode()
        web.unlink()
        return f'<figure><img src="data:image/png;base64,{data}" alt="{alt}"></figure>'

    body = markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"])
    body = re.sub(r'<p><img alt="([^"]*)" src="([^"]+)"\s*/?></p>', inline_image, body)
    body = body.replace("<table>", '<div class="tablewrap"><table>')
    body = body.replace("</table>", "</table></div>")

    out = Path(args.out)
    out.write_text(f"<title>{args.title}</title>\n<style>{STYLE}</style>\n<main>\n{body}\n</main>\n")
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
