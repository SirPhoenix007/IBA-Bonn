# extract_pixe_pdfplumber_v2.py
from pathlib import Path
import re
import pdfplumber
import pandas as pd

PDF_PATH = Path("PIXE_ECPSSR.pdf")
OUT_DIR = Path("output_tables")
OUT_DIR.mkdir(exist_ok=True)

FLOAT = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?"

TARGET_RE = re.compile(r"^(?P<z>\d{1,3})-(?P<sym>[A-Za-z]{1,3})$")

# Target optional, danach 6 Zahlen.
TABLE_I_WITH_TARGET_RE = re.compile(
    rf"^(?P<target>\d{{1,3}}-[A-Za-z]{{1,3}})\s+"
    rf"(?P<energy>{FLOAT})\s+"
    rf"(?P<logxi>{FLOAT})\s+"
    rf"(?P<ecpssr>{FLOAT})\s+"
    rf"(?P<ref>{FLOAT})\s+"
    rf"(?P<rnd>{FLOAT})\s+"
    rf"(?P<tot>{FLOAT})$"
)

TABLE_I_NO_TARGET_RE = re.compile(
    rf"^(?P<energy>{FLOAT})\s+"
    rf"(?P<logxi>{FLOAT})\s+"
    rf"(?P<ecpssr>{FLOAT})\s+"
    rf"(?P<ref>{FLOAT})\s+"
    rf"(?P<rnd>{FLOAT})\s+"
    rf"(?P<tot>{FLOAT})$"
)

# Tabelle II: Z2 ENERGY SIGMA ERROR REFERENCE
TABLE_II_RE = re.compile(
    rf"^(?P<z2>\d{{1,3}})\s+"
    rf"(?P<energy>{FLOAT})\s+"
    rf"(?P<sigma>{FLOAT})\s+"
    rf"(?P<error>{FLOAT})\s+"
    rf"(?P<reference>[A-Z][A-Z0-9]+)$"
)

PAGE_NUMBER_LINE_RE = re.compile(r"^\d+\s*$")


def to_float(s: str):
    return float(s.replace("D", "E"))


def normalize_line(line: str) -> str:
    line = line.replace("−", "-").replace("–", "-")
    line = line.replace(",", ".")
    line = re.sub(r"\s+", " ", line.strip())
    return line


def split_possible_multicolumn_line(line: str):
    """
    Zerlegt Zeilen, in denen pdfplumber linke und rechte Tabellenspalte
    hintereinander geschrieben hat.
    """
    line = normalize_line(line)
    if not line:
        return []

    # Häufig beginnt ein neuer Tabellenblock mitten in der Zeile mit z.B. "13-Al".
    # Vor solche Targets wird ein Trenner gesetzt.
    line = re.sub(r"\s+(?=\d{1,3}-[A-Za-z]{1,3}\s)", " ||| ", line)

    # Für Tabelle II beginnt ein neuer Datensatz oft mit: Leerzeichen + Z2 + Energie + Sigma ...
    # Das ist schwieriger eindeutig; daher zusätzlich später Sliding-Window-Parsing.
    parts = [p.strip() for p in line.split("|||") if p.strip()]
    return parts


def parse_table_i_from_tokens(tokens, current_target):
    """
    Parse Tabelle I auch dann, wenn mehrere Records in einer Zeile stehen.
    Record mit Target: target + 6 numerische Tokens.
    Record ohne Target: 6 numerische Tokens und current_target muss existieren.
    """
    rows = []
    i = 0

    while i < len(tokens):
        token = tokens[i]

        target_match = TARGET_RE.match(token)
        if target_match and i + 6 < len(tokens):
            nums = tokens[i + 1:i + 7]
            if all(re.fullmatch(FLOAT, x) for x in nums):
                current_target = {
                    "target_Z": int(target_match.group("z")),
                    "target_symbol": target_match.group("sym"),
                }
                rows.append({
                    **current_target,
                    "energy_MeV": to_float(nums[0]),
                    "log_xi": to_float(nums[1]),
                    "sigma_ecpssr_barn": to_float(nums[2]),
                    "sigma_reference_barn": to_float(nums[3]),
                    "error_rnd_percent": to_float(nums[4]),
                    "error_tot_percent": to_float(nums[5]),
                })
                i += 7
                continue

        if current_target and i + 5 < len(tokens):
            nums = tokens[i:i + 6]
            if all(re.fullmatch(FLOAT, x) for x in nums):
                rows.append({
                    **current_target,
                    "energy_MeV": to_float(nums[0]),
                    "log_xi": to_float(nums[1]),
                    "sigma_ecpssr_barn": to_float(nums[2]),
                    "sigma_reference_barn": to_float(nums[3]),
                    "error_rnd_percent": to_float(nums[4]),
                    "error_tot_percent": to_float(nums[5]),
                })
                i += 6
                continue

        i += 1

    return rows, current_target


def parse_table_ii_from_tokens(tokens):
    """
    Tabelle II Record: z2, energy, sigma, error, reference.
    """
    rows = []
    i = 0

    while i + 4 < len(tokens):
        z2, energy, sigma, error, ref = tokens[i:i + 5]

        if (
            re.fullmatch(r"\d{1,3}", z2)
            and re.fullmatch(FLOAT, energy)
            and re.fullmatch(FLOAT, sigma)
            and re.fullmatch(FLOAT, error)
            and re.fullmatch(r"[A-Z][A-Z0-9]+", ref)
        ):
            rows.append({
                "z2": int(z2),
                "energy_MeV": to_float(energy),
                "sigma_barn": to_float(sigma),
                "error_percent": to_float(error),
                "reference": ref,
            })
            i += 5
        else:
            i += 1

    return rows


def classify_page(text: str):
    """
    Grobe Klassifikation anhand sichtbarer Inhalte.
    """
    if "REFERENCES FOR TABLE II" in text:
        return "refs"
    if "Experimental Cross Sections" in text or "REFERENCE" in text and "SIGMA" in text and "Z2" in text:
        return "table_ii"
    if "Empirical Reference Cross Sections" in text or "ECPSSR" in text and "REFERENCE" in text and "LOG XI" in text:
        return "table_i"
    return "unknown"


def extract():
    table_i_rows = []
    table_ii_rows = []

    current_target = None
    mode = "unknown"

    with pdfplumber.open(PDF_PATH) as pdf:
        for page_idx, page in enumerate(pdf.pages, start=1):
            text = page.extract_text(x_tolerance=1, y_tolerance=3) or ""
            page_kind = classify_page(text)

            # Modus-Logik: nach erkannter Tabelle I/II im jeweiligen Modus bleiben.
            if page_kind == "table_i":
                mode = "table_i"
            elif page_kind == "table_ii":
                mode = "table_ii"
            elif page_kind == "refs":
                mode = "refs"

            for raw_line in text.splitlines():
                line = normalize_line(raw_line)

                if not line or PAGE_NUMBER_LINE_RE.match(line):
                    continue

                # Kopf-/Fußzeilen ignorieren.
                if any(skip in line for skip in [
                    "H. PAUL", "K-Shell Ionization", "Atomic Data",
                    "TABLE I", "TABLE II", "See page", "TARGET ENERGY",
                    "LOG XI", "SIGMA", "ERROR", "REFERENCES FOR TABLE II",
                    "EXPLANATION OF TABLES"
                ]):
                    continue

                parts = split_possible_multicolumn_line(line)

                for part in parts:
                    tokens = part.split()

                    if mode == "table_i":
                        rows, current_target = parse_table_i_from_tokens(tokens, current_target)
                        for r in rows:
                            r["page"] = page_idx
                        table_i_rows.extend(rows)

                    elif mode == "table_ii":
                        rows = parse_table_ii_from_tokens(tokens)
                        for r in rows:
                            r["page"] = page_idx
                        table_ii_rows.extend(rows)

    df_i = pd.DataFrame(table_i_rows)
    df_ii = pd.DataFrame(table_ii_rows)

    return df_i, df_ii


def clean_table_i(df):
    if df.empty:
        return df

    # Plausibilitätsfilter für Tabelle I.
    df = df[
        (df["target_Z"].between(4, 92))
        & (df["energy_MeV"] > 0)
        & (df["sigma_ecpssr_barn"] > 0)
        & (df["sigma_reference_barn"] > 0)
        & (df["error_rnd_percent"] >= 0)
        & (df["error_tot_percent"] >= 0)
    ].copy()

    df = df.drop_duplicates(
        subset=[
            "target_Z", "target_symbol", "energy_MeV", "log_xi",
            "sigma_ecpssr_barn", "sigma_reference_barn",
            "error_rnd_percent", "error_tot_percent"
        ]
    )

    return df.sort_values(["target_Z", "energy_MeV"]).reset_index(drop=True)


def clean_table_ii(df):
    if df.empty:
        return df

    df = df[
        (df["z2"].between(4, 92))
        & (df["energy_MeV"] > 0)
        & (df["sigma_barn"] > 0)
        & (df["error_percent"] >= 0)
    ].copy()

    df = df.drop_duplicates(
        subset=["z2", "energy_MeV", "sigma_barn", "error_percent", "reference"]
    )

    return df.sort_values(["z2", "energy_MeV", "reference"]).reset_index(drop=True)


def main():
    df_i, df_ii = extract()
    df_i = clean_table_i(df_i)
    df_ii = clean_table_ii(df_ii)

    out_i_csv = OUT_DIR / "table_I_reference_cross_sections.csv"
    out_ii_csv = OUT_DIR / "table_II_experimental_cross_sections.csv"

    df_i.to_csv(out_i_csv, index=False)
    df_ii.to_csv(out_ii_csv, index=False)

    df_i.to_json(OUT_DIR / "table_I_reference_cross_sections.json", orient="records", indent=2)
    df_ii.to_json(OUT_DIR / "table_II_experimental_cross_sections.json", orient="records", indent=2)

    print(f"Tabelle I Zeilen:  {len(df_i)} -> {out_i_csv}")
    print(f"Tabelle II Zeilen: {len(df_ii)} -> {out_ii_csv}")

    print("\nTabelle I Vorschau:")
    print(df_i.head(10).to_string(index=False))

    print("\nTabelle II Vorschau:")
    print(df_ii.head(10).to_string(index=False))


if __name__ == "__main__":
    main()