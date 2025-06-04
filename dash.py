import streamlit as st
import os.path as path
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scipy.stats as stats

# Local path
if getattr(sys, 'frozen', False):
    this_path = path.dirname(sys.executable)
elif '_file_' in globals():
    this_path = path.abspath(path.dirname(_file_))


# -------------------------
# FUNKCJE POMOCNICZE
# -------------------------

def replace_zero_dose_with_replicas(df: pd.DataFrame, dose_column: str = 'Dawka (MBq/ml)', new_doses: list = [50, 100, 200]) -> pd.DataFrame:
    zero_dose_rows = df[df[dose_column] == 0]
    replicated_rows = pd.concat([
        zero_dose_rows.assign(**{dose_column: dose}) for dose in new_doses
    ], ignore_index=True)
    df_without_zero = df[df[dose_column] != 0]
    return pd.concat([df_without_zero, replicated_rows], ignore_index=True)


def draw_detailed_barplot(data, cell_line="MDA-MB-231", hypothesis_group=None):
    st.write("## Wykres słupkowy z oznaczeniem istotności (test jednostronny)")
    dose_col = "Dawka (MBq/ml)"
    viability_col = "Aktywność metaboliczna (%)"
    data["Grupa"] = data["Grupa"].str.strip()
    data_line = data[data["Linia"] == cell_line]
    times = sorted(data_line["Czas (h)"].dropna().unique())
    groups = sorted(data_line["Grupa"].dropna().unique())
    non_zero_doses = sorted([d for d in data_line[dose_col].dropna().unique() if d != 0.0])

    fig, axs = plt.subplots(1, len(non_zero_doses), figsize=(6 * len(non_zero_doses), 6), sharey=True)
    if len(non_zero_doses) == 1:
        axs = [axs]

    for idx_d, dose in enumerate(non_zero_doses):
        data_dose = data_line[(data_line[dose_col] == dose) | (data_line[dose_col] == 0)]
        ax = axs[idx_d]
        bar_width = 0.15
        x = np.arange(len(times))
        bar_positions = {}
        bar_heights = {}

        for idx_g, group in enumerate(groups):
            center = x + (idx_g - len(groups) / 2) * bar_width + bar_width / 2
            bar_positions[group] = center
            values = []
            errors = []
            for time in times:
                subset = data_dose[
                    (data_dose["Grupa"] == group) &
                    (data_dose["Czas (h)"] == time)
                    ][viability_col]
                values.append(subset.mean())
                errors.append(subset.sem())
            ax.bar(center, values, yerr=errors, width=bar_width, label=group, capsize=3)
            bar_heights[group] = values

        for time_idx, time in enumerate(times):
            hypo_vals = data_dose[
                (data_dose[dose_col] == dose) &
                (data_dose["Czas (h)"] == time) &
                (data_dose["Grupa"] == hypothesis_group)
                ][viability_col]

            offset_base = data_dose[viability_col].max() * 0.03 if not data_dose.empty else 5
            level = 0

            for group in groups:
                if group == hypothesis_group:
                    continue
                comp_vals = data_dose[
                    (data_dose[dose_col] == dose) &
                    (data_dose["Czas (h)"] == time) &
                    (data_dose["Grupa"] == group)
                    ][viability_col]

                if len(hypo_vals) < 2 or len(comp_vals) < 2:
                    continue

                t_stat, p_two_sided = stats.ttest_ind(hypo_vals, comp_vals, equal_var=False)
                p_one_sided = p_two_sided / 2 if t_stat < 0 else 1 - (p_two_sided / 2)

                stars = (
                    "***" if p_one_sided < 0.001 else
                    "**" if p_one_sided < 0.01 else
                    "*" if p_one_sided < 0.05 else ""
                )
                annotation = stars if stars else "n.s."

                try:
                    x1 = bar_positions[hypothesis_group][time_idx]
                    x2 = bar_positions[group][time_idx]
                    y = max(bar_heights[hypothesis_group][time_idx], bar_heights[group][time_idx]) + (level + 1) * offset_base

                    ax.plot([x1, x2], [y, y], lw=1.5, c='black')
                    ax.plot([x1, x1], [y, y - 1], lw=1.0, c='black')
                    ax.plot([x2, x2], [y, y - 1], lw=1.0, c='black')
                    ax.text((x1 + x2) / 2, y + 1.0, annotation, ha='center', va='bottom', fontsize=12, color='gray' if annotation == "n.s." else 'black')
                    level += 1
                except Exception as e:
                    st.warning(f"Błąd rysowania: {e}")
                    continue

        ax.set_title(f"{chr(65 + idx_d)}. {dose} MBq/ml (z kontrolą 0)")
        ax.set_xticks(x)
        ax.set_xticklabels([str(t) for t in times])
        ax.set_xlabel("Czas [h]")
        ax.set_ylabel("Aktywność metaboliczna [%]")
        ax.legend(title="Grupa")
        ax.grid(True, axis='y')

    plt.tight_layout()
    st.pyplot(fig)

    st.markdown("---")
    st.markdown("### Interpretacja oznaczeń istotności:")
    st.markdown("""
        - ✴️ n.s. – brak istotności (p ≥ 0.05)
        - ⭐ * – p < 0.05
        - ⭐⭐ ** – p < 0.01
        - ⭐⭐⭐ * – p < 0.001  
        Test jednostronny: czy grupa hipotezy ma niższą aktywność niż pozostałe grupy.
    """)

# -------------------------
# GŁÓWNA APLIKACJA
# -------------------------

st.title("Wersja 7: Analiza cytotoksyczności – koniugat vs komponenty")

uploaded_file = st.file_uploader("Wgraj plik DANE.xlsx", type=["xlsx"])
if uploaded_file:
    mda_detailed = pd.read_excel(uploaded_file, sheet_name="MDA detailed data")
    t98g_detailed = pd.read_excel(uploaded_file, sheet_name="T98G detailed data")

    for df in [mda_detailed, t98g_detailed]:
        df["Grupa"] = df["Grupa"].str.strip()
        df["Czas (h)"] = df["Czas (h)"].astype(int)
        df["Dawka (MBq/ml)"] = df["Dawka (MBq/ml)"].astype(float)
        df["Viability (%)"] = df["Viability (%)"].astype(float)
        df["Powtórzenie"] = df["Powtórzenie"].astype(int)

    mda_detailed["Linia"] = "MDA-MB-231"
    mda_detailed["Aktywność metaboliczna (%)"] = mda_detailed["Viability (%)"]

    t98g_detailed["Linia"] = "T98G"
    t98g_detailed["Aktywność metaboliczna (%)"] = t98g_detailed["Viability (%)"]

    data = pd.concat([mda_detailed, t98g_detailed], ignore_index=True)
    data = replace_zero_dose_with_replicas(data)

    cell_line = st.selectbox("Wybierz linię komórkową", data["Linia"].unique())
    sorted_list_of_groups = sorted(data["Grupa"].unique(), key=len, reverse=True)
    hypothesis_group = st.selectbox("Wybierz grupę do hipotezy", sorted_list_of_groups)

    selected_data = data[data["Linia"] == cell_line]
    st.write("### Podgląd danych")
    st.dataframe(selected_data)

    st.write("## Analiza statystyczna (dla każdej dawki osobno)")
    for dose in sorted(selected_data["Dawka (MBq/ml)"].unique()):
        st.markdown(f"### Dawka: {dose} MBq/ml")
        dose_data = selected_data[selected_data["Dawka (MBq/ml)"] == dose]
        groups = dose_data["Grupa"].unique()
        data_list = [dose_data[dose_data["Grupa"] == g]["Aktywność metaboliczna (%)"] for g in groups]

        if any(len(x) < 2 for x in data_list):
            st.warning("Za mało danych dla niektórych grup — pominięto analizę.")
            continue

        f_val, p_val = stats.f_oneway(*data_list)
        st.write(f"*ANOVA:* F = {f_val:.2f}, p = {p_val:.4f}")
        st.markdown("- Testuje, czy średnia aktywność różni się w którejkolwiek z grup.")
        if p_val < 0.05:
            st.success("Istnieją statystycznie istotne różnice między grupami.")
        else:
            st.info("Brak istotnych różnic między grupami.")

        tukey = pairwise_tukeyhsd(endog=dose_data["Aktywność metaboliczna (%)"], groups=dose_data["Grupa"], alpha=0.05)
        st.write("*Test post-hoc: Tukey HSD*")
        st.text(tukey.summary())
        st.markdown("- Pokazuje, które dokładnie grupy różnią się między sobą.")

    st.write("## Test jednostronny: Czy grupa hipotezy ma niższą aktywność?")
    for dose in sorted(selected_data["Dawka (MBq/ml)"].unique()):
        dose_data = selected_data[selected_data["Dawka (MBq/ml)"] == dose]
        st.markdown(f"### Dawka: {dose} MBq/ml")
        other_groups = [g for g in dose_data["Grupa"].unique() if g != hypothesis_group]

        for other in other_groups:
            hypo_vals = dose_data[dose_data["Grupa"] == hypothesis_group]["Aktywność metaboliczna (%)"]
            comp_vals = dose_data[dose_data["Grupa"] == other]["Aktywność metaboliczna (%)"]

            if len(hypo_vals) < 2 or len(comp_vals) < 2:
                st.markdown(f"- {hypothesis_group} vs {other}: za mało danych.")
                continue

            t_stat, p_two_sided = stats.ttest_ind(hypo_vals, comp_vals, equal_var=False)
            p_one_sided = p_two_sided / 2 if t_stat < 0 else 1 - (p_two_sided / 2)
            significant = "istotna" if p_one_sided < 0.05 else "nieistotna"
            st.markdown(f"- {hypothesis_group} vs {other}: t = {t_stat:.2f}, p = {p_one_sided:.4f} → różnica {significant}")

    draw_detailed_barplot(selected_data, cell_line=cell_line, hypothesis_group=hypothesis_group)