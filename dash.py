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

# Ścieżka lokalna
if getattr(sys, 'frozen', False):
    this_path = path.dirname(sys.executable)
elif __file__:
    this_path = path.abspath(path.dirname(__file__))

st.title("Analiza cytotoksyczności: koniugat vs komponenty")

uploaded_file = st.file_uploader("Wgraj plik DANE.xlsx", type=["xlsx"])
if uploaded_file:
    # Wczytanie tylko danych detailed
    mda_detailed = pd.read_excel(uploaded_file, sheet_name="MDA detailed data")
    t98g_detailed = pd.read_excel(uploaded_file, sheet_name="T98G detailed data")

    # Dodanie kolumn pomocniczych
    mda_detailed["Linia"] = "MDA-MB-231"
    mda_detailed["Aktywność metaboliczna (%)"] = mda_detailed["Viability (%)"]

    t98g_detailed["Linia"] = "T98G"
    t98g_detailed["Aktywność metaboliczna (%)"] = t98g_detailed["Viability (%)"]

    dane = pd.concat([mda_detailed, t98g_detailed], ignore_index=True)

    linia_wybor = st.selectbox("Wybierz linię komórkową", dane["Linia"].unique())

    dane_wybrane = dane[dane["Linia"] == linia_wybor]

    st.write("### Podgląd danych")
    st.dataframe(dane_wybrane)

    # Statystyka + wykres
    st.write("## Analiza statystyczna")

    grupy = dane_wybrane["Grupa"].unique()
    dane_list = [dane_wybrane[dane_wybrane["Grupa"] == g]["Aktywność metaboliczna (%)"] for g in grupy]

    # ANOVA
    f_val, p_val = stats.f_oneway(*dane_list)
    st.write(f"Test ANOVA: F = {f_val:.2f}, p = {p_val:.4f}")

    # Post-hoc: Tukey HSD
    tukey = pairwise_tukeyhsd(
        endog=dane_wybrane["Aktywność metaboliczna (%)"],
        groups=dane_wybrane["Grupa"],
        alpha=0.05
    )
    st.write("### Post-hoc test Tukey HSD:")
    st.text(tukey.summary())

    # Wykres 1: boxplot + scatter
    st.write("## Wykres z porównaniem grup")

    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(data=dane_wybrane, x="Grupa", y="Aktywność metaboliczna (%)")
    sns.stripplot(data=dane_wybrane, x="Grupa", y="Aktywność metaboliczna (%)", color="black", alpha=0.4)

    # Gwiazdki z Tukey HSD
    def add_stat_annotation(ax, tukey):
        ymax = dane_wybrane["Aktywność metaboliczna (%)"].max() + 5
        for row in tukey.summary().data[1:]:
            grup1 = row[0]
            grup2 = row[1]
            p_adj = row[4]
            x1, x2 = grupy.tolist().index(grup1), grupy.tolist().index(grup2)
            stars = ""
            if p_adj < 0.001:
                stars = "*"
            elif p_adj < 0.01:
                stars = "**"
            elif p_adj < 0.05:
                stars = "*"
            if stars:
                ax.plot([x1, x1, x2, x2], [ymax, ymax + 1, ymax + 1, ymax], lw=1.5, c='k')
                ax.text((x1 + x2) / 2, ymax + 1.5, stars, ha='center', va='bottom', color='k')
                ymax += 6

    add_stat_annotation(ax, tukey)
    plt.xticks(rotation=45)
    plt.tight_layout()
    st.pyplot(plt.gcf())

    # Wykres 2: słupkowy z podziałem na czas i dawkę
    def rysuj_wykres_szczegolowy(dane_we, linia="MDA-MB-231"):
        st.write("## Wykres słupkowy z porównaniem grup (jak w publikacjach)")

        # Upewniamy się, że dane zawierają kolumnę 'Dawka (MBq/ml)' i 'Czas (h)'
        if "Dawka (MBq/ml)" not in dane_we.columns or "Czas (h)" not in dane_we.columns:
            st.error("Dane muszą zawierać kolumny 'Dawka (MBq/ml)' oraz 'Czas (h)'.")
            return

        kol_dawka = "Dawka (MBq/ml)"
        kol_viability = "Aktywność metaboliczna (%)"

        dane_linia = dane_we[dane_we["Linia"] == linia]
        dawki = sorted(dane_linia[kol_dawka].dropna().unique())
        czasy = sorted(dane_linia["Czas (h)"].dropna().unique())
        grupy = sorted(dane_linia["Grupa"].dropna().unique())

        fig, axs = plt.subplots(1, len(dawki), figsize=(6 * len(dawki), 6), sharey=True)

        for idx_d, dawka in enumerate(dawki):
            dane_dawka = dane_linia[dane_linia[kol_dawka] == dawka]
            ax = axs[idx_d] if len(dawki) > 1 else axs

            szerokosc_słupka = 0.15
            x = np.arange(len(czasy))

            for idx_g, grupa in enumerate(grupy):
                srodek = x + (idx_g - len(grupy) / 2) * szerokosc_słupka + szerokosc_słupka / 2
                wartosci = []
                bledy = []
                for czas in czasy:
                    podzbior = dane_dawka[
                        (dane_dawka["Grupa"] == grupa) & (dane_dawka["Czas (h)"] == czas)
                        ][kol_viability]
                    wartosci.append(podzbior.mean())
                    bledy.append(podzbior.sem())
                ax.bar(srodek, wartosci, yerr=bledy, width=szerokosc_słupka, label=grupa, capsize=3)

            ax.set_title(f"{chr(65 + idx_d)}. {dawka} MBq/ml")
            ax.set_xticks(x)
            ax.set_xticklabels([str(c) for c in czasy])
            ax.set_xlabel("Czas [h]")
            ax.set_ylabel("Aktywność metaboliczna [%]")
            ax.legend(title="Grupa")
            ax.grid(True, axis='y')

        plt.tight_layout()
        st.pyplot(fig)

    # Wywołanie wykresu szczegółowego
    # używamy tylko danych detailed
    dane_detailed = pd.concat([mda_detailed, t98g_detailed], ignore_index=True)
    rysuj_wykres_szczegolowy(dane_detailed, linia=linia_wybor)