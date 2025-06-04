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
from itertools import combinations

# Ścieżka lokalna
if getattr(sys, 'frozen', False):
    this_path = path.dirname(sys.executable)
elif __file__:
    this_path = path.abspath(path.dirname(__file__))

st.title("Wersja 5: Analiza cytotoksyczności: koniugat vs komponenty")

uploaded_file = st.file_uploader("Wgraj plik DANE.xlsx", type=["xlsx"])
if uploaded_file:
    # Wczytanie tylko danych detailed
    mda_detailed = pd.read_excel(uploaded_file, sheet_name="MDA detailed data")
    t98g_detailed = pd.read_excel(uploaded_file, sheet_name="T98G detailed data")

    # Normalizacja wartości dla porównań
    for df_set in [mda_detailed, t98g_detailed]:
        df_set["Grupa"] = df_set["Grupa"].str.strip()  # usuń spacje
        df_set["Czas (h)"] = df_set["Czas (h)"].astype(int)  # równe porównania
        df_set["Dawka (MBq/ml)"] = df_set["Dawka (MBq/ml)"].astype(float)
        df_set["Viability (%)"] = df_set["Viability (%)"].astype(float)
        df_set["Powtórzenie"] = df_set["Powtórzenie"].astype(int)

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

    def rysuj_wykres_szczegolowy2(dane_we, linia="MDA-MB-231"):
        st.write("## Wykres słupkowy z porównaniem grup (z testami istotności i kontrolą 0 MBq/ml)")

        kol_dawka = "Dawka (MBq/ml)"
        kol_viability = "Aktywność metaboliczna (%)"

        if kol_dawka not in dane_we.columns or "Czas (h)" not in dane_we.columns:
            st.error("Brakuje kolumn 'Dawka (MBq/ml)' lub 'Czas (h)' w danych.")
            return

        dane_linia = dane_we[dane_we["Linia"] == linia]
        czasy = sorted(dane_linia["Czas (h)"].dropna().unique())
        grupy = sorted(dane_linia["Grupa"].dropna().unique())

        dawki_bez_0 = sorted([d for d in dane_linia[kol_dawka].dropna().unique() if d != 0.0])

        fig, axs = plt.subplots(1, len(dawki_bez_0), figsize=(6 * len(dawki_bez_0), 6), sharey=True)

        if len(dawki_bez_0) == 1:
            axs = [axs]

        for idx_d, dawka in enumerate(dawki_bez_0):
            dane_dawka = dane_linia[(dane_linia[kol_dawka] == dawka) | (dane_linia[kol_dawka] == 0)]
            ax = axs[idx_d]

            szerokosc_słupka = 0.15
            x = np.arange(len(czasy))

            pozycje_słupków = {}
            wysokosci_słupków = {}
            for idx_g, grupa in enumerate(grupy):
                srodek = x + (idx_g - len(grupy) / 2) * szerokosc_słupka + szerokosc_słupka / 2
                pozycje_słupków[grupa] = srodek
                wartosci = []
                bledy = []
                for czas in czasy:
                    podzbior = dane_dawka[
                        (dane_dawka["Grupa"] == grupa) & (dane_dawka["Czas (h)"] == czas)
                        ][kol_viability]
                    wartosci.append(podzbior.mean())
                    bledy.append(podzbior.sem())
                ax.bar(srodek, wartosci, yerr=bledy, width=szerokosc_słupka, label=grupa, capsize=3)
                wysokosci_słupków[grupa] = wartosci

            # TESTY T: dawka 0 vs dawka X (dla tej samej grupy i czasu)
            offset = dane_dawka[kol_viability].max() * 0.05

            for czas_idx, czas in enumerate(czasy):
                for grupa in grupy:
                    kontrola = dane_linia[
                        (dane_linia[kol_dawka] == 0) &
                        (dane_linia["Czas (h)"] == czas) &
                        (dane_linia["Grupa"] == grupa)
                        ][kol_viability]

                    badana = dane_linia[
                        (dane_linia[kol_dawka] == dawka) &
                        (dane_linia["Czas (h)"] == czas) &
                        (dane_linia["Grupa"] == grupa)
                        ][kol_viability]

                    if len(kontrola) < 2 or len(badana) < 2:
                        continue

                    _, p = stats.ttest_ind(kontrola, badana, equal_var=False)

                    if p < 0.05:
                        poziom = "*" if p < 0.001 else "*" if p < 0.01 else ""
                        x_pos = pozycje_słupków[grupa][czas_idx]
                        y_pos = max(kontrola.mean(), badana.mean(), wysokosci_słupków[grupa][czas_idx]) + offset * 1.5
                        ax.text(x_pos, y_pos, poziom, ha='center', va='bottom', fontsize=14, color='black')

            ax.set_title(f"{chr(65 + idx_d)}. {dawka} MBq/ml (z kontrolą 0)")
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
    rysuj_wykres_szczegolowy2(dane_detailed, linia=linia_wybor)

    st.write("## Tabela statystyk opisowych")

    # Filtrowanie tylko dla wybranej linii
    dane_stat = dane_detailed[dane_detailed["Linia"] == linia_wybor].copy()

    # Normalizacja, żeby uniknąć błędów
    dane_stat["Grupa"] = dane_stat["Grupa"].str.strip()
    dane_stat["Czas (h)"] = dane_stat["Czas (h)"].astype(int)
    dane_stat["Dawka (MBq/ml)"] = dane_stat["Dawka (MBq/ml)"].astype(float)

    # Grupowanie i obliczanie statystyk
    tabela = (
        dane_stat.groupby(["Grupa", "Czas (h)", "Dawka (MBq/ml)"])["Aktywność metaboliczna (%)"]
        .agg([
            ("Średnia", "mean"),
            ("Min", "min"),
            ("Max", "max"),
            ("Odch. std", "std"),
            ("n", "count")
        ])
        .reset_index()
        .sort_values(by=["Grupa", "Czas (h)", "Dawka (MBq/ml)"])
    )

    # Wyświetlenie tabeli
    st.dataframe(tabela.style.format({
        "Średnia": "{:.2f}",
        "Min": "{:.2f}",
        "Max": "{:.2f}",
        "Odch. std": "{:.2f}",
        "n": "{:.0f}"
    }))

    # === Raport tekstowy z testów t dla każdej grupy ===
    st.write("## Raport statystyczny (t-test: porównanie każdej grupy z kontrolą 0 MBq/ml)")

    raport = ""

    for linia in dane_detailed["Linia"].unique():
        dane_linia = dane_detailed[dane_detailed["Linia"] == linia]
        dawki_testowe = sorted([d for d in dane_linia["Dawka (MBq/ml)"].dropna().unique() if d != 0.0])
        czasy = sorted(dane_linia["Czas (h)"].dropna().unique())
        grupy = sorted(dane_linia["Grupa"].dropna().unique())

        raport += f"\n### Linia komórkowa: {linia}\n"
        for dawka in dawki_testowe:
            for czas in czasy:
                raport += f"\n**Dawka: {dawka} MBq/ml | Czas: {czas} h**\n"
                for grupa in grupy:
                    kontrola = dane_linia[(dane_linia["Dawka (MBq/ml)"] == 0) &
                                          (dane_linia["Czas (h)"] == czas) &
                                          (dane_linia["Grupa"] == grupa)]["Aktywność metaboliczna (%)"]

                    badana = dane_linia[(dane_linia["Dawka (MBq/ml)"] == dawka) &
                                        (dane_linia["Czas (h)"] == czas) &
                                        (dane_linia["Grupa"] == grupa)]["Aktywność metaboliczna (%)"]

                    if len(kontrola) >= 2 and len(badana) >= 2:
                        _, p = stats.ttest_ind(kontrola, badana, equal_var=False)
                        stars = "*" if p < 0.001 else "*" if p < 0.01 else "" if p < 0.05 else ""
                        interpretacja = "istotna statystycznie" if p < 0.05 else "nieistotna"
                        raport += (
                            f"- Grupa *{grupa}*: "
                            f"średnia kontrola = {kontrola.mean():.2f}, "
                            f"średnia testowa = {badana.mean():.2f}, "
                            f"p = {p:.4f} → {interpretacja} {stars}\n"
                        )
                    else:
                        raport += f"- Grupa *{grupa}*: brak danych (za mało pomiarów)\n"

    # Wyświetl raport
    st.markdown(raport)


