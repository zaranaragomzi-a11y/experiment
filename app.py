# ==============================
# app.py — Streamlit Acid–Base Titration Lab (v2, fixed)
# ==============================
# 신규 기능
# 1) 지시약 추천 모드: 등가점 pH 기반 TOP 추천 & 적용 버튼
# 2) 실험 CSV 업로드: (volume_mL, pH) 오버레이 & RMSE 계산

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from math import log10
from scipy.optimize import brentq
import io
import pandas as pd

st.set_page_config(page_title="중화적정 가상 실험실", layout="wide")

# ---------- 상수 ----------
Kw = 1.0e-14  # 25 ℃
H2O_pKw = 14.0

# ---------- 유틸 ----------
def clamp(x, a, b):
    return max(a, min(b, x))

def mix_color(c1, c2, f):
    """두 색상 hex를 비율 f(0~1)로 선형 혼합"""
    c1 = c1.lstrip('#'); c2 = c2.lstrip('#')
    r1,g1,b1 = int(c1[0:2],16), int(c1[2:4],16), int(c1[4:6],16)
    r2,g2,b2 = int(c2[0:2],16), int(c2[2:4],16), int(c2[4:6],16)
    r = int(r1*(1-f)+r2*f); g = int(g1*(1-f)+g2*f); b = int(b1*(1-f)+b2*f)
    return f"#{r:02x}{g:02x}{b:02x}"

# ---------- 지시약 데이터 ----------
# name: (pH_low, pH_high, acid_hex, base_hex)
INDICATORS = {
    "메틸 오렌지": (3.1, 4.4, "#ff6b00", "#ffff00"),
    "메틸 레드":   (4.2, 6.3, "#ff0018", "#ffff00"),
    "브로모티몰 블루": (6.0, 7.6, "#ffff00", "#0000ff"),
    "페놀프탈레인": (8.2, 10.0, "#ffffff", "#ff00ff"),
    "팀올 블루(전이2)": (8.0, 9.6, "#00ffff", "#0000ff"),
}
IND_LIST = list(INDICATORS.keys())

# ---------- 화학 계산 유틸 ----------

def Ka_list_from_pKa(pKas):
    return [10**(-p) for p in pKas]

def amphiprotic_pH(pKa1, pKa2):
    return 0.5*(pKa1 + pKa2)

def alpha_polyacid(H, Kas):
    n = len(Kas)
    K = [1.0]
    for i in range(n):
        K.append(K[-1]*Kas[i])
    D = 0.0
    for i in range(n+1):
        D += K[i]*(H**(n-i))
    alphas = []
    for i in range(n+1):
        alphas.append(K[i]*(H**(n-i))/D)
    return alphas

def electroneutrality_root_H(total_acid_moles, Kas, total_base_moles, Kbs, V_L, Na_moles=0.0, Cl_moles=0.0):
    V = V_L
    Kw_val = Kw
    def f(logH):
        H = 10**(-logH)
        OH = Kw_val/H
        Na = Na_moles / V
        Cl = Cl_moles / V
        Ca = total_acid_moles / V if total_acid_moles > 0 else 0.0
        na = 0.0
        if total_acid_moles > 0 and len(Kas) > 0:
            al = alpha_polyacid(H, Kas)
            n = len(Kas)
            for i in range(1, n+1):
                na += i * (al[i]*Ca)
        Cb = total_base_moles / V if total_base_moles > 0 else 0.0
        pb = 0.0
        if total_base_moles > 0 and len(Kbs) > 0:
            Kb = Kbs[0]
            OHc = OH
            try:
                x = (-(OHc - Cb) + np.sqrt((OHc - Cb)**2 + 4*Kb*Cb)) / 2.0
                pb += x
            except Exception:
                pass
        lhs = H + Na + pb
        rhs = OH + Cl + na
        return lhs - rhs
    try:
        root = brentq(lambda x: f(x), -2, 16)
        return 10**(-root)
    except ValueError:
        grid = np.linspace(-2, 16, 361)
        vals = [f(x) for x in grid]
        idx = int(np.argmin(np.abs(vals)))
        return 10**(-grid[idx])

def compute_pH_curve(scn, Ca, Va_mL, Cb, Vb_max_mL, pKas_analyte, pKas_conj_acid_of_base,
                     analyte_type, titrant_type, valency_analyte, valency_titrant,
                     allow_weak_titrant=False):
    Va = Va_mL/1000.0
    Vb_list = np.arange(0.0, Vb_max_mL+1e-9, 1.0)
    n_analyte = Ca * Va
    Na_m = 0.0
    Cl_m = 0.0
    eq_points = []
    valA = valency_analyte
    valT = valency_titrant

    if analyte_type == 'acid':
        if titrant_type in ['strong-base','weak-base']:
            for k in range(1, valA+1):
                Veq = (n_analyte * k) / (Cb * valT) * 1000.0
                if 0 <= Veq <= Vb_max_mL*1.2:
                    eq_points.append(Veq)
    else:
        if titrant_type in ['strong-acid','weak-acid']:
            for k in range(1, valA+1):
                Veq = (n_analyte * k) / (Cb * valT) * 1000.0
                if 0 <= Veq <= Vb_max_mL*1.2:
                    eq_points.append(Veq)

    pH_list = []

    for Vb_mL in Vb_list:
        Vb = Vb_mL/1000.0
        Vtot = Va + Vb
        total_acid_moles = 0.0
        total_base_moles = 0.0
        Kas = []
        Kbs = []
        if analyte_type == 'acid':
            if len(pKas_analyte) > 0:
                Kas = Ka_list_from_pKa(pKas_analyte)
                total_acid_moles = n_analyte
            else:
                Cl_m += n_analyte*valA
        else:
            if len(pKas_conj_acid_of_base) > 0:
                Kas = Ka_list_from_pKa(pKas_conj_acid_of_base)
                total_acid_moles = n_analyte
            else:
                Na_m += n_analyte*valA
        if titrant_type == 'strong-base':
            Na_m += Cb * Vb * valT
        elif titrant_type == 'strong-acid':
            Cl_m += Cb * Vb * valT
        elif titrant_type == 'weak-base' and allow_weak_titrant:
            total_base_moles += Cb*Vb
            Kbs = [1.8e-5]
        elif titrant_type == 'weak-acid' and allow_weak_titrant:
            total_acid_moles += Cb*Vb
            Kas = Ka_list_from_pKa([4.76])
        H = electroneutrality_root_H(total_acid_moles, Kas, total_base_moles, Kbs, Vtot, Na_moles=Na_m, Cl_moles=Cl_m)
        pH_list.append(-np.log10(H))

    eq_pH = []
    for Veq in eq_points:
        Vb = Veq/1000.0
        Vtot = Va + Vb
        if analyte_type == 'acid':
            if len(pKas_analyte) == 0 and titrant_type == 'strong-base':
                eq_pH.append(7.00)
            elif len(pKas_analyte) == 1 and titrant_type == 'strong-base':
                Ka = 10**(-pKas_analyte[0])
                Kb = Kw/Ka
                C = (Ca*Va)/Vtot
                x = np.sqrt(max(Kb*C, 1e-30))
                pOH = -np.log10(x)
                eq_pH.append(H2O_pKw - pOH)
            elif len(pKas_analyte) >= 2 and titrant_type == 'strong-base':
                eq_pH.append(amphiprotic_pH(pKas_analyte[0], pKas_analyte[1]))
                if len(pKas_analyte) >= 3 and valency_analyte >= 3:
                    eq_pH.append(amphiprotic_pH(pKas_analyte[1], pKas_analyte[2]))
            else:
                H = electroneutrality_root_H((Ca*Va), Ka_list_from_pKa(pKas_analyte), 0.0, [], Vtot, Na_moles=0.0, Cl_moles=0.0)
                eq_pH.append(-np.log10(H))
        else:
            if len(pKas_conj_acid_of_base) == 0 and titrant_type == 'strong-acid':
                eq_pH.append(7.00)
            elif len(pKas_conj_acid_of_base) == 1 and titrant_type == 'strong-acid':
                Ka = 10**(-pKas_conj_acid_of_base[0])
                C = (Ca*Va)/Vtot
                x = np.sqrt(max(Ka*C, 1e-30))
                eq_pH.append(-np.log10(x))
            else:
                H = electroneutrality_root_H((Ca*Va), Ka_list_from_pKa(pKas_conj_acid_of_base), 0.0, [], Vtot)
                eq_pH.append(-np.log10(H))

    notes = "교육용 근사식을 포함합니다. 경계 영역(매우 묽은 용액, 특이 pKa 간격 등)에서는 완전한 수치해법이 더 정확합니다."
    return {"V": Vb_list, "pH": pH_list, "eq_points": eq_points, "eq_pH": eq_pH, "notes": notes}

# ... (중략: UI 코드 동일) ...

# 등가점/지시약 적합성 — ▶▶ 반드시 결과 계산 이후 위치
if len(eq_points) > 0:
    st.subheader("중화점(등가점) 정보")
    rows = []
    for i, Veq in enumerate(eq_points):
        e_pH = eq_pH[i] if i < len(eq_pH) else np.nan
        ok = (ind_low <= e_pH <= ind_high)
        verdict = "✅ 적합" if ok else "⚠️ 부적합"
        msg = f"중화점의 pH는 {e_pH:.2f}이므로 {ind_name}은(는) {'적합합니다.' if ok else '적합하지 않습니다.'}"
        rows.append((i+1, Veq, e_pH, verdict, msg))
    st.table({"EQ#":[r[0] for r in rows], "V_eq (mL)": [f"{r[1]:.1f}" for r in rows], "pH_eq": [f"{r[2]:.2f}" for r in rows], "판정": [r[3] for r in rows]})
    st.info("\n".join([r[4] for r in rows]))
else:
    st.info("등가점을 계산할 수 없는 설정입니다. 농도/부피/종류를 확인하세요.")
