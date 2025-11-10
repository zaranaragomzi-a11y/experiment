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

# 다프로톤 약산 HnA 분포함수 alpha

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

# 전하평형 해 찾기

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

# 시뮬레이터 본체

def compute_pH_curve(scn, Ca, Va_mL, Cb, Vb_max_mL, pKas_analyte, pKas_conj_acid_of_base,
                     analyte_type, titrant_type, valency_analyte, valency_titrant,
                     allow_weak_titrant=False):
    Va = Va_mL/1000.0
    Vb_list = np.arange(0.0, Vb_max_mL+1e-9, 1.0)  # 1 mL 간격
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

# ---------- UI ----------
st.title("🧪 중화적정 곡선 가상 실험실")
st.caption("GitHub + Streamlit 배포용 | 25℃ 가정 | 1 mL 간격 실시간 시뮬레이션")

colA, colB = st.columns([1,1])
with colA:
    st.subheader("분석 물질 (Analyte)")
    analyte_kind = st.selectbox("종류", ["강산", "약산(1가)", "약산(2가)", "약산(3가)", "강염기", "약염기(1가)", "약염기(2가)", "약염기(3가)"])
    Ca = st.number_input("분석 물질 농도 (M)", min_value=0.0001, max_value=5.0, value=0.100, step=0.001, format="%.3f")
    Va_mL = st.number_input("분석 물질 부피 (mL)", min_value=1.0, max_value=1000.0, value=25.0, step=1.0, format="%.1f")

    pKas_analyte = []
    pKas_conj_acid_of_base = []
    analyte_type = 'acid'
    valency_analyte = 1

    if analyte_kind == "강산":
        analyte_type = 'acid'; valency_analyte = 1; pKas_analyte = []
    elif analyte_kind == "약산(1가)":
        analyte_type = 'acid'; valency_analyte = 1
        pKas_analyte = [st.number_input("pKa(약산)", 0.0, 14.0, 4.76, 0.01)]
    elif analyte_kind == "약산(2가)":
        analyte_type = 'acid'; valency_analyte = 2
        p1 = st.number_input("pKa1", 0.0, 14.0, 2.15, 0.01)
        p2 = st.number_input("pKa2", 0.0, 14.0, 7.20, 0.01)
        pKas_analyte = [p1,p2]
    elif analyte_kind == "약산(3가)":
        analyte_type = 'acid'; valency_analyte = 3
        p1 = st.number_input("pKa1", 0.0, 14.0, 2.15, 0.01)
        p2 = st.number_input("pKa2", 0.0, 14.0, 7.20, 0.01)
        p3 = st.number_input("pKa3", 0.0, 14.0, 12.35, 0.01)
        pKas_analyte = [p1,p2,p3]
    elif analyte_kind == "강염기":
        analyte_type = 'base'; valency_analyte = 1; pKas_conj_acid_of_base = []
    elif analyte_kind == "약염기(1가)":
        analyte_type = 'base'; valency_analyte = 1
        pKb_or_pKa = st.number_input("pKb 또는 짝산 pKa(=14-pKb)", 0.0, 14.0, 4.75, 0.01)
        pKas_conj_acid_of_base = [pKb_or_pKa]
    elif analyte_kind == "약염기(2가)":
        analyte_type = 'base'; valency_analyte = 2
        pKa_conj1 = st.number_input("짝산 pKa1", 0.0, 14.0, 6.35, 0.01)
        pKa_conj2 = st.number_input("짝산 pKa2", 0.0, 14.0, 10.33, 0.01)
        pKas_conj_acid_of_base = [pKa_conj1, pKa_conj2]
    elif analyte_kind == "약염기(3가)":
        analyte_type = 'base'; valency_analyte = 3
        pKa_conj1 = st.number_input("짝산 pKa1", 0.0, 14.0, 2.15, 0.01)
        pKa_conj2 = st.number_input("짝산 pKa2", 0.0, 14.0, 7.20, 0.01)
        pKa_conj3 = st.number_input("짝산 pKa3", 0.0, 14.0, 12.35, 0.01)
        pKas_conj_acid_of_base = [pKa_conj1, pKa_conj2, pKa_conj3]

with colB:
    st.subheader("적정 용액 (Titrant)")
    if analyte_type == 'acid':
        titrant_choice = st.selectbox("종류", ["강염기", "약염기(단염기, 교육용)"])
        titrant_type = 'strong-base' if titrant_choice=="강염기" else 'weak-base'
        valency_titrant = 1
    else:
        titrant_choice = st.selectbox("종류", ["강산", "약산(단염기, 교육용)"])
        titrant_type = 'strong-acid' if titrant_choice=="강산" else 'weak-acid'
        valency_titrant = 1
    Cb = st.number_input("적정 용액 농도 (M)", min_value=0.0001, max_value=5.0, value=0.100, step=0.001, format="%.3f")
    Vb_max_mL = st.number_input("적정 용액 최대 주입량 (mL)", min_value=5.0, max_value=200.0, value=50.0, step=1.0, format="%.1f")

# 지시약 영역
st.subheader("지시약 (Indicator)")
ind_name = st.selectbox("선택", IND_LIST, index=IND_LIST.index("브로모티몰 블루") if "브로모티몰 블루" in IND_LIST else 0)
ind_low, ind_high, ind_acid_hex, ind_base_hex = INDICATORS[ind_name]
ind_mid = 0.5*(ind_low+ind_high)

# 계산 실행
res = compute_pH_curve(
    scn="main",
    Ca=Ca, Va_mL=Va_mL, Cb=Cb, Vb_max_mL=Vb_max_mL,
    pKas_analyte=pKas_analyte,
    pKas_conj_acid_of_base=pKas_conj_acid_of_base,
    analyte_type=analyte_type, titrant_type=titrant_type,
    valency_analyte=valency_analyte, valency_titrant=valency_titrant,
    allow_weak_titrant=True,
)
V_list = res["V"]; pH_list = res["pH"]; eq_points = res["eq_points"]; eq_pH = res["eq_pH"]

# 지시약 추천 모드
st.markdown("---")
with st.expander("🔎 지시약 추천 모드 (등가점 기반)", expanded=True):
    if len(eq_pH) == 0:
        st.info("등가점 pH를 계산할 수 없어 추천이 어렵습니다. 조건을 조정하세요.")
    else:
        targets = eq_pH
        rows = []
        for name, (lo, hi, *_colors) in INDICATORS.items():
            score_list = []
            cover_flags = []
            for p in targets:
                mid = 0.5*(lo+hi)
                if lo <= p <= hi:
                    score = 100 - abs(p - mid)*10
                    cover = True
                else:
                    score = 100 - (min(abs(p-lo), abs(p-hi))*20 + 50)
                    cover = False
                score_list.append(score)
                cover_flags.append(cover)
            score_avg = float(np.mean(score_list))
            coverage = all(cover_flags) if len(cover_flags)>0 else False
            rows.append((name, score_avg, coverage, (lo,hi)))
        rows.sort(key=lambda x: (x[2], x[1]), reverse=True)
        rec_names = [r[0] for r in rows[:3]]
        st.write("**추천 지시약 TOP 3:** ", ", ".join([f"{n}" for n in rec_names]))
        st.dataframe(pd.DataFrame({
            "지시약": [r[0] for r in rows],
            "점수(가중)": [f"{r[1]:.1f}" for r in rows],
            "전이구간": [f"{r[3][0]:.1f}–{r[3][1]:.1f}" for r in rows],
            "모든 등가점 커버": ["✅" if r[2] else "—" for r in rows],
        }))
        apply_name = st.selectbox("추천 적용", [r[0] for r in rows])
        if st.button("이 지시약 적용하기"):
            ind_name = apply_name
            ind_low, ind_high, ind_acid_hex, ind_base_hex = INDICATORS[ind_name]
            st.success(f"{ind_name} 적용 완료!")

# 인터랙티브 주입량
st.subheader("적정 진행")
V_now = st.slider("현재 적정 용액 주입량 (mL)", min_value=float(V_list[0]), max_value=float(V_list[-1]), value=float(V_list[0]), step=1.0)
idx_now = int((V_now - V_list[0]) / 1.0)
idx_now = clamp(idx_now, 0, len(pH_list)-1)
pH_now = pH_list[idx_now]

# 지시약 색상
if pH_now <= ind_low:
    sol_hex = ind_acid_hex
elif pH_now >= ind_high:
    sol_hex = ind_base_hex
else:
    f = (pH_now - ind_low) / (ind_high - ind_low)
    sol_hex = mix_color(ind_acid_hex, ind_base_hex, f)

# 그래프 + 색상 카드 + CSV 업로드/오버레이
col1, col2 = st.columns([2,1])
with col1:
    fig, ax = plt.subplots(figsize=(7,4))
    ax.plot(V_list, pH_list, linewidth=2, label="시뮬 곡선")
    ax.set_xlabel("적정 용액 주입량 (mL)")
    ax.set_ylabel("pH")
    ax.set_title("적정 곡선")
    ax.axvline(V_now, linestyle='--', linewidth=1)
    ax.scatter([V_now], [pH_now], label="현재", zorder=3)
    for i, Veq in enumerate(eq_points):
        ax.axvline(Veq, linestyle=':', linewidth=1)
        if i < len(eq_pH):
            ax.text(Veq, clamp(eq_pH[i],0,14)+0.2, f"EQ{i+1}", rotation=90, va='bottom', ha='right', fontsize=8)

    # ===== CSV 업로드 =====
    uploaded = st.file_uploader("실험 CSV 업로드 (열: volume_mL, pH)", type=["csv"])
    if uploaded is not None:
        try:
            df = pd.read_csv(uploaded)
        except Exception:
            uploaded.seek(0)
            data = uploaded.read().decode('utf-8')
            df = pd.read_csv(io.StringIO(data))
        cols = {c.lower().strip(): c for c in df.columns}
        vol_col = None; ph_col = None
        for key in cols:
            if 'vol' in key or 'ml' in key or '부피' in key:
                vol_col = cols[key]
            if 'ph' in key:
                ph_col = cols[key]
        if vol_col is None and len(df.columns)>=2:
            vol_col = df.columns[0]
        if ph_col is None and len(df.columns)>=2:
            ph_col = df.columns[1]
        expV = df[vol_col].astype(float).to_numpy()
        expP = df[ph_col].astype(float).to_numpy()
        ax.scatter(expV, expP, marker='x', label="실험 데이터", zorder=4)
        simP_interp = np.interp(expV, V_list, pH_list)
        valid = np.isfinite(simP_interp) & np.isfinite(expP)
        if valid.sum() > 0:
            rmse = float(np.sqrt(np.mean((simP_interp[valid]-expP[valid])**2)))
            st.toast(f"실험-시뮬 RMSE = {rmse:.3f} pH", icon="📐")
    ax.legend()
    st.pyplot(fig)

with col2:
    st.markdown("**용액 색상 (선택 지시약 기준)**")
    st.markdown(
        f"""
        <div style='width:100%;height:140px;border-radius:12px;border:1px solid #ccc;background:{sol_hex};'></div>
        <div style='margin-top:8px;'>현재 pH = <b>{pH_now:.2f}</b></div>
        """,
        unsafe_allow_html=True,
    )

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
    st.info("
".join([r[4] for r in rows]))
else:
    st.info("등가점을 계산할 수 없는 설정입니다. 농도/부피/종류를 확인하세요.")

st.markdown("---")
st.markdown("**계산 참고**")
st.caption(res["notes"])

st.markdown(
"""
**주요 가정**
- 25 ℃ 에서 \(K_w = 1.0	imes 10^{-14}\), 활동도계수 = 1.
- 약산/약염기: 완충구간 및 등가점에서 표준 근사(헨더슨–하셀발흐, 양쪽성 근사)를 사용하고, 전체 곡선은 전하평형 방정식을 브렌트 방법으로 수치해결.
- 다염기성 분석물질의 적정은 강산/강염기를 권장합니다. 단염기성 약산/약염기를 적정제로 선택하는 모드는 교육용 예시로만 제공됩니다.
"""
)

# ==============================
# README.md (요약)
# ==============================
README = r"""
# 중화적정 곡선 가상 실험실 (v2)

## 새 기능
- **지시약 추천 모드**: 등가점 pH와 전이구간을 비교해 가중 점수화 → 추천 TOP 3 & 적용 버튼
- **실험 CSV 오버레이**: `volume_mL, pH` 열을 가진 CSV를 업로드 → 시뮬 곡선과 비교, RMSE 계산

## 실행
```bash
pip install -r requirements.txt
streamlit run app.py
```

## CSV 포맷 예시
```csv
volume_mL,pH
0,2.90
1,2.95
...
```

## 수업 활용 팁
- 동일 조건에서 추천 지시약과 임의 지시약을 비교하며 **탐구 질문**을 유도하세요.
- 실험 데이터 업로드 후 RMSE를 줄이는 방향으로 **모수 역추정(농도, pKa)** 과제를 부여해 보세요.
"""

st.sidebar.markdown("**README (요약)**")
st.sidebar.code(README)

# ==============================
# requirements.txt (미니멀)
# ==============================
REQ = """
streamlit
numpy
scipy
matplotlib
pandas
"""

st.sidebar.markdown("**requirements.txt**")
st.sidebar.code(REQ)
